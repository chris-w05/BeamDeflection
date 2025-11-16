// frontend/src/App.jsx
import React, { useEffect, useState, useRef } from "react";
import { runBeamSimulation } from "./api";
import { Line } from "react-chartjs-2";
import {
  Chart as ChartJS,
  CategoryScale,
  LinearScale,
  PointElement,
  LineElement,
  Title,
  Tooltip,
  Legend,
} from "chart.js";

import { evaluate, parse } from "mathjs";



ChartJS.register(CategoryScale, LinearScale, PointElement, LineElement, Title, Tooltip, Legend);

const Button = ({ children, ...p }) => (
  <button className="px-3 py-1 bg-blue-500 text-white rounded hover:bg-blue-600" {...p}>
    {children}
  </button>
);
const Card = ({ children, ...p }) => <div className="border rounded-lg shadow p-6 bg-white" {...p}>{children}</div>;
const Input = (props) => <input className="border p-2 rounded w-full" {...props} />;
const Label = (props) => <label className="block font-semibold mb-1" {...props} />;

export default function App() {
  const [lengthStr, setLengthStr] = useState("");
  const [EStr, setEStr] = useState("");
  const [IStr, setIStr] = useState("");
  const [nElementsStr, setNElementsStr] = useState("");
  const [leftType, setLeftType] = useState("PIN");
  const [rightType, setRightType] = useState("PIN");
  const [interiorSupports, setInteriorSupports] = useState([]);
  const [loadItems, setLoadItems] = useState([]);
  const [chartData, setChartData] = useState({ labels: [], datasets: [] });
  const [reactions, setReactions] = useState(null);
  const [scale, setScale] = useState(1.0);
  const [errorMsg, setErrorMsg] = useState(null);
  const [loading, setLoading] = useState(false);

  const debounceRef = useRef(null);
  const idCounter = useRef(0);

  // Default: mid-span moment
  useEffect(() => {
    if (loadItems.length === 0) addLoad("moment");
  }, []);

  function safeParseFloat(s, min = -Infinity, max = Infinity) {
    if (!s?.trim()) return null;
    const v = parseFloat(s);
    if (!isFinite(v)) return null;
    return Math.max(min, Math.min(max, v)); // clamp to [min, max]
  }

  function clampValue(value, min = -Infinity, max = Infinity) {
    let v = parseFloat(value);
    if (!isFinite(v)) return null;
    return Math.max(min, Math.min(max, v));
  }


  //For parsing distributed loads
  function safeQ(x, expr, L) {
    try {
      // Create a mathjs scope with variables x and L, and math functions
      const scope = { x, L, e: Math.E, pi: Math.PI };
      return evaluate(expr, scope);
    } catch (err) {
      console.error("Error evaluating q(x):", expr, err);
      return 0; // fallback
    }
  }


  function buildPayload() {
    const L = safeParseFloat(lengthStr, 0.001) || 6.0;       // length > 0
    const E = safeParseFloat(EStr, 1e-6) || 210e9;           // E > 0
    const I = safeParseFloat(IStr, 1e-12) || 8.1e-6;         // I > 0
    const nElements = safeParseFloat(nElementsStr, 1, 1000) || 200; // 1 <= nElements <= 1000

    if ([L, E, I, nElements].some(x => x === null)) {
      setErrorMsg("Invalid input: all fields must be numeric and within allowed ranges.");
      return null;
    }

    let constraints = [
      [0.0, leftType.toUpperCase()],
      [L, rightType.toUpperCase()],
    ];

    interiorSupports.forEach(item => {
      let x = safeParseFloat(item.xStr, 0, L);
      if (x !== null && x > 0 && x < L) {
        constraints.push([x, item.type.toUpperCase()]);
      }
    });

    constraints.sort((a, b) => a[0] - b[0]);

    const loads = [];
    loadItems.forEach(item => {
      if (item.type === "point") {
        let x = safeParseFloat(item.xStr, 0, L);
        const Fy = safeParseFloat(item.valStr);
        if (x !== null && Fy !== null) loads.push(["point", x, Fy]);
      } else if (item.type === "moment") {
        let x = safeParseFloat(item.xStr, 0, L);
        const M = safeParseFloat(item.valStr);
        if (x !== null && M !== null) loads.push(["moment", x, M]);
      } else if (item.type === "distributed_load") {
        const x0 = safeParseFloat(item.x0Str, 0, L);
        const x1 = safeParseFloat(item.x1Str, 0, L);
        if (x0 === null || x1 === null || x1 <= x0) return;

        if (item.isVariable) {
          const expr = item.exprStr?.trim();
          if (!expr) return;
          // q(x) will be processed safely later
          loads.push(["dist", x0, x1, expr]);
        } else {
          const q = safeParseFloat(item.qStr);
          if (q === null) return;
          loads.push(["dist", x0, x1, q]);
        }
      }
    });

    return { length: L, E, I, nElements, constraints, loads };
  }

  async function runSim() {
    const payload = buildPayload();
    if (!payload) return;
    setLoading(true);
    setErrorMsg(null);
    try {
      const res = await runBeamSimulation(payload);
      setReactions(res.reactions);
      setScale(res.scale);
      setChartData({
        labels: res.x.map(xi => xi.toFixed(3)),
        datasets: [
          { label: "Deformed (scaled)", data: res.v_scaled, borderColor: "rgb(75,192,192)", tension: 0.3 },
          { label: "Undeformed", data: res.undeformed, borderColor: "#000", borderDash: [5, 5] },
        ],
      });
    } catch (err) {
      setErrorMsg(err.message || "API error");
    } finally {
      setLoading(false);
    }
  }

  useEffect(() => {
    if (debounceRef.current) clearTimeout(debounceRef.current);
    debounceRef.current = setTimeout(runSim, 500);
    return () => clearTimeout(debounceRef.current);
  }, [lengthStr, EStr, IStr, nElementsStr, leftType, rightType, JSON.stringify(interiorSupports), JSON.stringify(loadItems)]);

  const addLoad = (type) => {
    const currentL = safeParseFloat(lengthStr);
    const mid = currentL ? (currentL / 2).toFixed(2) : "";
    const end = currentL ? currentL.toFixed(2) : "";
    setLoadItems(prev => [...prev, {
      id: idCounter.current++,
      type,
      xStr: type === "distributed_load" ? "" : mid,
      x0Str: "0.0",
      x1Str: end,
      qStr: "-5000",
      exprStr: "-5000 + 10000 * (x / L)",
      isVariable: false,
      valStr: type === "point" ? "-10000" : (type === "moment" ? "10000" : ""),
    }]);
  };

  const updateLoad = (id, field, value) => setLoadItems(prev => prev.map(i => i.id === id ? { ...i, [field]: value } : i));
  const deleteLoad = (id) => setLoadItems(prev => prev.filter(i => i.id !== id));
  const addInteriorSupport = () => {
    const mid = (safeParseFloat(lengthStr) || 6) / 2;
    setInteriorSupports(prev => [...prev, { id: idCounter.current++, xStr: mid.toFixed(1), type: "PIN" }]);
  };
  const updateInterior = (id, field, value) => setInteriorSupports(prev => prev.map(i => i.id === id ? { ...i, [field]: value } : i));
  const deleteInterior = (id) => setInteriorSupports(prev => prev.filter(i => i.id !== id));

  return (
    <div className="p-6 grid grid-cols-1 lg:grid-cols-2 gap-6">
      {/* LEFT PANEL: Inputs */}
      <Card>
        <h2 className="text-2xl font-bold mb-4">Beam Parameters & Supports</h2>
        <div className="grid grid-cols-2 gap-4">
          <div>
            <Label>Length (m)</Label>
            <Input
              value={lengthStr}
              onChange={(e) => {
                let val = parseFloat(e.target.value);
                if (isNaN(val) || val <= 0) val = 0.1; // minimum length
                setLengthStr(val.toString());
              }}
              placeholder="6.0"
            />
          </div>
          <div>
            <Label>E (Pa)</Label>
            <Input
              value={EStr}
              onChange={(e) => {
                let val = parseFloat(e.target.value);
                if (isNaN(val) || val <= 0) val = 1e3; // minimum E
                setEStr(val.toString());
              }}
              placeholder="210e9"
            />
          </div>
          <div>
            <Label>I (m⁴)</Label>
            <Input
              value={IStr}
              onChange={(e) => {
                let val = parseFloat(e.target.value);
                if (isNaN(val) || val <= 0) val = 1e-12; // minimum I
                setIStr(val.toString());
              }}
              placeholder="8.1e-6"
            />
          </div>
          <div>
            <Label>Mesh elements</Label>
            <Input
              value={nElementsStr}
              onChange={(e) => {
                let val = parseInt(e.target.value);
                if (isNaN(val) || val < 2) val = 2; // minimum 2 elements
                if (val > 1000) val = 1000; // maximum 1000 elements
                setNElementsStr(val.toString());
              }}
              placeholder="200"
            />
          </div>
          <div>
            <Label>Left support (x=0)</Label>
            <select
              className="border p-2 w-full"
              value={leftType}
              onChange={(e) => setLeftType(e.target.value)}
            >
              <option>PIN</option>
              <option>FIXED</option>
              <option>FREE</option>
            </select>
          </div>
          <div>
            <Label>Right support (x={lengthStr || "L"})</Label>
            <select
              className="border p-2 w-full"
              value={rightType}
              onChange={(e) => setRightType(e.target.value)}
            >
              <option>PIN</option>
              <option>FIXED</option>
              <option>FREE</option>
            </select>
          </div>
        </div>

        {/* Interior Supports */}
        <h3 className="mt-6 font-semibold">Interior Supports</h3>
        {interiorSupports.map((item) => (
          <div key={item.id} className="flex gap-2 mt-2">
            <Input
              placeholder="x (m)"
              value={item.xStr}
              onChange={(e) => {
                let val = parseFloat(e.target.value);
                const L = parseFloat(lengthStr) || 6.0;
                if (isNaN(val) || val < 0) val = 0;
                if (val > L) val = L;
                updateInterior(item.id, "xStr", val.toString());
              }}
            />
            <select
              value={item.type}
              onChange={(e) => updateInterior(item.id, "type", e.target.value)}
            >
              <option>PIN</option>
              <option>FIXED</option>
              <option>FREE</option>
            </select>
            <Button onClick={() => deleteInterior(item.id)}>Delete</Button>
          </div>
        ))}
        <Button className="mt-2 w-full" onClick={addInteriorSupport}>
          Add Interior Support
        </Button>

        {/* Loads */}
        <h3 className="mt-6 font-semibold">Loads</h3>
        {loadItems.map((item) => (
          <div key={item.id} className="border rounded p-4 mt-4">
            <select
              value={item.type}
              onChange={(e) => updateLoad(item.id, "type", e.target.value)}
            >
              <option value="point">Point Load (+ up)</option>
              <option value="moment">Moment (+ CCW)</option>
              <option value="distributed_load">Distributed Load</option>
            </select>
            {(item.type === "point" || item.type === "moment") && (
              <>
                <Label>Position (m)</Label>
                <Input
                  value={item.xStr}
                  onChange={(e) => {
                    let val = parseFloat(e.target.value);
                    const L = parseFloat(lengthStr) || 6.0;
                    if (isNaN(val) || val < 0) val = 0;
                    if (val > L) val = L;
                    updateLoad(item.id, "xStr", val.toString());
                  }}
                />
                <Label>Value {item.type === "point" ? "(N)" : "(Nm)"}</Label>
                <Input
                  value={item.valStr}
                  onChange={(e) => {
                    let val = parseFloat(e.target.value);
                    if (isNaN(val)) val = 0;
                    updateLoad(item.id, "valStr", val.toString());
                  }}
                />
              </>
            )}
            {item.type === "distributed_load" && (
              <>
                <Label>Start x (m)</Label>
                <Input
                  value={item.x0Str}
                  onChange={(e) => {
                    let val = parseFloat(e.target.value);
                    const L = parseFloat(lengthStr) || 6.0;
                    if (isNaN(val) || val < 0) val = 0;
                    if (val > L) val = L;
                    updateLoad(item.id, "x0Str", val.toString());
                  }}
                />
                <Label>End x (m)</Label>
                <Input
                  value={item.x1Str}
                  onChange={(e) => {
                    let val = parseFloat(e.target.value);
                    const L = parseFloat(lengthStr) || 6.0;
                    if (isNaN(val) || val < 0) val = 0;
                    if (val > L) val = L;
                    updateLoad(item.id, "x1Str", val.toString());
                  }}
                />
                <label className="flex items-center gap-2">
                  <input
                    type="checkbox"
                    checked={item.isVariable || false}
                    onChange={(e) => updateLoad(item.id, "isVariable", e.target.checked)}
                  />
                  Variable q(x)
                </label>
                {item.isVariable ? (
                  <>
                    <Label>q(x) expression</Label>
                    <Input
                      value={item.exprStr}
                      onChange={(e) => updateLoad(item.id, "exprStr", e.target.value)}
                    />
                    <small>Example: 5000*(L-x)</small>
                  </>
                ) : (
                  <>
                    <Label>Constant q (N/m)</Label>
                    <Input
                      value={item.qStr}
                      onChange={(e) => {
                        let val = parseFloat(e.target.value);
                        if (isNaN(val)) val = 0;
                        updateLoad(item.id, "qStr", val.toString());
                      }}
                    />
                  </>
                )}
              </>
            )}
            <Button className="mt-2 w-full" onClick={() => deleteLoad(item.id)}>
              Delete Load
            </Button>
          </div>
        ))}
        <div className="mt-4 flex flex-col gap-2">
          <Button onClick={() => addLoad("point")}>Add Point Load</Button>
          <Button onClick={() => addLoad("moment")}>Add Moment</Button>
          <Button onClick={() => addLoad("distributed_load")}>Add Distributed Load</Button>
        </div>
        <div className="mt-6 text-sm text-gray-600">
          Auto-solving on change (500ms debounce). Empty fields allowed while typing.
          {loading && <p className="text-blue-600">Solving...</p>}
          {errorMsg && <p className="text-red-600">Error: {errorMsg}</p>}
        </div>
      </Card>

      {/* RIGHT PANEL: Graph */}
      <Card>
        <h2 className="text-2xl font-bold mb-4">Beam Deformation {scale !== 1.0 && `(scale: ${scale})`}</h2>
        {loading && <p className="text-blue-600">Solving...</p>}
        {errorMsg && <p className="text-red-600">{errorMsg}</p>}
        {!loading && chartData.labels?.length > 0 ? (
          <Line data={chartData} options={{ responsive: true, plugins: { legend: { position: "top" } } }} />
        ) : (
          <p className="text-gray-500 italic">
            Enter beam length, E, I, add supports/loads → deflection plot will appear automatically.
          </p>
        )}
        {reactions && chartData.labels?.length > 0 && (
          <div className="mt-4 text-sm">
            <b>Reactions (first 8):</b> {reactions.slice(0, 8).map(r => typeof r === 'number' ? r.toFixed(3) : r).join(", ")}
          </div>
        )}
      </Card>
    </div>
  );
}
