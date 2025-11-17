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
  const [shearData, setShearData] = useState({ labels: [], datasets: [] });
  const [momentData, setMomentData] = useState({ labels: [], datasets: [] });
  const [reactions, setReactions] = useState(null);
  const [scale, setScale] = useState(1.0);
  const [errorMsg, setErrorMsg] = useState(null);
  const [loading, setLoading] = useState(false);
  const SUPPORT_TYPES = ["PIN", "FIXED", "FREE", "ELASTIC", "PRESCRIBED"];
  const [leftExtra, setLeftExtra] = useState({});      // ky, ktheta, v, theta
  const [rightExtra, setRightExtra] = useState({});

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

    let constraints = [];

    // Helper to parse extra fields safely
    function parseExtra(obj, keys) {
      const result = {};
      for (const k of keys) {
        if (obj[k] !== undefined && obj[k] !== "") {
          const val = safeParseFloat(obj[k]);
          if (val !== null) result[k] = val;
        }
      }
      return Object.keys(result).length ? result : null;
    }

    // Left support
    if (["PIN", "FIXED", "FREE"].includes(leftType)) {
      constraints.push([0.0, leftType.toUpperCase()]);
    } else if (leftType === "ELASTIC") {
      const extra = parseExtra(leftExtra, ["ky", "ktheta"]);
      if (extra) constraints.push([0.0, "ELASTIC", extra]);
    } else if (leftType === "PRESCRIBED") {
      const extra = parseExtra(leftExtra, ["v", "theta"]);
      if (extra) constraints.push([0.0, "PRESCRIBED", extra]);
    }

    // Right support
    if (["PIN", "FIXED", "FREE"].includes(rightType)) {
      constraints.push([L, rightType.toUpperCase()]);
    } else if (rightType === "ELASTIC") {
      const extra = parseExtra(rightExtra, ["ky", "ktheta"]);
      if (extra) constraints.push([L, "ELASTIC", extra]);
    } else if (rightType === "PRESCRIBED") {
      const extra = parseExtra(rightExtra, ["v", "theta"]);
      if (extra) constraints.push([L, "PRESCRIBED", extra]);
    }

    // Interior supports
    interiorSupports.forEach(item => {
      let x = safeParseFloat(item.xStr, 0, L);
      if (x === null || x <= 0 || x >= L) return;

      if (["PIN", "FIXED", "FREE"].includes(item.type)) {
        constraints.push([x, item.type.toUpperCase()]);
      } else if (item.type === "ELASTIC") {
        const extra = parseExtra(item, ["ky", "ktheta"]);
        if (extra) constraints.push([x, "ELASTIC", extra]);
      } else if (item.type === "PRESCRIBED") {
        const extra = parseExtra(item, ["v", "theta"]);
        if (extra) constraints.push([x, "PRESCRIBED", extra]);
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
          loads.push([
            "dist",
            x0,
            x1,
            (x) => safeQ(x, expr, L)
          ]);
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

      // Deflection / undeformed
      setChartData({
        labels: res.x.map(xi => xi.toFixed(3)),
        datasets: [
          { label: "Deformed (scaled)", data: res.v_scaled, borderColor: "rgb(75,192,192)", tension: 0.3 },
          { label: "Undeformed", data: res.undeformed, borderColor: "#000", borderDash: [5, 5] },
        ],
      });

      // Shear and Moment: use the returned sampled arrays (x_samples, shear_samples, moment_samples)
      if (res.x_samples && res.shear_samples && res.moment_samples) {
        const labels = res.x_samples.map(xi => xi.toFixed(3));
        setShearData({
          labels,
          datasets: [{ label: "Shear V(x) (N)", data: res.shear_samples, tension: 0.1 }],
        });
        setMomentData({
          labels,
          datasets: [{ label: "Moment M(x) (Nm)", data: res.moment_samples, tension: 0.1 }],
        });
      } else {
        setShearData({ labels: [], datasets: [] });
        setMomentData({ labels: [], datasets: [] });
      }
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
  }, [
    lengthStr,
    EStr,
    IStr,
    nElementsStr,
    leftType,
    rightType,
    JSON.stringify(leftExtra),     
    JSON.stringify(rightExtra),     
    JSON.stringify(interiorSupports), // already correct, but now includes ky/kθ/v/θ
    JSON.stringify(loadItems)
  ]);

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
    setInteriorSupports(prev => [...prev, {
      id: idCounter.current++,
      xStr: mid.toFixed(1),
      type: "PIN",
      ky: "", ktheta: "", v: "", theta: ""   // extra fields
    }]);
  };
  const updateInterior = (id, field, value) => setInteriorSupports(prev =>
    prev.map(i => i.id === id ? { ...i, [field]: value } : i)
  );
  const deleteInterior = (id) => setInteriorSupports(prev => prev.filter(i => i.id !== id));

  return (
    <div className="p-6 grid grid-cols-1 lg:grid-cols-2 gap-6">
      {/* LEFT PANEL: Inputs */}
      <Card>
        <h2 className="text-2xl font-bold mb-4">Beam Parameters & Supports</h2>
        <div className="grid grid-cols-2 gap-4">
          <div>
            <Label>Length (m)</Label>
            <Input value={lengthStr} onChange={(e) => setLengthStr(e.target.value)} placeholder="6.0" />
          </div>
          <div>
            <Label>E (Pa)</Label>
            <Input value={EStr} onChange={(e) => setEStr(e.target.value)} placeholder="210e9" />
          </div>
          <div>
            <Label>I (m⁴)</Label>
            <Input value={IStr} onChange={(e) => setIStr(e.target.value)} placeholder="8.1e-6" />
          </div>
          <div>
            <Label>Mesh elements</Label>
            <Input value={nElementsStr} onChange={(e) => setNElementsStr(e.target.value)} placeholder="200" />
          </div>
          <div>
            <Label>Left support (x=0)</Label>
            <select className="border p-2 w-full" value={leftType} onChange={e => setLeftType(e.target.value)}>
              {SUPPORT_TYPES.map(t => <option key={t}>{t}</option>)}
            </select>
            {["ELASTIC", "PRESCRIBED"].includes(leftType) && (
              <div className="mt-2 pl-4 space-y-2">
                {leftType === "ELASTIC" && (
                  <>
                    <Input placeholder="ky (N/m)" value={leftExtra?.ky || ""} onChange={e => setLeftExtra({ ...leftExtra, ky: e.target.value })} />
                    <Input placeholder="kθ (Nm/rad)" value={leftExtra?.ktheta || ""} onChange={e => setLeftExtra({ ...leftExtra, ktheta: e.target.value })} />
                  </>
                )}
                {leftType === "PRESCRIBED" && (
                  <>
                    <Input placeholder="Prescribed v (m)" value={leftExtra?.v || ""} onChange={e => setLeftExtra({ ...leftExtra, v: e.target.value })} />
                    <Input placeholder="Prescribed θ (rad)" value={leftExtra?.theta || ""} onChange={e => setLeftExtra({ ...leftExtra, theta: e.target.value })} />
                  </>
                )}
              </div>
            )}
          </div>

          <div>
            <Label>Right support (x={lengthStr || "L"})</Label>
            <select className="border p-2 w-full" value={rightType} onChange={e => setRightType(e.target.value)}>
              {SUPPORT_TYPES.map(t => <option key={t}>{t}</option>)}
            </select>
            {["ELASTIC", "PRESCRIBED"].includes(rightType) && (
              <div className="mt-2 pl-4 space-y-2">
                {rightType === "ELASTIC" && (
                  <>
                    <Input placeholder="ky (N/m)" value={rightExtra?.ky || ""} onChange={e => setRightExtra({ ...rightExtra, ky: e.target.value })} />
                    <Input placeholder="kθ (Nm/rad)" value={rightExtra?.ktheta || ""} onChange={e => setRightExtra({ ...rightExtra, ktheta: e.target.value })} />
                  </>
                )}
                {rightType === "PRESCRIBED" && (
                  <>
                    <Input placeholder="Prescribed v (m)" value={rightExtra?.v || ""} onChange={e => setRightExtra({ ...rightExtra, v: e.target.value })} />
                    <Input placeholder="Prescribed θ (rad)" value={rightExtra?.theta || ""} onChange={e => setRightExtra({ ...rightExtra, theta: e.target.value })} />
                  </>
                )}
              </div>
            )}
          </div>
        </div>

        <h3 className="mt-6 font-semibold">Interior Supports</h3>
        {interiorSupports.map(item => (
          <div key={item.id} className="border rounded p-4 mt-4 bg-gray-50">
            <div className="flex gap-2 items-end">
              <div className="flex-1">
                <Label>Position x (m)</Label>
                <Input value={item.xStr} onChange={e => updateInterior(item.id, 'xStr', e.target.value)} />
              </div>
              <select value={item.type} onChange={e => updateInterior(item.id, 'type', e.target.value)}>
                {SUPPORT_TYPES.map(t => <option key={t}>{t}</option>)}
              </select>
              <Button onClick={() => deleteInterior(item.id)}>Delete</Button>
            </div>

            {["ELASTIC", "PRESCRIBED"].includes(item.type) && (
              <div className="mt-3 pl-4 grid grid-cols-2 gap-3">
                {item.type === "ELASTIC" && (
                  <>
                    <div><Label>ky (N/m)</Label><Input value={item.ky || ""} onChange={e => updateInterior(item.id, 'ky', e.target.value)} /></div>
                    <div><Label>kθ (Nm/rad)</Label><Input value={item.ktheta || ""} onChange={e => updateInterior(item.id, 'ktheta', e.target.value)} /></div>
                  </>
                )}
                {item.type === "PRESCRIBED" && (
                  <>
                    <div><Label>Prescribed v (m)</Label><Input value={item.v || ""} onChange={e => updateInterior(item.id, 'v', e.target.value)} /></div>
                    <div><Label>Prescribed θ (rad)</Label><Input value={item.theta || ""} onChange={e => updateInterior(item.id, 'theta', e.target.value)} /></div>
                  </>
                )}
              </div>
            )}
          </div>
        ))}
        <Button className="mt-2 w-full" onClick={addInteriorSupport}>Add Interior Support</Button>

        <h3 className="mt-6 font-semibold">Loads</h3>
        {loadItems.map(item => (
          <div key={item.id} className="border rounded p-4 mt-4">
            <select value={item.type} onChange={e => updateLoad(item.id, 'type', e.target.value)}>
              <option value="point">Point Load (+ up)</option>
              <option value="moment">Moment (+ CCW)</option>
              <option value="distributed_load">Distributed Load</option>
            </select>
            {(item.type === "point" || item.type === "moment") && (
              <>
                <Label>Position (m)</Label>
                <Input value={item.xStr} onChange={e => updateLoad(item.id, 'xStr', e.target.value)} />
                <Label>Value {item.type === "point" ? "(N)" : "(Nm)"}</Label>
                <Input value={item.valStr} onChange={e => updateLoad(item.id, 'valStr', e.target.value)} />
              </>
            )}
            {item.type === "distributed_load" && (
              <>
                <Label>Start x (m)</Label><Input value={item.x0Str} onChange={e => updateLoad(item.id, 'x0Str', e.target.value)} />
                <Label>End x (m)</Label><Input value={item.x1Str} onChange={e => updateLoad(item.id, 'x1Str', e.target.value)} />
                <label className="flex items-center gap-2">
                  <input type="checkbox" checked={item.isVariable || false} onChange={e => updateLoad(item.id, 'isVariable', e.target.checked)} />
                  Variable q(x)
                </label>
                {item.isVariable ? (
                  <>
                    <Label>q(x) expression</Label>
                    <Input value={item.exprStr} onChange={e => updateLoad(item.id, 'exprStr', e.target.value)} />
                    <small>Example: 5000*(L-x)</small>
                  </>
                ) : (
                  <>
                    <Label>Constant q (N/m)</Label>
                    <Input value={item.qStr} onChange={e => updateLoad(item.id, 'qStr', e.target.value)} />
                  </>
                )}
              </>
            )}
            <Button className="mt-2 w-full" onClick={() => deleteLoad(item.id)}>Delete Load</Button>
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
          <>
            <div className="mb-4">
              <Line data={chartData} options={{ responsive: true, plugins: { legend: { position: "top" } } }} />
            </div>

            {/* Shear diagram */}
            <div className="mb-4">
              <h3 className="font-semibold">Shear Diagram</h3>
              {shearData.labels?.length ? (
                <Line data={shearData} options={{ responsive: true, plugins: { legend: { position: "top" } }, scales: { y: { beginAtZero: false } } }} />
              ) : <p className="text-gray-500 italic">No shear data available.</p>}
            </div>

            {/* Moment diagram */}
            <div className="mb-4">
              <h3 className="font-semibold">Moment Diagram</h3>
              {momentData.labels?.length ? (
                <Line data={momentData} options={{ responsive: true, plugins: { legend: { position: "top" } }, scales: { y: { beginAtZero: false } } }} />
              ) : <p className="text-gray-500 italic">No moment data available.</p>}
            </div>
          </>
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
