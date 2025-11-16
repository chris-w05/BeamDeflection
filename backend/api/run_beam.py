# backend/api/run_beam.py (full updated file - replace entirely)
from fastapi import APIRouter, HTTPException, Request
from fastapi.responses import JSONResponse
import numpy as np
from .BeamDeflection import solve_beam

router = APIRouter()

def sanitize_number(x, name):
    try:
        if x is None:
            raise ValueError("Missing")
        v = float(x)
        if not np.isfinite(v):
            raise ValueError("Not finite")
        return v
    except Exception:
        raise HTTPException(status_code=400, detail=f"Invalid numeric value for {name}: {x}")

@router.post("/run_beam")
async def run_beam(request: Request):
    data = await request.json()

    L = sanitize_number(data.get("length"), "length")
    E = sanitize_number(data.get("E"), "E")
    I = sanitize_number(data.get("I"), "I")
    n_elements = int(sanitize_number(data.get("nElements", 200), "nElements"))  # fixed typo

    # Constraints
    raw_constraints = data.get("constraints", [])
    constraints = []
    for c in raw_constraints:
        if not isinstance(c, (list, tuple)) or len(c) < 2:
            continue
        try:
            xpos = float(c[0])
            xpos = max(0.0, min(L, xpos))
            ctype = str(c[1]).upper()
            constraints.append((xpos, ctype))
        except:
            raise HTTPException(400, detail=f"Invalid constraint: {c}")

    # Loads - flexible length check
    raw_loads = data.get("loads", [])
    loads = []
    for ld in raw_loads:
        if not isinstance(ld, (list, tuple)):
            raise HTTPException(400, detail=f"Load not list/tuple: {ld}")
        if len(ld) < 3:
            raise HTTPException(400, detail=f"Load tuple too short (min 3 items): {ld}")
        
        typ = str(ld[0]).lower()
        if typ in ("point", "point_load"):
            if len(ld) != 3:
                raise HTTPException(400, detail=f"Point load expects exactly 3 items [type, x, Fy]: {ld}")
            x = sanitize_number(ld[1], "point x")
            Fy = sanitize_number(ld[2], "point Fy")
            x = max(0.0, min(L, x))
            loads.append(("point", x, Fy))
        
        elif typ in ("moment", "moment_load"):
            if len(ld) != 3:
                raise HTTPException(400, detail=f"Moment load expects exactly 3 items [type, x, M]: {ld}")
            x = sanitize_number(ld[1], "moment x")
            M = sanitize_number(ld[2], "moment M")
            x = max(0.0, min(L, x))
            loads.append(("moment", x, M))
        
        elif typ in ("dist", "distributed", "distributed_load"):
            if len(ld) != 4:
                raise HTTPException(400, detail=f"Distributed load expects exactly 4 items [type, x0, x1, q]: {ld}")
            x0 = sanitize_number(ld[1], "dist x0")
            x1 = sanitize_number(ld[2], "dist x1")
            raw_q = ld[3]
            # order-independent
            start = min(x0, x1)
            end = max(x0, x1)
            start = max(0.0, start)
            end = min(L, end)
            if end <= start:
                continue  # skip zero/negative length
            
            if isinstance(raw_q, str):
                expr = raw_q.strip()
                if not expr:
                    raise HTTPException(400, detail="Empty q(x) string")
                # Very safe eval - np + bare math
                q = lambda xx: eval(expr, {"__builtins__": {}}, {
                    "x": xx, "L": L, "np": np,
                    "sin": np.sin, "cos": np.cos, "tan": np.tan,
                    "pi": np.pi, "sqrt": np.sqrt, "exp": np.exp, "log": np.log
                })
                try:
                    _ = q(0.0)  # test call
                except Exception as e:
                    raise HTTPException(400, detail=f"Invalid q(x)='{expr}': {e}")
            else:
                q = sanitize_number(raw_q, "dist q")
            loads.append(("dist", start, end, q))
        
        else:
            raise HTTPException(400, detail=f"Unknown load type '{typ}' in {ld}")

    # Solve
    try:
        xs, u, reactions = solve_beam(L, n_elements, E, I, constraints, loads)
        vs = u[0::2]
        max_def = np.max(np.abs(vs))
        scale = 0.2 * L / max_def if max_def > 1e-8 else 1.0
        return JSONResponse({
            "x": xs.tolist(),
            "v_scaled": (vs).tolist(),
            "undeformed": [0.0] * len(xs),
            "reactions": reactions.tolist(),
            "scale": 1 #round(scale, 3)
        })
    except np.linalg.LinAlgError:
        raise HTTPException(400, detail="Beam unstable - add more supports (need at least two PIN/FIXED or equivalent)")
    except Exception as e:
        raise HTTPException(500, detail=f"Solver crash: {e}")