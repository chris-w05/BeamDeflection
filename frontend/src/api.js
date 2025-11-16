import { solveBeam } from "./beamSolver";

export async function runBeamSimulation(params) {
    // Convert distributed load expressions into JS functions
    const loads = params.loads.map(ld => {
        if (ld[0] === "dist" && typeof ld[3] === "string") {
            const expr = ld[3];
            const L = params.length;
            const f = new Function("x", "L", `return ${expr};`);
            return ["dist", ld[1], ld[2], (x) => f(x, L)];
        }
        return ld;
    });

    const result = solveBeam(params.length, params.nElements, params.E, params.I, params.constraints, loads);
    return result;
}
