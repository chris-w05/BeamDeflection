// frontend/src/beamSolver.js

// Utilities
function linspace(start, end, num) {
    const arr = [];
    if (num === 1) return [start];
    const step = (end - start) / (num - 1);
    for (let i = 0; i < num; i++) arr.push(start + step * i);
    return arr;
}

// --- Matrix math helpers (minimal) ---
function zeros(rows, cols) {
    const arr = [];
    for (let i = 0; i < rows; i++) arr.push(Array(cols).fill(0));
    return arr;
}

function addMatrixSlice(K, ke, dofs) {
    for (let i = 0; i < dofs.length; i++) {
        for (let j = 0; j < dofs.length; j++) {
            K[dofs[i]][dofs[j]] += ke[i][j];
        }
    }
}

// Solve linear system Ax = b using simple Gaussian elimination
function solveLinearSystem(A, b) {
    const n = b.length;
    const x = Array(n).fill(0);
    const M = A.map((row, i) => row.concat(b[i]));

    // Forward elimination with partial pivoting
    for (let k = 0; k < n; k++) {
        let i_max = k;
        for (let i = k + 1; i < n; i++) {
            if (Math.abs(M[i][k]) > Math.abs(M[i_max][k])) i_max = i;
        }
        [M[k], M[i_max]] = [M[i_max], M[k]];

        if (Math.abs(M[k][k]) < 1e-10) {
            console.warn("Singular matrix detected");
            return x; // return zeros
        }

        for (let i = k + 1; i < n; i++) {
            const f = M[i][k] / M[k][k];
            for (let j = k; j <= n; j++) M[i][j] -= M[k][j] * f;
        }
    }

    // Back substitution
    for (let i = n - 1; i >= 0; i--) {
        x[i] = M[i][n];
        for (let j = i + 1; j < n; j++) {
            x[i] -= M[i][j] * x[j];
        }
        x[i] /= M[i][i];
    }
    return x;
}

// --- FEM Beam functions ---
function beamElementStiffness(EI, L) {
    const k = [
        [12, 6 * L, -12, 6 * L],
        [6 * L, 4 * L * L, -6 * L, 2 * L * L],
        [-12, -6 * L, 12, -6 * L],
        [6 * L, 2 * L * L, -6 * L, 4 * L * L]
    ];
    return k.map(row => row.map(v => v * (EI / (L * L * L))));
}

// Trapezoidal integration helper
function trapzArrayUniform(y, dx) {
    if (!y || y.length <= 1) return 0;
    let sum = 0;
    for (let i = 0; i < y.length - 1; i++) sum += 0.5 * (y[i] + y[i + 1]);
    return sum * dx;
}

// Main solver
export function solveBeam(L, nElements, E, I, constraints = [], loads = []) {
    const nNodes = nElements + 1;
    const xs = linspace(0, L, nNodes);
    const totalDofs = nNodes * 2;
    const K = zeros(totalDofs, totalDofs);
    const F = Array(totalDofs).fill(0);

    // --- Assemble stiffness ---
    for (let e = 0; e < nElements; e++) {
        const Le = xs[e + 1] - xs[e];
        const ke = beamElementStiffness(E * I, Le);
        const dofs = [2 * e, 2 * e + 1, 2 * (e + 1), 2 * (e + 1) + 1];
        addMatrixSlice(K, ke, dofs);
    }

    // --- Apply loads ---
    const pointLoads = [];
    const pointMoments = [];
    const distLoads = [];

    for (const load of loads) {
        if (load[0] === "point") {
            const [_, x, Fy] = load;
            let idx = 0, minDist = Infinity;
            xs.forEach((xi, i) => {
                const d = Math.abs(xi - x);
                if (d < minDist) { minDist = d; idx = i; }
            });
            F[2 * idx] += Fy;
            pointLoads.push({ x, Fy });
        } else if (load[0] === "moment") {
            const [_, x, M] = load;
            let idx = 0, minDist = Infinity;
            xs.forEach((xi, i) => {
                const d = Math.abs(xi - x);
                if (d < minDist) { minDist = d; idx = i; }
            });
            F[2 * idx + 1] += M;
            pointMoments.push({ x, M });
        } else if (load[0] === "dist") {
            const [_, x0, x1, q] = load;
            const qFunc = typeof q === "function" ? q : () => q;
            distLoads.push({ x0, x1, qFunc });

            const nInt = 25;
            for (let e = 0; e < nElements; e++) {
                const xe1 = xs[e], xe2 = xs[e + 1];
                const Le = xe2 - xe1;
                const a = Math.max(xe1, x0);
                const b = Math.min(xe2, x1);
                if (b <= a) continue;

                const xi = Array.from({ length: nInt }, (_, i) => a + (b - a) * i / (nInt - 1));
                const s = xi.map(x => (x - xe1) / Le);
                const N1 = s.map(si => 1 - 3 * si ** 2 + 2 * si ** 3);
                const N2 = s.map(si => Le * (si - 2 * si ** 2 + si ** 3));
                const N3 = s.map(si => 3 * si ** 2 - 2 * si ** 3);
                const N4 = s.map(si => Le * (si ** 3 - si ** 2));
                const qvals = xi.map(qFunc);
                const dx = (b - a) / (nInt - 1);

                const f1 = trapzArrayUniform(qvals.map((qv, i) => qv * N1[i]), dx);
                const f2 = trapzArrayUniform(qvals.map((qv, i) => qv * N2[i]), dx);
                const f3 = trapzArrayUniform(qvals.map((qv, i) => qv * N3[i]), dx);
                const f4 = trapzArrayUniform(qvals.map((qv, i) => qv * N4[i]), dx);

                const fe = [f1, f2, f3, f4];
                const dofs = [2 * e, 2 * e + 1, 2 * (e + 1), 2 * (e + 1) + 1];
                for (let i = 0; i < 4; i++) F[dofs[i]] += fe[i];
            }
        }
    }

    // --- Constraints ---
    const constrained = new Set();        // DOFs with Dirichlet condition
    const prescribed = new Map();         // DOF → prescribed value (non-zero)

    for (const c of constraints) {
        let xpos, typ, data = {};
        if (Array.isArray(c)) {
            [xpos, typ, data = {}] = c;
            typ = String(typ).toUpperCase();
        } else continue;

        let idx = 0, minDist = Infinity;
        xs.forEach((xi, i) => {
            const d = Math.abs(xi - xpos);
            if (d < minDist) { minDist = d; idx = i; }
        });
        const vDof = 2 * idx;
        const tDof = 2 * idx + 1;

        if (["FIXED", "PIN", "ROLLER", "SIMPLE"].includes(typ)) {
            if (["FIXED"].includes(typ)) {
                constrained.add(vDof);
                constrained.add(tDof);
            } else {
                constrained.add(vDof);
            }
            // Zero displacement by default
            prescribed.set(vDof, 0);
            if (constrained.has(tDof)) prescribed.set(tDof, 0);
        } else if (typ === "ELASTIC") {
            if (data.ky) K[vDof][vDof] += data.ky;
            if (data.ktheta) K[tDof][tDof] += data.ktheta;
        } else if (typ === "PRESCRIBED") {
            if (data.v !== undefined) {
                constrained.add(vDof);
                prescribed.set(vDof, data.v);
            }
            if (data.theta !== undefined) {
                constrained.add(tDof);
                prescribed.set(tDof, data.theta);
            }
        }
    }

    // --- Solve system with non-homogeneous Dirichlet BCs ---
    const allDofs = Array.from({ length: totalDofs }, (_, i) => i);
    const freeDofs = allDofs.filter(d => !constrained.has(d));
    const consDofs = allDofs.filter(d => constrained.has(d));

    const Kff = freeDofs.map(i => freeDofs.map(j => K[i][j]));
    const Kfc = freeDofs.map(i => consDofs.map(j => K[i][j]));

    const Ff = freeDofs.map(i => F[i]);

    // Prescribed displacements
    const uc = consDofs.map(d => prescribed.get(d) ?? 0);

    // RHS = Ff - Kfc * uc
    const rhs = Ff.map((f, i) => {
        let val = f;
        for (let j = 0; j < consDofs.length; j++) {
            val -= Kfc[i][j] * uc[j];
        }
        return val;
    });

    // Solve for free DOFs
    let u = Array(totalDofs).fill(0);
    if (freeDofs.length > 0 && Kff.length > 0) {
        const uf = solveLinearSystem(Kff, rhs);
        freeDofs.forEach((d, i) => u[d] = uf[i]);
    }

    // Apply prescribed values
    consDofs.forEach((d, i) => u[d] = uc[i]);

    // --- Reactions (including elastic supports) ---
    const reactions = Array(totalDofs).fill(0);

    // For hard constraints (PIN/FIXED/PRESCRIBED): R = K * u - F
    consDofs.forEach((d, i) => {
        let r = -F[d];
        for (let j = 0; j < totalDofs; j++) {
            r += K[d][j] * u[j];
        }
        reactions[d] = r;
    });

    // For elastic supports: add spring forces
    for (const c of constraints) {
        if (c[1] === "ELASTIC" && c[2]) {
            let idx = 0, minDist = Infinity;
            xs.forEach((xi, i) => {
                const d = Math.abs(xi - c[0]);
                if (d < minDist) { minDist = d; idx = i; }
            });
            const vDof = 2 * idx;
            const tDof = 2 * idx + 1;
            if (c[2].ky) reactions[vDof] += -c[2].ky * u[vDof];
            if (c[2].ktheta) reactions[tDof] += -c[2].ktheta * u[tDof];
        }
    }

    // --- Post-processing ---
    const vs = u.filter((_, i) =>  i % 2 === 0);
    const scale = 1;

    // Internal forces (shear & moment)
    function q_total(x) {
        let s = 0;
        for (const d of distLoads) {
            if (x >= d.x0 && x <= d.x1) s += d.qFunc(x);
        }
        return s;
    }

    const samplesPerElement = 12;
    const totalSamples = nElements * samplesPerElement + 1;
    const x_samples = linspace(0, L, totalSamples);
    const dx = L / (totalSamples - 1);
    const eps = 1e-12;

    const q_vals = x_samples.map(q_total);
    const cumQ = [];
    let cum = 0; cumQ.push(0);
    for (let i = 1; i < x_samples.length; i++) {
        cum += 0.5 * (q_vals[i - 1] + q_vals[i]) * dx;
        cumQ.push(cum);
    }

    const pointLoadSum = x_samples.map(x => {
        let s = 0;
        for (const p of pointLoads) if (p.x < x - eps) s += p.Fy;
        return s;
    });

    const supportReactionSum = x_samples.map(x => {
        let s = 0;
        for (let i = 0; i < nNodes; i++) {
            if (xs[i] < x - eps && constrained.has(2 * i)) {
                s += reactions[2 * i] || 0;
            }
        }
        return s;
    });

    const shear_samples = x_samples.map((_, i) =>
        supportReactionSum[i] + pointLoadSum[i] - cumQ[i]
    );

    const cumV = [];
    cum = 0; cumV.push(0);
    for (let i = 1; i < x_samples.length; i++) {
        cum += 0.5 * (shear_samples[i - 1] + shear_samples[i]) * dx;
        cumV.push(cum);
    }

    const moment_samples = x_samples.map((_, i) => cumV[i]);

    return {
        x: xs,
        v_scaled: vs.map(v => v * scale),
        undeformed: Array(xs.length).fill(0),
        reactions,
        scale,
        x_samples,
        shear_samples,
        moment_samples
    };
}