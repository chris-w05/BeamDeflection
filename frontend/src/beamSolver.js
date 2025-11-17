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

    // Forward elimination
    for (let k = 0; k < n; k++) {
        let i_max = k;
        for (let i = k + 1; i < n; i++) if (Math.abs(M[i][k]) > Math.abs(M[i_max][k])) i_max = i;
        [M[k], M[i_max]] = [M[i_max], M[k]];

        for (let i = k + 1; i < n; i++) {
            const f = M[i][k] / M[k][k];
            for (let j = k; j <= n; j++) M[i][j] -= M[k][j] * f;
        }
    }

    // Back substitution
    for (let i = n - 1; i >= 0; i--) {
        x[i] = M[i][n] / M[i][i];
        for (let k = i - 1; k >= 0; k--) M[k][n] -= M[k][i] * x[i];
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

// Consistent load vector for uniform load
function uniformlyDistributedLoadVector(q, L) {
    const f = [6, L, 6, -L].map(v => v * q * L / 12);
    return f;
}

// Assemble global stiffness K and force F
export function solveBeam(L, nElements, E, I, constraints = [], loads = []) {
    const nNodes = nElements + 1;
    const xs = linspace(0, L, nNodes);
    const totalDofs = nNodes * 2;
    const K = zeros(totalDofs, totalDofs);
    const F = Array(totalDofs).fill(0);

    // --- assemble stiffness ---
    for (let e = 0; e < nElements; e++) {
        const Le = xs[e + 1] - xs[e];
        const ke = beamElementStiffness(E * I, Le);
        const dofs = [2 * e, 2 * e + 1, 2 * (e + 1), 2 * (e + 1) + 1];
        addMatrixSlice(K, ke, dofs);
    }

    // --- apply loads ---
    for (const load of loads) {
        if (load[0] === "point") {
            const [_, x, Fy] = load;
            let idx = 0;
            let minDist = Infinity;
            xs.forEach((xi, i) => { const d = Math.abs(xi - x); if (d < minDist) { minDist = d; idx = i; } });
            F[2 * idx] += Fy;
        }
        else if (load[0] === "moment") {
            const [_, x, M] = load;
            let idx = 0;
            let minDist = Infinity;
            xs.forEach((xi, i) => { const d = Math.abs(xi - x); if (d < minDist) { minDist = d; idx = i; } });
            F[2 * idx + 1] += M;
        }
        // Inside your solver, replace the dist load block with:
        else if (load[0] === "dist") {
            const [_, x0, x1, q] = load;
            const qFunc = typeof q === "function" ? q : (x) => q;

            const nIntegrationPoints = 25; // number of points per element
            for (let e = 0; e < nElements; e++) {
                const xe1 = xs[e], xe2 = xs[e + 1];
                const Le = xe2 - xe1;

                const a = Math.max(xe1, x0);
                const b = Math.min(xe2, x1);
                if (b <= a) continue;

                // Local coordinate array along overlap
                const xi = Array.from({ length: nIntegrationPoints }, (_, i) => a + (b - a) * i / (nIntegrationPoints - 1));
                const s = xi.map(xiVal => (xiVal - xe1) / Le);

                // Hermite shape functions
                const N1 = s.map(si => 1 - 3 * si ** 2 + 2 * si ** 3);
                const N2 = s.map(si => Le * (si - 2 * si ** 2 + si ** 3));
                const N3 = s.map(si => 3 * si ** 2 - 2 * si ** 3);
                const N4 = s.map(si => Le * (si ** 3 - si ** 2));

                // Evaluate q(x) at all integration points
                const qvals = xi.map(xiVal => qFunc(xiVal));

                // Approximate integral using trapezoidal rule
                const trapz = (y, dx) => {
                    let sum = 0;
                    for (let i = 0; i < y.length - 1; i++) {
                        sum += 0.5 * (y[i] + y[i + 1]) * dx;
                    }
                    return sum;
                };
                const dx = (b - a) / (nIntegrationPoints - 1);
                const f1 = trapz(qvals.map((qv, i) => qv * N1[i]), dx);
                const f2 = trapz(qvals.map((qv, i) => qv * N2[i]), dx);
                const f3 = trapz(qvals.map((qv, i) => qv * N3[i]), dx);
                const f4 = trapz(qvals.map((qv, i) => qv * N4[i]), dx);

                const fe = [f1, f2, f3, f4];
                const dofs = [2 * e, 2 * e + 1, 2 * (e + 1), 2 * (e + 1) + 1];
                for (let i = 0; i < 4; i++) F[dofs[i]] += fe[i];
            }
        }
    }

    // --- constraints ---
    const constrained = new Set();
    for (const [xpos, typ] of constraints) {
        let idx = 0;
        let minDist = Infinity;
        xs.forEach((xi, i) => { const d = Math.abs(xi - xpos); if (d < minDist) { minDist = d; idx = i; } });
        const t = typ.toUpperCase();
        if (["FIXED", "RIGID", "CLAMPED"].includes(t)) { constrained.add(2 * idx); constrained.add(2 * idx + 1); }
        else if (["PIN", "ROLLER", "SIMPLE"].includes(t)) { constrained.add(2 * idx); }
    }

    const allDofs = Array.from({ length: totalDofs }, (_, i) => i);
    const freeDofs = allDofs.filter(i => !constrained.has(i));

    // Partition matrices (Kff, Ff)
    const Kff = freeDofs.map(i => freeDofs.map(j => K[i][j]));
    const Ff = freeDofs.map(i => F[i]);

    let u = Array(totalDofs).fill(0);
    if (Kff.length > 0) {
        const uFree = solveLinearSystem(Kff, Ff);
        freeDofs.forEach((d, i) => u[d] = uFree[i]);
    }

    // reactions
    const reactions = K.map((row, i) => row.reduce((sum, v, j) => sum + v * u[j], 0)).map((v, i) => v - F[i]);

    const vs = u.filter((_, i) => i % 2 === 0); // vertical displacement
    let maxDef = Math.max(...vs.map(v => Math.abs(v)));

    // Scale set to 1 
    let scale = 1;//maxDef > 1e-8 ? 0.2 * L / maxDef : 1;

    return { x: xs, v_scaled: vs.map(v => v * scale), undeformed: Array(xs.length).fill(0), reactions, scale };
}
