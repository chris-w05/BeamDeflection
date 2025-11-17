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

// small trapezoidal helper over an array of y values with uniform dx
function trapzArrayUniform(y, dx) {
    if (!y || y.length <= 1) return 0;
    let sum = 0;
    for (let i = 0; i < y.length - 1; i++) sum += 0.5 * (y[i] + y[i + 1]);
    return sum * dx;
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

    // --- apply loads (build applied force vector F) ---
    const pointLoads = [];   // {x, Fy}
    const pointMoments = []; // {x, M}
    const distLoads = [];    // {x0, x1, qFunc}

    for (const load of loads) {
        if (load[0] === "point") {
            const [_, x, Fy] = load;
            let idx = 0, minDist = Infinity;
            xs.forEach((xi, i) => { const d = Math.abs(xi - x); if (d < minDist) { minDist = d; idx = i; } });
            F[2 * idx] += Fy;
            pointLoads.push({ x, Fy });
        } else if (load[0] === "moment") {
            const [_, x, M] = load;
            let idx = 0, minDist = Infinity;
            xs.forEach((xi, i) => { const d = Math.abs(xi - x); if (d < minDist) { minDist = d; idx = i; } });
            F[2 * idx + 1] += M;
            pointMoments.push({ x, M });
        } else if (load[0] === "dist") {
            const [_, x0, x1, q] = load;
            const qFunc = typeof q === "function" ? q : (x) => q;
            distLoads.push({ x0, x1, qFunc });

            const nIntegrationPoints = 25;
            for (let e = 0; e < nElements; e++) {
                const xe1 = xs[e], xe2 = xs[e + 1];
                const Le = xe2 - xe1;
                const a = Math.max(xe1, x0);
                const b = Math.min(xe2, x1);
                if (b <= a) continue;

                const xi = Array.from({ length: nIntegrationPoints }, (_, i) => a + (b - a) * i / (nIntegrationPoints - 1));
                const s = xi.map(xiVal => (xiVal - xe1) / Le);
                const N1 = s.map(si => 1 - 3 * si ** 2 + 2 * si ** 3);
                const N2 = s.map(si => Le * (si - 2 * si ** 2 + si ** 3));
                const N3 = s.map(si => 3 * si ** 2 - 2 * si ** 3);
                const N4 = s.map(si => Le * (si ** 3 - si ** 2));
                const qvals = xi.map(xiVal => qFunc(xiVal));
                const dx = (b - a) / (nIntegrationPoints - 1);

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

    const F_applied = F.slice();

    // --- constraints: PIN/ROLLER/SIMPLE -> v fixed; FIXED -> v, theta fixed ---
    const constrained = new Set(); // Dirichlet DOFs with value 0
    for (const [xpos, typRaw] of constraints) {
        const typ = String(typRaw || "").toUpperCase();
        let idx = 0, minDist = Infinity;
        xs.forEach((xi, i) => { const d = Math.abs(xi - xpos); if (d < minDist) { minDist = d; idx = i; } });
        const vDof = 2 * idx, tDof = 2 * idx + 1;
        if (["FIXED", "RIGID", "CLAMPED"].includes(typ)) { constrained.add(vDof); constrained.add(tDof); }
        else if (["PIN", "ROLLER", "SIMPLE"].includes(typ)) { constrained.add(vDof); }
    }

    // --- partition and solve indeterminate systems ---
    const allDofs = Array.from({ length: totalDofs }, (_, i) => i);
    const freeDofs = allDofs.filter(i => !constrained.has(i));
    const consDofs = allDofs.filter(i => constrained.has(i));

    const Kff = freeDofs.map(i => freeDofs.map(j => K[i][j]));
    const Kfc = freeDofs.map(i => consDofs.map(j => K[i][j]));
    const Kcf = consDofs.map(i => freeDofs.map(j => K[i][j]));
    const Kcc = consDofs.map(i => consDofs.map(j => K[i][j]));

    const Ff = freeDofs.map(i => F[i]);
    const Fc = consDofs.map(i => F[i]); // usually 0 unless you placed loads at constrained nodes

    // solve Kff * uf = Ff - Kfc * uc, with uc = 0
    let u = Array(totalDofs).fill(0);
    if (Kff.length > 0) {
        const rhs = Ff.slice(); // uc = 0 => rhs = Ff
        const uf = solveLinearSystem(Kff, rhs);
        freeDofs.forEach((d, i) => u[d] = uf[i]);
    }
    // uc = 0 at constrained DOFs (already)
    consDofs.forEach(d => u[d] = 0);

    // reactions at constrained DOFs: Rc = Kcf*uf + Kcc*uc - Fc = Kcf*uf - Fc
    const reactions = Array(totalDofs).fill(0);
    consDofs.forEach((d, iRow) => {
        let r = -Fc[iRow];
        for (let j = 0; j < freeDofs.length; j++) r += Kcf[iRow][j] * u[freeDofs[j]];
        // Kcc*uc is 0 since uc=0
        reactions[d] = r;
    });

    // displacements for plotting
    const vs = u.filter((_, i) => i % 2 === 0);
    const scale = 1;

    // --- INTERNAL FORCE (shear & moment) using your fixed signs ---
    function q_total(x) {
        let s = 0;
        for (const d of distLoads) if (x >= d.x0 && x <= d.x1) s += d.qFunc(x);
        return s;
    }

    const samplesPerElement = 12;
    const totalSamples = nElements * samplesPerElement + 1;
    const x_samples = linspace(0, L, totalSamples);
    const dx_global = L / (totalSamples - 1);
    const eps = 1e-12;

    const q_vals = x_samples.map(xi => q_total(xi));
    const cumulative_q = [];
    {
        let cum = 0; cumulative_q.push(0);
        for (let i = 1; i < x_samples.length; i++) {
            const area = 0.5 * (q_vals[i - 1] + q_vals[i]) * dx_global;
            cum += area; cumulative_q.push(cum);
        }
    }

    const pLoads = pointLoads.map(p => ({ x: p.x, val: p.Fy }));
    const pMoms = pointMoments.map(p => ({ x: p.x, val: p.M }));

    function isAtConstrainedNode(x) {
        for (let j = 0; j < nNodes; j++) {
            if (Math.abs(x - xs[j]) < eps && (constrained.has(2 * j) || constrained.has(2 * j + 1))) return true;
        }
        return false;
    }

    const cumulativeAppliedPointLoadsAtSample = x_samples.map(xi => {
        let s = 0;
        for (const pl of pLoads) if (!isAtConstrainedNode(pl.x) && pl.x < xi - eps) s += pl.val;
        return s;
    });
    const cumulativeAppliedMomAtSample = x_samples.map(xi => {
        let s = 0;
        for (const pm of pMoms) if (!isAtConstrainedNode(pm.x) && pm.x < xi - eps) s += pm.val;
        return s;
    });

    const cumulativeSupportReactionAtSample = x_samples.map(xi => {
        let s = 0;
        for (let i = 0; i < nNodes; i++) {
            if (xs[i] < xi - eps && constrained.has(2 * i)) s += reactions[2 * i] || 0;
        }
        return s;
    });
    const cumulativeSupportReactionMomAtSample = x_samples.map(xi => {
        let s = 0;
        for (let i = 0; i < nNodes; i++) {
            if (xs[i] < xi - eps && constrained.has(2 * i + 1)) s += reactions[2 * i + 1] || 0;
        }
        return s;
    });

    // Use the sign convention that made your diagrams correct
    const shear_samples = [];
    for (let i = 0; i < x_samples.length; i++) {
        const Vx = cumulativeSupportReactionAtSample[i]
            + cumulativeAppliedPointLoadsAtSample[i]
            - cumulative_q[i];
        shear_samples.push(Vx);
    }

    const cumulativeIntegralV = [];
    {
        let cumV = 0; cumulativeIntegralV.push(0);
        for (let i = 1; i < x_samples.length; i++) {
            const areaV = 0.5 * (shear_samples[i - 1] + shear_samples[i]) * dx_global;
            cumV += areaV; cumulativeIntegralV.push(cumV);
        }
    }

    const moment_samples = [];
    for (let i = 0; i < x_samples.length; i++) {
        const Mx = - cumulativeSupportReactionMomAtSample[i]
            - cumulativeAppliedMomAtSample[i]
            + cumulativeIntegralV[i];
        moment_samples.push(Mx);
    }

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


