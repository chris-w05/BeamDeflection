# backend/api/BeamDeflection.py

# Finite-element Euler-Bernoulli 2D beam solver and plot
# - Uses 2-node Hermitian beam elements (2 DOF per node: vertical displacement v and rotation theta)
# - Supports constraints placed at arbitrary x (snapped to nearest node): 'FIXED' (v=0,theta=0), 'PIN' (v=0), 'FREE' (none)
# - Loads supported:
#    ('point', x, Fy)     -> point vertical load (positive upward)
#    ('moment', x, M)     -> point moment (positive counter-clockwise)
#    ('dist', x0, x1, q)  -> uniform distributed load q between x0 and x1 (q positive upward)
# - Example usage at bottom demonstrates a simply supported beam with a mid-span point load
#
# Note: This is a linear small-deflection Euler-Bernoulli beam FEM solver (static). For interior supports, we
# snap supports to the nearest node; increase `n_elements` to place supports more accurately.

import numpy as np
import matplotlib.pyplot as plt
from math import sin, cos, pi, log, sqrt, tan, asin, acos, atan
from scipy.integrate import simpson

def beam_element_stiffness(EI, L):
    """4x4 element stiffness for Euler-Bernoulli beam (v,theta at node1 and node2)"""
    k = (EI / L**3) * np.array([
        [12, 6*L, -12, 6*L],
        [6*L, 4*L**2, -6*L, 2*L**2],
        [-12, -6*L, 12, -6*L],
        [6*L, 2*L**2, -6*L, 4*L**2]
    ], dtype=float)
    return k

def uniformly_distributed_load_vector(q, L):
    """Consistent load vector for uniform q (force per length), positive upward"""
    # f_e = q * L / 12 * [6, L, 6, -L]
    f = (q * L / 12.0) * np.array([6.0, L, 6.0, -L])
    return f

def assemble_beam(L, n_elements, E, I, constraints, loads):
    """
    Assemble global stiffness matrix K and force vector F for 2D Euler-Bernoulli beam.
    
    - constraints: list of (x_pos, type) tuples; type in {'FIXED','PIN','FREE','CLAMPED','RIGID','ROLLER','SIMPLE'}
    - loads: list of 
        ('point', x, Fy)          # Fy positive = upward
        ('moment', x, M)          # M positive = CCW
        ('dist', x0, x1, q)       # q float (uniform) or callable q(x); positive upward; applied only on [x0,x1]
    Returns: xs, K, F, constrained_dofs_list
    """
    n_nodes = n_elements + 1
    xs = np.linspace(0, L, n_nodes)
    total_dofs = n_nodes * 2
    K = np.zeros((total_dofs, total_dofs))
    F = np.zeros(total_dofs)

    # --- Element stiffness assembly ---
    for e in range(n_elements):
        x1, x2 = xs[e], xs[e+1]
        Le = x2 - x1
        ke = beam_element_stiffness(E * I, Le)
        dofs = [2*e, 2*e + 1, 2*(e + 1), 2*(e + 1) + 1]
        K[np.ix_(dofs, dofs)] += ke  # Vectorized addition

    # --- Distributed loads (unified numerical integration - accurate for uniform & variable q(x)) ---
    for load in loads:
        if load[0] == 'dist':
            _, x0, x1, q = load
            q_func = q if callable(q) else lambda x: float(q)  # Ensure callable
            for e in range(n_elements):
                xe1, xe2 = xs[e], xs[e+1]
                Le = xe2 - xe1
                a = max(xe1, x0)
                b = min(xe2, x1)
                if b <= a:
                    continue

                # Local coordinate array over *overlap only*
                local_a = a - xe1
                local_b = b - xe1
                xi = np.linspace(local_a, local_b, 25)  # 25 points → excellent accuracy for polynomials ≤ degree 3
                xvals = xe1 + xi
                s = xi / Le

                # Hermite shape functions (standard convention)
                N1 = 1 - 3*s**2 + 2*s**3
                N2 = Le * (s - 2*s**2 + s**3)
                N3 = 3*s**2 - 2*s**3
                N4 = Le * (s**3 - s**2)  # Corrected sign & formula

                qvals = q_func(xvals)

                # Integrate (Simpson's rule if available, else trapz)
                try:
                    f1 = simpson(qvals * N1, xi)
                    f2 = simpson(qvals * N2, xi)
                    f3 = simpson(qvals * N3, xi)
                    f4 = simpson(qvals * N4, xi)
                except NameError:
                    f1 = np.trapz(qvals * N1, xi)
                    f2 = np.trapz(qvals * N2, xi)
                    f3 = np.trapz(qvals * N3, xi)
                    f4 = np.trapz(qvals * N4, xi)

                fe = np.array([f1, f2, f3, f4])
                dofs = [2*e, 2*e + 1, 2*(e + 1), 2*(e + 1) + 1]
                F[dofs] += fe

    # --- Point loads & moments (snapped to nearest node) ---
    for load in loads:
        if load[0] == 'point':
            _, x, Fy = load
            idx = np.argmin(np.abs(xs - x))
            F[2 * idx] += Fy
        elif load[0] == 'moment':
            _, x, M = load
            idx = np.argmin(np.abs(xs - x))
            F[2 * idx + 1] += M

    # --- Constraints (snapped to nearest node) ---
    constrained = []
    seen = set()
    for xpos, typ in constraints:
        idx = np.argmin(np.abs(xs - xpos))
        key = 2 * idx  # Use even DOF as key to avoid duplicate full constraints at same node
        if idx in seen:
            continue
        typ = typ.upper()
        if typ in {'FIXED', 'CLAMPED', 'RIGID'}:
            constrained += [2*idx, 2*idx + 1]
            seen.add(idx)
        elif typ in {'PIN', 'ROLLER', 'SIMPLE'}:
            constrained += [2*idx]
            seen.add(idx)
        # FREE or unknown → no DOFs constrained

    constrained = sorted(set(constrained))
    return xs, K, F, constrained

def solve_beam(L, n_elements, E, I, constraints, loads):
    xs, K, F, constrained = assemble_beam(L, n_elements, E, I, constraints, loads)
    total_dofs = K.shape[0]
    all_dofs = np.arange(total_dofs)
    free = np.setdiff1d(all_dofs, constrained)
    # Partition and solve
    Kff = K[np.ix_(free, free)]
    Kfc = K[np.ix_(free, constrained)]
    Ff = F[free]
    # For homogeneous prescribed displacements = 0 so no extra terms
    # Solve for free DOFs
    u = np.zeros(total_dofs)
    if Kff.size == 0:
        print("All DOFs constrained.")
        return xs, u, np.zeros_like(u)
    sol = np.linalg.solve(Kff, Ff)
    u[free] = sol
    reactions = K @ u - F
    return xs, u, reactions

def plot_beam(xs, u, scale=1.0, ax=None, show=True):
    """
    Plot undeformed and deformed beam.
    u: global DOF vector; vertical DOF at even indices, rotations at odd indices
    scale: factor to magnify deformation for visualization (if None, auto-scale)
    """
    vs = u[0::2]
    thetas = u[1::2]
    if scale is None:
        # auto scale such that max deflection is ~10% of beam length
        max_def = np.max(np.abs(vs))
        if max_def < 1e-12:
            scale = 1.0
        else:
            scale = 0.1 * (xs[-1] - xs[0]) / max_def

    x_dense = np.linspace(xs[0], xs[-1], max(200, len(xs)*20))
    # build cubic Hermite interpolation from nodal v and theta
    v_interp = np.zeros_like(x_dense)
    for e in range(len(xs)-1):
        x1 = xs[e]; x2 = xs[e+1]; Le = x2 - x1
        # nodal values
        v1 = vs[e]; t1 = thetas[e]
        v2 = vs[e+1]; t2 = thetas[e+1]
        xe = (x_dense >= x1) & (x_dense <= x2)
        xi = x_dense[xe] - x1
        s = xi / Le
        H1 = 1 - 3*s**2 + 2*s**3
        H2 = xi * (1 - 2*s + s**2)  # multiplied by Le for proper scaling
        H3 = 3*s**2 - 2*s**3
        H4 = xi * ( -s**2 + s)      # *Le scaling
        v_local = H1*v1 + H2*(t1/1.0) + H3*v2 + H4*(t2/1.0)
        v_interp[xe] = v_local

    if ax is None:
        fig, ax = plt.subplots(figsize=(8,3))
    ax.plot(xs, np.zeros_like(xs), 'k--', label='undeformed')
    ax.plot(x_dense, v_interp*scale, '-', linewidth=2, label=f'deformed scale: {scale:.2g}')
    ax.set_xlabel('x')
    ax.set_ylabel('vertical displacement (scaled)')
    ax.legend()
    ax.grid(True)
    if show:
        plt.show()

# Example usage
if __name__ == "__main__":
    # beam properties
    L = 6.0            # length (m)
    E = 210e9          # Young's modulus (Pa)
    I = 8.1e-6         # second moment area (m^4) - adjust per cross-section
    n_elements = 200    # mesh refinement

    # supports: simply supported at x=0 and x=L (PINs)
    constraints = [
        (0.0, 'rigid'),
        # (L, 'Pin')
    ]

    # loads: mid-span downward point load of 10 kN, plus uniform downward load 0.0
    loads = [
        ('moment', L/2.0, 10000.0),
        # ('point', 1.9*L, -10000.0),# downward 10kN at midspan (negative is downward)
        #('dist', 0.0, L, -2000.0)         # optional: uniform downward load -2000 N/m along whole beam
    ]

    xs, u, reactions = solve_beam(L, n_elements, E, I, constraints, loads)
    # print("Nodal vertical displacements (m):")
    # for i, x in enumerate(xs):
    #     print(f"x={x:.3f} m: v={u[2*i]:.6e} m, theta={u[2*i+1]:.6e} rad")

    # auto-scale: choose scale so that max deflection visible
    max_def = np.max(np.abs(u[0::2]))
    if max_def < 1e-12:
        scale = 1.0
    else:
        scale = 0.2 * L / max_def
    plot_beam(xs, u, scale=scale)
