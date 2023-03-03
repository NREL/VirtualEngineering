import numpy as np

try:
    import dolfinx
    import ufl
    from petsc4py import PETSc
except:
    pass


def smoothstep(in_vec, lo_edge, hi_edge, lo_val, hi_val):
    out_vec = np.zeros(np.shape(in_vec))

    for k, x in enumerate(in_vec):
        if x < lo_edge:
            out_vec[k] = lo_val
        elif x >= hi_edge:
            out_vec[k] = hi_val
        else:
            xs = (x - lo_edge) / (hi_edge - lo_edge)

            out_vec[k] = (hi_val - lo_val) * xs * xs * xs * (
                xs * (xs * 6.0 - 15.0) + 10.0
            ) + lo_val

    return out_vec


def linstep(in_vec, lo_edge, hi_edge, lo_val, hi_val):
    out_vec = np.zeros(np.shape(in_vec))

    for k, x in enumerate(in_vec):
        if x <= lo_edge:
            out_vec[k] = lo_val
        elif x > hi_edge:
            out_vec[k] = hi_val
        else:
            xs = (x - lo_edge) / (hi_edge - lo_edge)
            out_vec[k] = lo_val + (hi_val - lo_val) * xs

    return out_vec


def project(v, target_func, bcs=[]):
    # Ensure we have a mesh and attach to measure
    V = target_func.function_space
    dx = ufl.dx(V.mesh)

    # Define variational problem for projection
    w = ufl.TestFunction(V)
    Pv = ufl.TrialFunction(V)
    a = dolfinx.fem.form(ufl.inner(Pv, w) * dx)
    L = dolfinx.fem.form(ufl.inner(v, w) * dx)

    # Assemble linear system
    A = dolfinx.fem.petsc.assemble_matrix(a, bcs)
    A.assemble()
    b = dolfinx.fem.petsc.assemble_vector(L)
    dolfinx.fem.petsc.apply_lifting(b, [a], [bcs])
    b.ghostUpdate(addv=PETSc.InsertMode.ADD, mode=PETSc.ScatterMode.REVERSE)
    dolfinx.fem.petsc.set_bc(b, bcs)

    # Solve linear system
    solver = PETSc.KSP().create(A.getComm())
    solver.setOperators(A)
    solver.solve(b, target_func.vector)

    return target_func
