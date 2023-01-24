# fun comment
from fenics import *
import numpy as np
import matplotlib.pyplot as plt
from PretreatmentVisualizer import update_figure_1, update_figure_2, finalize_figures
from Utilities import linstep, smoothstep


class Pretreatment:
    def __init__(self):

        set_log_level(LogLevel.ERROR)

        self._init_constants()

    # ================================================================

    def _init_constants(self):

        # Physical Constants
        self.k_b = 1.380649e-23
        self.eta = 2.8e-4  # dynamic viscosity of solvent (water)
        self.R = 8.31446261815324

        # Diffusion coefficient parameters (originally nano-meters)
        self.d_pore = 15.0 / 1.0e9  # pore diameter
        self.r_xy = 0.34 / 1.0e9  # xylose molecular radius
        self.r_xo = 1.02 / 1.0e9  # xylooligomer molecular radius
        self.r_f = 0.28 / 1.0e9  # furfural molecular radius

        # Thermal conductivities
        self.k_s = 0.12  # solid
        self.k_l = 0.63  # liquid
        self.k_g = 0.0197  # gas
        self.h = 5000.0  # convective heat transfer coeff

        # Specific heat capacities
        self.C_s = 2000.0  # solid
        self.C_l = 4350.0  # liquid
        self.C_g = 1970.0  # gas

        # Mass densities
        self.rho_s = 1000.0  # solid
        self.rho_l = 1000.0  # liquid
        self.rho_g = 2.01928202521  # gas
        self.rho_x = 730.0  # xlyan
        self.rho_gl = 1500.0  # glucan
        self.rho_li = 1300.0  # lignin
        self.rho_os = 0.5 * (self.rho_gl + self.rho_li)  # other solids

        # Molecular weights (originally g/mol)
        self.M_x = 132.0 / 1.0e3  # xylan
        self.M_xy = 150.0 / 1.0e3  # xylose
        self.M_xo = 450.0 / 1.0e3  # xylooligomer
        self.M_f = 96.0 / 1.0e3  # furfural
        self.M_w = 18.01528 / 1.0e3  # water

        # Latent heat of steam
        self.L_cond = 2200.0e3

        # Condensation and evaporation rate (s^-1)
        self.k_bar = 10.0

    # ================================================================

    def generate_mesh(self, nn=15):
        self.length = 5e-3 / 2.0
        self.diam = 2.0 * self.length / 20.0

        self.mesh = IntervalMesh(nn, 0.0, self.length)
        self.ndim = self.mesh.topology().dim()

    # ================================================================

    def build_functions(self, degree=2):
        d1p3 = FiniteElement("P", interval, degree)
        self.S = FunctionSpace(self.mesh, d1p3)

        self.num_comp = 7

        d1p3_mix = MixedElement([d1p3 for k in range(self.num_comp)])
        self.Q = FunctionSpace(self.mesh, d1p3_mix)

        self.v = TestFunctions(self.Q)

        self.u = Function(self.Q)
        self.u_n = Function(self.Q)

        self.c_s_id = 0
        self.eps_l_id = 1
        self.fstar_x_id = 2
        self.cstar_xo_id = 3
        self.cstar_xy_id = 4
        self.cstar_f_id = 5
        self.T_id = 6

    # ================================================================

    def build_problem(self):

        self._define_reaction_rates()

        self._set_initial_conditions()
        self._set_boundary_conditions()

        self._build_general_forms()

        self._build_convection_terms()
        self._build_diffusion_terms()
        self._build_reaction_terms()
        self._build_forcing_terms()

        self._build_nonlinear_form()

    # ================================================================

    def _define_reaction_rates(self):
        # Reaction rates (originally in A:cm^3/mol/s and E_A:kJ/mol)
        self.k_xo = 1.0e9 * exp(-110.0e3 / (self.R * self.u[self.T_id]))
        self.k_x1 = 8.0e11 * exp(-130.0e3 / (self.R * self.u[self.T_id]))
        self.k_x2 = 2.5e8 * exp(-110.0e3 / (self.R * self.u[self.T_id]))
        self.k_f = 7.0e5 * exp(-98.0e3 / (self.R * self.u[self.T_id]))

    # ================================================================

    def _set_initial_conditions(self):
        self.c_acid0 = 0.1 * 1e3  # Convert mol/self.length to mol/m^3
        self.c_sbulk = 0.14 * 1e3

        self.f_x0 = 0.26

        self.eps_p0 = 0.8
        self.eps_lt = self.eps_p0
        self.eps_l0 = 0.25

        self.f_is0 = 0.44

        self.T_0 = 300.0
        self.T_s = 423.0

        ic = Function(self.S)

        # u[self.eps_l_id]
        ic.vector()[:] = self.eps_l0
        assign(self.u.sub(self.eps_l_id), ic)
        assign(self.u_n.sub(self.eps_l_id), ic)

        # u[self.fstar_x_id]
        fstar_x0 = self.f_x0 * (1.0 - self.eps_p0)
        ic.vector()[:] = fstar_x0
        assign(self.u.sub(self.fstar_x_id), ic)
        assign(self.u_n.sub(self.fstar_x_id), ic)

        # u[self.T_id]
        ic.vector()[:] = self.T_0
        assign(self.u.sub(self.T_id), ic)
        assign(self.u_n.sub(self.T_id), ic)

    # ================================================================

    def _set_boundary_conditions(self):

        # TODO: there's probably a better way of making this accessible from subdomain
        domain_length = self.length
        tol = 1e-12

        class x0_lo(SubDomain):
            def inside(self, x, on_boundary):
                return on_boundary and x[0] < tol

        class x0_hi(SubDomain):
            def inside(self, x, on_boundary):
                return on_boundary and x[0] > (domain_length - tol)

        self.bcs = []
        self.bcs.append(
            DirichletBC(self.Q.sub(self.c_s_id), Constant(self.c_sbulk), x0_hi(123))
        )
        self.bcs.append(
            DirichletBC(self.Q.sub(self.T_id), Constant(self.T_s), x0_hi())
        )  # TODO: remove this testing term once heat flux is working

        boundary_markers = MeshFunction("size_t", self.mesh, self.ndim - 1)
        boundary_markers.set_all(2)
        x0_lo().mark(boundary_markers, 0)
        x0_hi().mark(boundary_markers, 1)

        self.ds = Measure("ds", domain=self.mesh, subdomain_data=boundary_markers)

    # ================================================================

    def _build_general_forms(self):

        aa = self.rho_x / self.rho_os
        # aa = 2.0 # TODO: This is set as 2.0 in the paper/OG code
        bb = self.f_x0 / (1.0 - self.f_x0) + aa
        cc = 1.0 - self.eps_p0

        ss1 = aa**2 * cc**2
        ss2 = 2.0 * (aa - 2.0) * bb * cc * self.u_n[self.fstar_x_id]
        ss3 = bb**2 * self.u_n[self.fstar_x_id] ** 2

        ss = sqrt(ss1 - ss2 + ss3)

        self.eps_p = -(aa * cc + bb * (self.u_n[self.fstar_x_id] - 2.0) + ss) / (
            2.0 * bb
        )
        self.eps_g = self.eps_p - self.u_n[self.eps_l_id]

        # --------------------------------

        self.c_acid = self.c_acid0 * self.eps_l0 / self.u_n[self.eps_l_id]

        # --------------------------------

        self.D_s = (
            self.eps_g
            * self.d_pore
            * sqrt((8.0 * self.R * self.u_n[self.T_id]) / (pi * self.M_w))
            / 3.0
        )

        self.D_xy_0 = (self.k_b * self.u_n[self.T_id]) / (
            6.0 * pi * self.eta * self.r_xy
        )
        self.D_xo_0 = (self.k_b * self.u_n[self.T_id]) / (
            6.0 * pi * self.eta * self.r_xo
        )
        self.D_f_0 = (self.k_b * self.u_n[self.T_id]) / (6.0 * pi * self.eta * self.r_f)

        self.D_xy = self.u_n[self.eps_l_id] * self.D_xy_0
        self.D_xo = self.u_n[self.eps_l_id] * self.D_xo_0
        self.D_f = self.u_n[self.eps_l_id] * self.D_f_0

        # --------------------------------

        self.rho_eff_c_eff = (
            (1.0 - self.eps_p) * self.rho_s * self.C_s
            + self.u_n[self.eps_l_id] * self.rho_l * self.C_l
            + self.eps_g * self.rho_g * self.C_g
        )
        self.k_eff = (
            (1.0 - self.eps_p) * self.k_s
            + self.u_n[self.eps_l_id] * self.k_l
            + self.eps_g * self.k_g
        )

        self.heat_flux = (
            self.h
            / self.k_eff
            * (self.u_n[self.T_id] - self.T_s)
            * self.v[self.T_id]
            * self.ds(1)
        )  # TODO: fix heat flux specification, currently not used in formulation

        # --------------------------------

        self.k_cond_ufl = conditional(
            le(self.u_n[self.T_id], self.T_s), self.k_bar, 0.0
        )
        self.k_evap_ufl = conditional(
            ge(self.u_n[self.T_id], self.T_s),
            conditional(gt(self.u_n[self.eps_l_id], self.eps_lt), self.k_bar, 0.0),
            0.0,
        )

        self.k_cond = Function(self.S)
        self.k_cond.vector()[:] = self.k_bar

        self.k_evap = Function(self.S)
        self.k_evap.vector()[:] = 0.0

    # ================================================================

    def _build_convection_terms(self):
        self.conv = [0 for k in range(self.num_comp)]

        self.conv[self.c_s_id] = 0
        self.conv[self.eps_l_id] = 0
        self.conv[self.fstar_x_id] = 0
        self.conv[self.cstar_xo_id] = (
            self.D_xo_0 / self.u_n[self.eps_l_id] * Dx(self.u_n[self.eps_l_id], 0)
        )
        self.conv[self.cstar_xy_id] = (
            self.D_xy_0 / self.u_n[self.eps_l_id] * Dx(self.u_n[self.eps_l_id], 0)
        )
        self.conv[self.cstar_f_id] = (
            self.D_f_0 / self.u_n[self.eps_l_id] * Dx(self.u_n[self.eps_l_id], 0)
        )
        # self.conv[self.cstar_xo_id] = self.D_xo_0 / self.u_n[self.eps_l_id] * grad(self.u_n[self.eps_l_id])
        # self.conv[self.cstar_xy_id] = self.D_xy_0 / self.u_n[self.eps_l_id] * grad(self.u_n[self.eps_l_id])
        # self.conv[self.cstar_f_id] = self.D_f_0 / self.u_n[self.eps_l_id] * grad(self.u_n[self.eps_l_id])
        self.conv[self.T_id] = 0

    # ================================================================

    def _build_diffusion_terms(self):
        self.diff = [0 for k in range(self.num_comp)]

        self.diff[self.c_s_id] = self.D_s
        self.diff[self.eps_l_id] = 0
        self.diff[self.fstar_x_id] = 0
        self.diff[self.cstar_xo_id] = self.D_xo_0
        self.diff[self.cstar_xy_id] = self.D_xy_0
        self.diff[self.cstar_f_id] = self.D_f_0
        self.diff[self.T_id] = self.k_eff / self.rho_eff_c_eff

    # ================================================================

    def _build_reaction_terms(self):
        self.reac = [0 for k in range(self.num_comp)]

        self.reac[self.c_s_id] = -self.k_cond * self.eps_g
        self.reac[self.eps_l_id] = (
            -self.M_w * self.k_cond * self.u_n[self.c_s_id] / self.rho_l - self.k_evap
        )
        self.reac[self.fstar_x_id] = -(self.k_xo + self.k_x1) * self.c_acid * self.eps_p
        self.reac[self.cstar_xo_id] = -self.k_x2 * self.c_acid
        self.reac[self.cstar_xy_id] = -self.k_f * self.c_acid
        self.reac[self.cstar_f_id] = 0
        self.reac[self.T_id] = -4.0 * self.h / (self.rho_eff_c_eff * self.diam)
        # self.reac[self.T_id] = -4.0 * self.h / self.diam

    # ================================================================

    def _build_forcing_terms(self):
        self.forc = [0 for k in range(self.num_comp)]

        self.forc[self.c_s_id] = (
            self.k_evap
            * (self.u_n[self.eps_l_id] - self.eps_lt)
            * self.rho_l
            / self.M_w
        )
        self.forc[self.eps_l_id] = (
            self.k_cond * self.u_n[self.c_s_id] * self.eps_p * self.M_w / self.rho_l
            + self.k_evap * self.eps_lt
        )
        self.forc[self.fstar_x_id] = 0
        self.forc[self.cstar_xo_id] = (
            self.k_xo
            * self.rho_s
            * self.c_acid
            * self.eps_p
            * self.u_n[self.fstar_x_id]
            / self.M_xo
        )
        self.forc[self.cstar_xy_id] = (
            self.k_x1
            * self.rho_s
            * self.c_acid
            * self.eps_p
            * self.u_n[self.fstar_x_id]
            / self.M_xy
            + self.k_x2
            * self.u_n[self.cstar_xo_id]
            * self.c_acid
            * self.M_xo
            / self.M_xy
        )
        self.forc[self.cstar_f_id] = (
            self.k_f * self.u_n[self.cstar_xy_id] * self.c_acid * self.M_xy / self.M_f
        )
        self.forc[self.T_id] = (
            self.L_cond
            / self.rho_eff_c_eff
            * (
                self.k_cond * self.u_n[self.c_s_id] * self.M_w * (self.eps_g)
                - self.k_evap * self.rho_l * (self.u_n[self.eps_l_id] - self.eps_lt)
            )
            + 4.0 * self.h / (self.rho_eff_c_eff * self.diam) * self.T_s
        )
        # self.forc[self.T_id] = (
        #     (self.L_cond / self.rho_eff_c_eff)
        #     * (
        #         self.k_cond * self.u_n[self.c_s_id] * self.M_w * (self.eps_g)
        #         - self.k_evap * self.rho_l * (self.u_n[self.eps_l_id] - self.eps_lt)
        #     )
        #     + 4.0 * self.h / self.diam * self.T_s
        # )

    # ================================================================

    def _build_nonlinear_form(self):
        self.F = 0

        # Use a placeholder value of dt to build the form, it will
        # be updated with the correct value in solve(...)
        self.dt_c = Constant(1.0)

        for k in range(self.num_comp):

            self.F += (1.0 / self.dt_c) * (self.u[k] - self.u_n[k]) * self.v[k] * dx

            if self.conv[k] != 0:
                # BEST ONE (move the grad onto the test space and swap sign)
                self.F -= self.conv[k] * self.u[k] * Dx(self.v[k], 0) * dx

                # self.F += Dx(self.conv[k]*self.u[k], 0) * self.v[k] * dx
                # self.F -= inner(self.conv[k]*self.u[k], grad(self.v[k]))*dx
                # self.F += dot(self.conv[k], grad(self.u[k]))*self.v[k]*dx

            if self.diff[k] != 0:
                # BEST ONE (move the grad onto the test space and swap sign)
                self.F += self.diff[k] * Dx(self.u[k], 0) * Dx(self.v[k], 0) * dx

                # self.F += self.diff[k] * inner(grad(self.u[k]), grad(self.v[k])) * dx
                # self.F += self.diff[k] * dot(grad(self.u[k]), grad(self.v[k])) * dx
                # self.F += inner(self.diff[k] * grad(self.u[k]), grad(self.v[k])) * dx

            if self.reac[k] != 0:
                self.F -= self.reac[k] * self.u[k] * self.v[k] * dx

            if self.forc[k] != 0:
                self.F -= self.forc[k] * self.v[k] * dx

        # TODO: fix heat flux specification, currently not used in formulation
        # self.F += self.heat_flux

    # ================================================================

    def solve(self, t_final=30.0, dt=1.0, save_every=2.0):

        # Time stepping parameters
        self.dt = dt
        self.dt_c.assign(self.dt)

        self.t_final = t_final
        t_steps = int(self.t_final / self.dt)
        save_every_n = int(save_every / self.dt)

        tt = 0
        self._save_solution(tt)

        plot_line_at = np.array([self.dt, 30.0, 300.0, 600.0, 1200.0])
        # plot_line_at = 60.0*np.array([1.7, 5.0, 10.0, 20.0])

        du = TrialFunction(self.Q)
        J  = derivative(self.F, self.u, du)

        nonlinear_problem = NonlinearVariationalProblem(self.F, self.u, self.bcs, J)
        nonlinear_solver  = NonlinearVariationalSolver(nonlinear_problem)

        # Set some of the solver options
        solver_parameters = nonlinear_solver.parameters
        solver_parameters["nonlinear_solver"] = "snes"
        # solver_parameters["snes_solver"]["linear_solver"] = "gmres"
        # solver_parameters["snes_solver"]["preconditioner"] = "jacobi"
        solver_parameters["snes_solver"]["report"] = False

        T_dofs = self.Q.sub(self.T_id).dofmap().dofs()
        eps_l_dofs = self.Q.sub(self.eps_l_id).dofmap().dofs()

        for k in range(t_steps):
            # Solve the nonlinear problem
            nonlinear_solver.solve()

            # Extract a numpy vector for temp and eps_l
            np_T = self.u.vector().vec()[T_dofs]
            np_eps_l = self.u.vector().vec()[eps_l_dofs]

            self._update_k_cond(np_T)
            self._update_k_evap(np_T, np_eps_l)

            # Update previous solution
            self.u_n.assign(self.u)

            tt = (k + 1) * self.dt

            if (k + 1) % save_every_n == 0:
                self._save_solution(tt)

            if np.amin(np.abs(tt - plot_line_at)) < 1e-6:
                line_id = np.argmin(np.abs(tt - plot_line_at))
                self = update_figure_1(self, tt, line_id)
                self = update_figure_2(self, tt, line_id)

        self = finalize_figures(self)

    # ================================================================

    def _update_k_cond(self, temp):

        T_tol = 0.001
        del_T = T_tol * self.T_s
        out_vec = linstep(temp, self.T_s - del_T, self.T_s + del_T, self.k_bar, 0.0)

        self.k_cond.vector()[:] = out_vec

    # ================================================================

    def _update_k_evap(self, temp, liq):

        T_tol = 0.001
        del_T = T_tol * self.T_s
        out_vec = linstep(temp, self.T_s - del_T, self.T_s + del_T, 0.0, self.k_bar)

        eps_lt_tol = 0.1
        del_eps_lt = eps_lt_tol * self.eps_lt
        out_vec_2 = linstep(liq, self.eps_lt, self.eps_lt+del_eps_lt, 0.0, 1.0)
        out_vec *= out_vec_2
        # TODO: needed to soften the transition to nonzero k_evap on
        # only one side of eps_lt

        # for k, v in enumerate(liq):
        #     if v < self.eps_lt:
        #         out_vec[k] = 0.0

        # TODO: Dividing by 1,000, i.e., using k_bar = 10/1000 = 0.01 for the
        # evaporation rate while still using 10.0 for the condensation rate
        # gives good agreement but why is this necessary at all?
        method = 3

        if method == 1:
            k_evap_scaling = 0.05
            k_evap_new = k_evap_scaling * out_vec

        elif method == 2:
            k_evap_old = self.k_evap.vector()[:]

            relaxation_time = 10.0
            delta_k_evap = 1.0/relaxation_time * (out_vec - k_evap_old)

            k_evap_new = k_evap_old + delta_k_evap

            # k_evap_new *= 0.4

        elif method == 3:
            k_evap_new = out_vec

        self.k_evap.vector()[:] = k_evap_new

    # ================================================================

    def _save_solution(self, tt):
        if not hasattr(self, "xdmf_file"):

            self.xdmf_file = XDMFFile("output/solution.xdmf")
            self.xdmf_file.parameters["flush_output"] = True
            self.xdmf_file.parameters["functions_share_mesh"] = True
            self.xdmf_file.parameters["rewrite_function_mesh"] = False

            self._u = self.u.split()

            self._u[self.c_s_id].rename("steam", "steam")
            self._u[self.eps_l_id].rename("liquid", "liquid")
            self._u[self.T_id].rename("temperature", "temperature")

            self.f_x = Function(self.S)
            self.f_x.rename("xylan", "xylan")

            self.c_xo = Function(self.S)
            self.c_xo.rename("xylooligomer", "xylooligomer")

            self.c_xy = Function(self.S)
            self.c_xy.rename("xylose", "xylose")

            self.c_f = Function(self.S)
            self.c_f.rename("furfural", "furfural")

            self.c_acid_save = Function(self.S)
            self.c_acid_save.rename("acid", "acid")

        f_x_new = project(
            self.u[self.fstar_x_id] / (1.0 - self.eps_p), self.S, solver_type="cg"
        )
        c_xo_new = project(
            self.u[self.cstar_xo_id] / self.u[self.eps_l_id], self.S, solver_type="cg"
        )
        c_xy_new = project(
            self.u[self.cstar_xy_id] / self.u[self.eps_l_id], self.S, solver_type="cg"
        )
        c_f_new = project(
            self.u[self.cstar_f_id] / self.u[self.eps_l_id], self.S, solver_type="cg"
        )
        c_acid_save_new = project(self.c_acid, self.S, solver_type="cg")

        self.f_x.assign(f_x_new)
        self.c_xo.assign(c_xo_new)
        self.c_xy.assign(c_xy_new)
        self.c_f.assign(c_f_new)
        self.c_acid_save.assign(c_acid_save_new)

        self.xdmf_file.write(self._u[self.c_s_id], tt)
        self.xdmf_file.write(self._u[self.eps_l_id], tt)
        self.xdmf_file.write(self.f_x, tt)
        self.xdmf_file.write(self.c_xo, tt)
        self.xdmf_file.write(self.c_xy, tt)
        self.xdmf_file.write(self.c_f, tt)
        self.xdmf_file.write(self._u[self.T_id], tt)
        self.xdmf_file.write(self.c_acid_save, tt)

        eps_l_bar = assemble(self.u[self.eps_l_id]*dx)
        eps_p_bar = assemble((1.0 - self.eps_p)*dx)

        fstar_x_bar = assemble(self.u[self.fstar_x_id]*dx)/eps_p_bar
        cstar_xo_bar = assemble(self.u[self.cstar_xo_id]*dx)/eps_l_bar
        cstar_xy_bar = assemble(self.u[self.cstar_xy_id]*dx)/eps_l_bar
        cstar_f_bar = assemble(self.u[self.cstar_f_id]*dx)/eps_l_bar

        fis = self.rho_s*eps_p_bar/(self.rho_s*eps_p_bar + self.rho_l*eps_l_bar)

        print(f"Time {tt} of {self.t_final}.")
        print(f"| Xylan (g/g):    {fstar_x_bar}")
        print(f"| Xylog (g/L):    {cstar_xo_bar*self.M_xo}")
        print(f"| Xylose (g/L):   {cstar_xy_bar*self.M_xy}")
        print(f"| Furfural (g/L): {cstar_f_bar*self.M_f}")
        print(f"| FIS (g/g):      {fis}")
        print()

    # ================================================================
