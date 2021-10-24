import numpy as np
from euler_model import Euler
from utils import convert_conserved_to_primitive, calc_speed_of_sound, calc_divided_difference


class SolverEno():
    def __init__(self, model: Euler, cfl=0.5):
        self.model = model
        self.cfl = cfl
        crj_table_3 = {-1: np.array([11 / 6, -7 / 6, 1 / 3]),
                       0: np.array([1 / 3, 5 / 6, -1 / 6]),
                       1: np.array([-1 / 6, 5 / 6, 1 / 3]),
                       2: np.array([1 / 3, -7 / 6, 11 / 6])}

        crj_table_5 = {-1: np.array([137 / 60, -163 / 60, 137 / 60, -21 / 20, 1 / 5]),
                       0: np.array([1 / 5, 77 / 60, -43 / 60, 17 / 60, -1 / 20]),
                       1: np.array([-1 / 20, 9 / 20, 47 / 60, -13 / 60, 1 / 30]),
                       2: np.array([1 / 30, -13 / 60, 47 / 60, 9 / 20, -1 / 20]),
                       3: np.array([-1 / 20, 17 / 60, -43 / 60, 77 / 60, 1 / 5]),
                       4: np.array([1 / 5, -21 / 20, 137 / 60, -163 / 60, 137 / 60]), }

        self.crj_table = {3: crj_table_3, 5: crj_table_5}

    def build_boundary_conditions_1p(self, u, ulr):
        if self.model.boundary_condition == 'transparent':
            u[0] = u[1]
            u[-1] = u[-2]
            ulr[0, :, 0] = u[0]
            ulr[-1, :, 1] = u[-1]
        indent = 0
        return u, ulr, indent

    @staticmethod
    def make_reconstruction_1p(u):
        return u[:], u[:]

    @staticmethod
    def calc_rhs(flux_left, flux_right, h):
        return -(flux_right - flux_left) / h

    def calc_delta_time(self, u, h):
        min_value = 1e9
        for i in range(len(u)):
            v_i = convert_conserved_to_primitive(u[i])
            velocity = v_i[1]
            down = calc_speed_of_sound(u[i]) + np.abs(velocity)
            min_value = min(min_value,
                            h / down)
        return self.cfl * min_value

    def build_boundary_conditions_3p(self, u, ulr):
        if self.model.boundary_condition == 'transparent':
            u[0] = u[1]
            u[-1] = u[-2]
            ulr[0, :, 0] = u[0]
            ulr[0, :, 1] = u[0]
            ulr[1, :, 0] = u[0]
            ulr[-1, :, 1] = u[-1]
            ulr[-1, :, 0] = u[-1]
            ulr[-2, :, 1] = u[-1]

        indent = 1
        return u, ulr, indent

    def build_boundary_conditions_5p(self, u, ulr):
        if self.model.boundary_condition == 'transparent':
            u[0] = u[2]
            u[1] = u[2]
            u[-1] = u[-3]
            u[-2] = u[-3]
            ulr[0, :, 0] = u[0]
            ulr[0, :, 1] = u[0]
            ulr[1, :, 0] = u[0]
            ulr[1, :, 1] = u[1]
            ulr[2, :, 0] = u[1]
            ulr[-1, :, 1] = u[-1]
            ulr[-1, :, 0] = u[-1]
            ulr[-2, :, 1] = u[-1]
            ulr[-2, :, 0] = u[-2]
            ulr[-3, :, 1] = u[-2]

        indent = 2
        return u, ulr, indent

    @staticmethod
    def calc_left_indent(u, position, cell_index, k=3):

        i_s = cell_index
        for m in range(2, k + 1):
            v_left = calc_divided_difference(u, i_s - 1.5, i_s - 1.5 + m, position)
            v_right = calc_divided_difference(u, i_s - 0.5, i_s - 0.5 + m, position)
            if np.abs(v_left) < np.abs(v_right):
                i_s -= 1

        return cell_index - i_s

    def make_reconstruction_for_component(self, u, position, cell_index, order):

        k = order
        most_left_point_in_the_stencil = self.calc_left_indent(u, position, cell_index, k)
        crj_coefficients = self.crj_table[k][most_left_point_in_the_stencil]
        crj_tilda_coefficients = self.crj_table[k][most_left_point_in_the_stencil - 1]

        left = 0
        right = 0

        for j in range(k):
            right += crj_coefficients[j] * u[cell_index - most_left_point_in_the_stencil + j]
            left += crj_tilda_coefficients[j] * u[cell_index - most_left_point_in_the_stencil + j]

        return left, right

    def make_reconstruction(self, u_data, position, cell_index, order=3):

        left = np.zeros(3)
        right = np.zeros_like(left)

        for i in range(3):
            left[i], right[i] = self.make_reconstruction_for_component(u_data[:, i], position, cell_index, order)

        return left, right

    def solve_problem(self, reconstructor, position_nodes_number, t_final, xl=0., xr=1.):
        u = np.zeros((1, position_nodes_number, 3))

        flux = np.zeros((position_nodes_number + 1, 3))
        ulr = np.zeros((position_nodes_number + 1, 3, 2))

        rhs = np.zeros((position_nodes_number, 3))
        h = (xr - xl) / position_nodes_number

        position = np.linspace(xl + h / 2, xr - h / 2, position_nodes_number)

        u[0] = self.model.build_initial_conditions(position)
        time = 0
        t = 0
        while time < t_final:

            u[t], ulr, indent = self.build_boundary_conditions_5p(u[t], ulr)
            delta_time = self.calc_delta_time(u[t], h)

            for cell_index in range(indent - 1, np.shape(position)[0] - indent):
                ulr[cell_index, :, 1], ulr[cell_index + 1, :, 0] = self.make_reconstruction(u[t], position,
                                                                                            cell_index, order=3)

                flux[cell_index, :] = self.model.flux_riemann_solver(ulr[cell_index, :, 0],
                                                                     ulr[cell_index, :, 1])
                if cell_index > 1:
                    rhs[cell_index - 1, :] = self.calc_rhs(flux[cell_index - 1, :], flux[cell_index, :], h)
            u_next = u[t] + delta_time * rhs
            u_next = np.reshape(u_next, (1, len(u_next), 3))
            u = np.concatenate((u, u_next), axis=0)
            t += 1
            time += delta_time
            print(time)
        return u, position, time

    def solve_problem_simple(self, position_nodes_number, t_final, xl=0., xr=1.):
        u = np.zeros((1, position_nodes_number, 3))

        flux = np.zeros((position_nodes_number + 1, 3))
        ulr = np.zeros((position_nodes_number + 1, 3, 2))

        rhs = np.zeros((position_nodes_number, 3))
        h = (xr - xl) / position_nodes_number

        position = np.linspace(xl + h / 2, xr - h / 2, position_nodes_number)

        u[0] = self.model.build_initial_conditions(position)
        time = 0
        t = 0
        while time < t_final:

            u[t], ulr, indent = self.build_boundary_conditions_1p(u[t], ulr)
            delta_time = self.calc_delta_time(u[t], h)

            for cell_index in range(indent, np.shape(position)[0] - indent):
                ulr[cell_index, :, 1], ulr[cell_index + 1, :, 0] = self.make_reconstruction_1p(u[t][cell_index, :])

                flux[cell_index, :] = self.model.flux_riemann_solver(ulr[cell_index, :, 0],
                                                                     ulr[cell_index, :, 1])
                if cell_index > 1:
                    rhs[cell_index - 1, :] = self.calc_rhs(flux[cell_index - 1, :], flux[cell_index, :], h)

            u_next = u[t] + delta_time * rhs
            u_next = np.reshape(u_next, (1, len(u_next), 3))
            u = np.concatenate((u, u_next), axis=0)
            t += 1
            time += delta_time
            print(time)
        return u, position, time
