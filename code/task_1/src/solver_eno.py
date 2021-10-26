import numpy as np
from src.euler_model import Euler
from src.utils import calc_divided_difference
from src.solver_bad import Solver


class SolverEno(Solver):
    def __init__(self, model: Euler, cfl=0.5):
        super().__init__(model, cfl)
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

    def reconstruct(self, u, position, cell_index, order):
        return self.make_reconstruction(u, position, cell_index, order)

    def build_boundary_conditions(self, u, ulr):
        return self.build_boundary_conditions_5p(u, ulr)
