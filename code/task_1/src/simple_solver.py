import numpy as np
from src.euler_model import Euler
from src.utils import convert_conserved_to_primitive, calc_speed_of_sound, calc_divided_difference
from src.solver import Solver


class SimpleSolver(Solver):

    def reconstruct(self, u):
        return self.make_reconstruction_1p(u)

    def build_boundary_conditions(self, u, ulr):
        return self.build_boundary_conditions_1p(u, ulr)
