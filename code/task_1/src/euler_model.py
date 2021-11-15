import numpy as np
from typing import List
from src.utils import convert_primitive_to_conserved, convert_conserved_to_primitive, calc_flux_x, calc_speed_of_sound, \
    calc_flux_y


class Euler:
    def __init__(
            self,
            initial_conditions_: List,
            boundary_condition_: str = 'transparent',
            riemann_solver_: str = 'hllc',
            wave_speed_estimator_: str = 'pressure_based',
            gamma_: np.double = 1.4,
            dim=1
    ):

        self.gamma: np.double = gamma_
        self.riemann_solver: str = riemann_solver_
        self.wave_speed_estimator: str = wave_speed_estimator_
        self.boundary_condition: str = boundary_condition_
        self.initial_conditions: List = initial_conditions_
        self.dim_sol = dim + 2

    def build_initial_conditions(self, x, y):
        u0 = np.zeros((len(x), len(y), self.dim_sol))

        """2 обл"""
        a = x < self.initial_conditions[4]
        b = y < self.initial_conditions[5]
        c2 = np.einsum('i,j', a, b)

        """4 обл"""
        a = x >= self.initial_conditions[4]
        b = y >= self.initial_conditions[5]
        c4 = np.einsum('i,j', a, b)

        """2 обл"""
        a = x < self.initial_conditions[4]
        b = y >= self.initial_conditions[5]
        c1 = np.einsum('i,j', a, b)

        """2 обл"""
        a = x >= self.initial_conditions[4]
        b = y < self.initial_conditions[5]
        c3 = np.einsum('i,j', a, b)

        # mask = x >= self.initial_conditions[4] and y >= self.initial_conditions[5]
        u0[c1, :] = convert_primitive_to_conserved(self.initial_conditions[0], self.gamma)

        # mask = x < self.initial_conditions[4] and y >= self.initial_conditions[5]
        u0[c2, :] = convert_primitive_to_conserved(self.initial_conditions[1], self.gamma)

        # mask = x < self.initial_conditions[4] and y < self.initial_conditions[5]
        u0[c3, :] = convert_primitive_to_conserved(self.initial_conditions[2], self.gamma)

        # mask = x >= self.initial_conditions[4] and y < self.initial_conditions[5]
        u0[c4, :] = convert_primitive_to_conserved(self.initial_conditions[3], self.gamma)

        return u0

    def calc_qk(self, p, p_star):
        # print(p, p_star)
        return np.sqrt(1 + (self.gamma + 1) * (p_star / p - 1) / 2 / self.gamma) if p_star <= p else 1

    def estimate_pressure_based_speed_of_wave_x(self, ul, ur):
        al = calc_speed_of_sound(ul)
        ar = calc_speed_of_sound(ur)
        vl = convert_conserved_to_primitive(ul)
        vr = convert_conserved_to_primitive(ur)
        aver_rho = (vl[0] + vr[0]) / 2
        aver_a = (al + ar) / 2
        p_pvrs = (vl[3] + vr[3] - (vr[1] - vl[1]) * aver_rho * aver_a) / 2
        p_star = max(0, p_pvrs)
        ql = self.calc_qk(vl[3], p_star)
        qr = self.calc_qk(vr[3], p_star)
        return vl[1] - al * ql, vr[1] + ar * qr

    def estimate_speed_of_wave_by_einfeldt_x(self, ul, ur):
        vl = convert_conserved_to_primitive(ul)  # V = [rho, u, p]
        vr = convert_conserved_to_primitive(ur)  # U = [rho, rho*u, E]

        hl = (ul[3] + vl[3]) / vl[0]
        hr = (ur[3] + vr[3]) / vr[0]

        u = (np.sqrt(vl[0]) * vl[1] + np.sqrt(vr[0]) * vr[1]) / (np.sqrt(vl[0]) + np.sqrt(vr[0]))
        h = (np.sqrt(vl[0]) * hl + np.sqrt(vr[0]) * hr) / (np.sqrt(vl[0]) + np.sqrt(vr[0]))
        a = np.sqrt((self.gamma - 1) * (h - 0.5 * u ** 2))
        # a = np.sqrt(np.abs((self.gamma - 1) * (h - 0.5 * u ** 2)))

        sl = u - a
        sr = u + a

        return sl, sr

    def estimate_speed_of_wave_by_einfeldt_y(self, ul, ur):
        vl = convert_conserved_to_primitive(ul)  # V = [rho, u, p]
        vr = convert_conserved_to_primitive(ur)  # U = [rho, rho*u, E]

        hl = (ul[3] + vl[3]) / vl[0]
        hr = (ur[3] + vr[3]) / vr[0]

        u = (np.sqrt(vl[0]) * vl[2] + np.sqrt(vr[0]) * vr[2]) / (np.sqrt(vl[0]) + np.sqrt(vr[0]))
        h = (np.sqrt(vl[0]) * hl + np.sqrt(vr[0]) * hr) / (np.sqrt(vl[0]) + np.sqrt(vr[0]))
        a = np.sqrt((self.gamma - 1) * (h - 0.5 * u ** 2))
        # a = np.sqrt(np.abs((self.gamma - 1) * (h - 0.5 * u ** 2)))

        sl = u - a
        sr = u + a

        return sl, sr

    @staticmethod
    def estimate_speed_of_wave_by_davis_x(ul, ur):
        vl = convert_conserved_to_primitive(ul)  # V = [rho, u, p]
        vr = convert_conserved_to_primitive(ur)

        al = calc_speed_of_sound(ul)
        ar = calc_speed_of_sound(ur)

        sl = min(vl[1] - al, vr[1] - ar)
        sr = max(vl[1] + al, vr[1] + ar)
        return sl, sr

    def estimate_speed_of_wave_x(self, ul, ur):
        if self.wave_speed_estimator == 'pressure_based':
            return self.estimate_pressure_based_speed_of_wave_x(ul, ur)
        if self.wave_speed_estimator == 'Einfeldt':
            return self.estimate_speed_of_wave_by_einfeldt_x(ul, ur)
        if self.wave_speed_estimator == 'Davis':
            return self.estimate_speed_of_wave_by_davis_x(ul, ur)

    @staticmethod
    def estimate_speed_of_wave_by_davis_y(ul, ur):
        vl = convert_conserved_to_primitive(ul)  # V = [rho, u, p]
        vr = convert_conserved_to_primitive(ur)

        al = calc_speed_of_sound(ul)
        ar = calc_speed_of_sound(ur)

        sl = min(vl[2] - al, vr[2] - ar)
        sr = max(vl[2] + al, vr[2] + ar)
        return sl, sr

    def estimate_speed_of_wave_y(self, ul, ur):
        if self.wave_speed_estimator == 'pressure_based':
            pass
            # return self.estimate_pressure_based_speed_of_wave_y(ul, ur)
        if self.wave_speed_estimator == 'Einfeldt':
            return self.estimate_speed_of_wave_by_einfeldt_y(ul, ur)
        if self.wave_speed_estimator == 'Davis':
            return self.estimate_speed_of_wave_by_davis_y(ul, ur)

    @staticmethod
    def calc_s_star_x(vl, vr, sr, sl):

        s_star = (vr[3] - vl[3] + vl[0] * vl[1] * (sl - vl[1]) - vr[0] * vr[1] * (sr - vr[1])) / (
                vl[0] * (sl - vl[1]) - vr[0] * (sr - vr[1]))

        return s_star

    @staticmethod
    def calc_s_star_y(vl, vr, sr, sl):

        s_star = (vr[3] - vl[3] + vl[0] * vl[2] * (sl - vl[2]) - vr[0] * vr[2] * (sr - vr[2])) / (
                vl[0] * (sl - vl[2]) - vr[0] * (sr - vr[2]))

        return s_star

    @staticmethod
    def calc_uk_star_x(uk, vk, s_star, sk):

        uk_star = np.zeros_like(vk)

        coefficient = vk[0] * (sk - vk[1]) / (sk - s_star)
        uk_star[0] = 1
        uk_star[1] = s_star
        uk_star[2] = vk[2]
        uk_star[3] = (uk[3] / vk[0] + (s_star - vk[1]) * (s_star + vk[3] / vk[0] / (sk - vk[1])))
        uk_star = coefficient * uk_star

        return uk_star

    @staticmethod
    def calc_uk_star_y(uk, vk, s_star, sk):

        uk_star = np.zeros_like(vk)

        coefficient = vk[0] * (sk - vk[2]) / (sk - s_star)
        uk_star[0] = 1
        uk_star[1] = s_star
        uk_star[2] = vk[1]
        # uk_star[1] = vk[1]
        # uk_star[2] = s_star
        uk_star[3] = (uk[3] / vk[0] + (s_star - vk[2]) * (s_star + vk[3] / vk[0] / (sk - vk[2])))
        uk_star = coefficient * uk_star

        return uk_star

    def calc_flux_hllc_x(self, ul, ur):  # TODO: Переписать под двумерный случай
        vl = convert_conserved_to_primitive(ul)
        vr = convert_conserved_to_primitive(ur)
        sl, sr = self.estimate_speed_of_wave_x(ul, ur)
        fl = calc_flux_x(ul)
        fr = calc_flux_x(ur)

        s_star = self.calc_s_star_x(vl, vr, sr, sl)

        ur_star = self.calc_uk_star_x(ur, vr, s_star, sr)
        ul_star = self.calc_uk_star_x(ul, vl, s_star, sl)

        if sl > 0:
            f = fl
        elif sr < 0:
            f = fr
        elif s_star > 0:
            f = fl + sl * (ul_star - ul)
        else:
            f = fr + sr * (ur_star - ur)

        return f

    def flux_riemann_solver_x(self, ul, ur):
        if self.riemann_solver == 'hllc':
            return self.calc_flux_hllc_x(ul, ur)

    def calc_flux_hllc_y(self, ul, ur):
        vl = convert_conserved_to_primitive(ul)
        vr = convert_conserved_to_primitive(ur)
        sl, sr = self.estimate_speed_of_wave_y(ul, ur)
        fl = calc_flux_y(ul)
        fr = calc_flux_y(ur)

        s_star = self.calc_s_star_y(vl, vr, sr, sl)

        ur_star = self.calc_uk_star_y(ur, vr, s_star, sr)
        ul_star = self.calc_uk_star_y(ul, vl, s_star, sl)

        if sl > 0:
            f = fl
        elif sr < 0:
            f = fr
        elif s_star > 0:
            f = fl + sl * (ul_star - ul)
        else:
            f = fr + sr * (ur_star - ur)

        a = 1
        # if sl < 0:
        #     f = fl
        # elif sr > 0:
        #     f = fr
        # elif s_star < 0:
        #     f = fl + sl * (ul_star - ul)
        # else:
        #     f = fr + sr * (ur_star - ur)

        return f

    def flux_riemann_solver_y(self, ul, ur):
        if self.riemann_solver == 'hllc':
            return self.calc_flux_hllc_y(ul, ur)
