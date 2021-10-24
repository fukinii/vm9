"""
Файл с функцией рассчета давления за ударной волной
"""
import abc
import numpy as np
from perfect_gas_state import PerfectGasState
from reflection_func import reflection_func, derivative_reflection_func
from shock_func import shock_func, derivative_shock_func


class InitialApproximation(abc.ABC):
    """
    Класс начальныых приближений
    """

    @staticmethod
    @abc.abstractmethod
    def calc_initial_approximation(
            left_gas: PerfectGasState,
            right_gas: PerfectGasState,
            left_velocity: float,
            right_velocity: float
    ) -> float:
        """
        Класс для вычисления начального приближения
        """


class InitialApproximationMean(InitialApproximation):
    """
    Считает начальное приближение средним арифметическм
    """

    @staticmethod
    def calc_initial_approximation(
            left_gas: PerfectGasState,
            right_gas: PerfectGasState,
            left_velocity: float,
            right_velocity: float
    ) -> float:
        """
        Класс для вычисления начального приближения
        """

        return 0.5 * (left_gas.pressure + right_gas.pressure)


def calc_function(gas: PerfectGasState, pressure: float) -> float:
    """
    Вычиляет функцию в зависимости от давления
    """
    if pressure > gas.pressure:
        return shock_func(gas, pressure)
    else:
        return reflection_func(gas, pressure)


def calc_derivative(gas: PerfectGasState, pressure: float) -> float:
    """
    Вычисляет производную
    """
    if pressure > gas.pressure:
        return derivative_shock_func(gas, pressure)
    else:
        return derivative_reflection_func(gas, pressure)


def pressure_solve(
        initial_approximation: InitialApproximation,
        left_gas: PerfectGasState,
        right_gas: PerfectGasState,
        left_velocity: float,
        right_velocity: float,
        tolerance: float = 1e-6
) -> float:
    """
    Вычиляет p* с заданной точностью
    """
    pressure = initial_approximation.calc_initial_approximation(
        left_gas,
        right_gas,
        left_velocity,
        right_velocity
    )
    error = 10
    while error > tolerance:
        delta_velocity = right_velocity - left_velocity
        function = calc_function(left_gas, pressure) + calc_function(right_gas, pressure) + delta_velocity
        derivative = calc_derivative(left_gas, pressure) + calc_derivative(right_gas, pressure)
        delta = - function / derivative
        new_pressure = pressure + delta
        error = abs(delta) / (abs(pressure) + abs(new_pressure))
        pressure = max(pressure + delta, tolerance)
    return pressure


def shock_wave_density(
        gas: PerfectGasState,
        pressure_star: float,
        u_star: float,
        velocity: float,
        x_coord: float,
        time: float
) -> float:
    density = gas.density * (
            pressure_star / gas.pressure + (gas.gamma - 1) / (gas.gamma + 1)
    ) / (
                      (gas.gamma - 1) / (gas.gamma + 1) * pressure_star / gas.pressure + 1
              )
    return density


def s_shock_wave_velocity(
        gas: PerfectGasState,
        pressure_star: float,
        u_star: float,
        velocity: float,
        x_coord: float,
        time: float
) -> float:
    left_velocity = velocity
    velocity_s = left_velocity - gas.sound_vel * (
            (
                    gas.gamma + 1
            ) / 2 / gas.gamma * pressure_star / gas.pressure + (
                    gas.gamma - 1
            ) / 2 / gas.gamma
    ) ** 0.5
    return velocity_s


def shock_wave_velocity(
        gas: PerfectGasState,
        pressure_star: float,
        u_star: float,
        velocity: float,
        x_coord: float,
        time: float
) -> float:
    return u_star


def shock_wave_pressure(
        gas: PerfectGasState,
        pressure_star: float,
        u_star: float,
        velocity: float,
        x_coord: float,
        time: float
) -> float:
    return pressure_star


def left_reflection_wave_density(
        gas: PerfectGasState,
        pressure_star: float,
        u_star: float,
        velocity: float,
        x_coord: float,
        time: float
) -> float:
    power = 2. / (gas.gamma - 1.)
    res = gas.density * (
            2. / (gas.gamma + 1) + (gas.gamma - 1) / (gas.gamma + 1) / gas.sound_vel * (velocity - x_coord / time)
    ) ** power
    return res


def left_reflection_wave_velocity(
        gas: PerfectGasState,
        pressure_star: float,
        u_star: float,
        velocity: float,
        x_coord: float,
        time: float
) -> float:
    res = 2. / (gas.gamma + 1.) * (
            gas.sound_vel + (gas.gamma - 1.) / 2. * velocity + x_coord / time
    )
    return res


def left_reflection_wave_pressure(
        gas: PerfectGasState,
        pressure_star: float,
        u_star: float,
        velocity: float,
        x_coord: float,
        time: float
) -> float:
    power = 2. * gas.gamma / (gas.gamma - 1)
    res = gas.pressure * (
            2. / (gas.gamma + 1) + (gas.gamma - 1) / (gas.gamma + 1) / gas.sound_vel * (
            velocity - x_coord / time
    )
    ) ** power
    return res


def right_reflection_wave_density(
        gas: PerfectGasState,
        pressure_star: float,
        u_star: float,
        velocity: float,
        x_coord: float,
        time: float
) -> float:
    power = 2. / (gas.gamma - 1)
    res = gas.density * (
            2. / (gas.gamma + 1) - (gas.gamma - 1) / (gas.gamma + 1) / gas.sound_vel * (
            velocity - x_coord / time
    )
    ) ** power
    return res


def right_reflection_wave_velocity(
        gas: PerfectGasState,
        pressure_star: float,
        u_star: float,
        velocity: float,
        x_coord: float,
        time: float
) -> float:
    res = 2. / (gas.gamma + 1.) * (
            -gas.sound_vel + (gas.gamma - 1.) / 2. * velocity + x_coord / time
    )
    return res


def right_reflection_wave_pressure(
        gas: PerfectGasState,
        pressure_star: float,
        u_star: float,
        velocity: float,
        x_coord: float,
        time: float
) -> float:
    power = 2. * gas.gamma / (gas.gamma - 1)
    res = gas.pressure * (
            2. / (gas.gamma + 1) - (gas.gamma - 1) / (gas.gamma + 1) / gas.sound_vel * (
            velocity - x_coord / time
    )
    ) ** power
    return res


def left_side_solution(
        gas: PerfectGasState,
        pressure_star: float,
        u_star: float,
        velocity: float,
        x_coord: float,
        time: float,
):
    if pressure_star > gas.pressure:
        velocity_s = velocity - gas.sound_vel * (
                (
                        gas.gamma + 1
                ) / 2 / gas.gamma * pressure_star / gas.pressure + (
                        gas.gamma - 1
                ) / 2 / gas.gamma
        ) ** 0.5
        if x_coord / time < velocity_s:
            return gas.pressure, gas.density, velocity
        else:
            res_pressure = shock_wave_pressure(gas, pressure_star, u_star, velocity, x_coord, time)
            res_density = shock_wave_density(gas, pressure_star, u_star, velocity, x_coord, time)
            res_velocity = shock_wave_velocity(gas, pressure_star, u_star, velocity, x_coord, time)
            return res_pressure, res_density, res_velocity
    else:
        velocity_s_hl = velocity - gas.sound_vel
        power = (gas.gamma - 1.) / 2. / gas.gamma
        velocity_s_tl = u_star - gas.sound_vel * (pressure_star / gas.pressure) ** power
        if x_coord / time < velocity_s_hl:
            return gas.pressure, gas.density, velocity
        elif velocity_s_hl <= x_coord / time <= velocity_s_tl:
            res_pressure = left_reflection_wave_pressure(gas, pressure_star, u_star, velocity, x_coord, time)
            res_density = left_reflection_wave_density(gas, pressure_star, u_star, velocity, x_coord, time)
            res_velocity = left_reflection_wave_velocity(gas, pressure_star, u_star, velocity, x_coord, time)
            return res_pressure, res_density, res_velocity
        else:
            res_density = gas.density * (pressure_star / gas.pressure) ** (1. / gas.gamma)
            return pressure_star, res_density, u_star


def right_side_solution(
        gas: PerfectGasState,
        pressure_star: float,
        u_star: float,
        velocity: float,
        x_coord: float,
        time: float,
):
    if pressure_star > gas.pressure:
        velocity_s = velocity + gas.sound_vel * (
                (
                        gas.gamma + 1
                ) / 2 / gas.gamma * pressure_star / gas.pressure + (
                        gas.gamma - 1
                ) / 2 / gas.gamma
        ) ** 0.5
        if x_coord / time >= velocity_s:
            return gas.pressure, gas.density, velocity
        else:
            res_pressure = shock_wave_pressure(gas, pressure_star, u_star, velocity, x_coord, time)
            res_density = shock_wave_density(gas, pressure_star, u_star, velocity, x_coord, time)
            res_velocity = shock_wave_velocity(gas, pressure_star, u_star, velocity, x_coord, time)
            return res_pressure, res_density, res_velocity
    else:
        velocity_s_hr = velocity + gas.sound_vel
        power = (gas.gamma - 1.) / 2. / gas.gamma
        velocity_s_tr = u_star + gas.sound_vel * (pressure_star / gas.pressure) ** power
        if x_coord / time > velocity_s_hr:
            return gas.pressure, gas.density, velocity
        elif velocity_s_tr <= x_coord / time <= velocity_s_hr:
            res_pressure = right_reflection_wave_pressure(gas, pressure_star, u_star, velocity, x_coord, time)
            res_density = right_reflection_wave_density(gas, pressure_star, u_star, velocity, x_coord, time)
            res_velocity = right_reflection_wave_velocity(gas, pressure_star, u_star, velocity, x_coord, time)
            return res_pressure, res_density, res_velocity
        else:
            res_density = gas.density * (pressure_star / gas.pressure) ** (1. / gas.gamma)
            return pressure_star, res_density, u_star


def solver_g(
        x_mesh: np.ndarray,
        times: np.ndarray,
        initial_approximation: InitialApproximation,
        left_gas: PerfectGasState,
        right_gas: PerfectGasState,
        left_velocity: float,
        right_velocity: float,
        x_0: float = 0.
) -> np.ndarray:
    """
    Произволдит расчет в сетке x_mesh во времена times
    """
    pressure_star = pressure_solve(
        initial_approximation,
        left_gas,
        right_gas,
        left_velocity,
        right_velocity
    )
    u_star = 0.5 * (left_velocity + right_velocity) + 0.5 * (
            calc_function(right_gas, pressure_star) - calc_function(left_gas, pressure_star)
    )
    density = np.ndarray((len(times), len(x_mesh)))
    pressure = np.ndarray((len(times), len(x_mesh)))
    velocity = np.ndarray((len(times), len(x_mesh)))
    greater_zero = x_mesh > x_0
    less_zero = np.logical_not(greater_zero)
    density[0, greater_zero] = right_gas.density
    pressure[0, greater_zero] = right_gas.pressure
    velocity[0, greater_zero] = right_velocity
    density[0, less_zero] = left_gas.density
    pressure[0, less_zero] = left_gas.pressure
    velocity[0, less_zero] = left_velocity
    for time_index, time in enumerate(times):
        if time_index > 0:
            for x_index, x_coord in enumerate(x_mesh):
                x_coord = x_coord - x_0
                if x_coord / time > u_star:
                    (
                        pressure[time_index, x_index],
                        density[time_index, x_index],
                        velocity[time_index, x_index]
                    ) = right_side_solution(right_gas, pressure_star, u_star, right_velocity, x_coord, time)
                else:
                    (
                        pressure[time_index, x_index],
                        density[time_index, x_index],
                        velocity[time_index, x_index]
                    ) = left_side_solution(left_gas, pressure_star, u_star, left_velocity, x_coord, time)
    return pressure, density, velocity
