"""
Файл с функциями для волны разрежения
"""
from perfect_gas_state import PerfectGasState


def reflection_func(gas: PerfectGasState, pressure: float) -> float:
    """
    Рассчитывает функцию для волны разрежения
    """
    power = (gas.gamma - 1) / 2 / gas.gamma
    return 2 * gas.sound_vel / (gas.gamma - 1) * ((pressure / gas.pressure)**power - 1)


def derivative_reflection_func(gas: PerfectGasState, pressure: float) -> float:
    """
    Вовзращает производную от функции для волны разрежения
    """
    power = -(gas.gamma + 1) / 2 / gas.gamma
    return (pressure / gas.pressure)**power / gas.sound_vel / gas.density
