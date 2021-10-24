"""
Файл с функцией левой ударной волны f_l
"""
from perfect_gas_state import PerfectGasState


def calc_a(gas: PerfectGasState) -> float:
    """
    Считает коэффициент А
    """
    return 2 / (gas.gamma + 1) / gas.density


def calc_b(gas: PerfectGasState) -> float:
    """
    Считает коэффициент B
    """
    return (gas.gamma - 1) / (gas.gamma + 1) * gas.pressure


def shock_func(gas: PerfectGasState, pressure: float) -> float:
    """
    Функция, считающая функцию для ударной волны
    """
    return (pressure - gas.pressure) * (calc_a(gas) / (calc_b(gas) + pressure))**0.5

def derivative_shock_func(gas: PerfectGasState, pressure: float) -> float:
    """
    Производная функция дял ударной волны
    """
    b_coef_plas_p = calc_b(gas) + pressure
    return (calc_a(gas) / b_coef_plas_p)**0.5 * (1 - (pressure - gas.pressure) / 2 / b_coef_plas_p)
