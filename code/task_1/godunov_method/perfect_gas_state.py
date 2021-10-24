"""
Файл с начальными данными идеального газа
"""
from consts import R
# R = 8.314  # Дж / кг / K


class PerfectGasState:
    """
    Класс состояния одномерного идеального газа
    """

    def __init__(self, pressure: float, temperature: float, gamma: float, mol_mass: float):
        """
        Конструктор класса

        :param pressure: давление в паскалях
        :type: pressure: float

        :param temperature: температура в кельвинах
        :type temperature: float

        :param gamma: показатель адиабаты
        :type gamma: float

        :param mol_mass: молярная масса
        :type mol_mass: float
        """
        self.pressure = pressure
        self.temperature = temperature
        self.density = mol_mass * pressure / R / temperature
        self.gamma = gamma
        self.mol_mass = mol_mass

    @staticmethod
    def build_from_pressure_density(
            pressure: float,
            density: float,
            gamma: float,
            mol_mass: float
    ) -> "PerfectGasState":
        """
        Конструктор класса

        :param pressure: давление в паскалях
        :type: pressure: float

        :param dencity: плотность в кг / м^3
        :type dencity: float

        :param gamma: показатель адиабаты
        :type gamma: float

        :param mol_mass: молярная масса
        :type mol_mass: float
        """
        temperature = pressure / density * mol_mass / R
        return PerfectGasState(pressure, temperature, gamma, mol_mass)
    
    @property
    def c_v(self) -> float:
        """
        Возвращает удельную теплоемкость при постоянном объеме в Дж / кг / K
        """
        return R / (self.gamma - 1) / self.mol_mass

    @property
    def c_p(self) -> float:
        """
        Возвращает удельную теплоемкость при постоянном давлении ДЖ / кг / K
        """
        return R * self.gamma / (self.gamma - 1) / self.mol_mass

    @property
    def sound_vel(self) -> float:
        """
        Возвращает скорость звука м / с
        """
        return (self.gamma * R * self.temperature / self.mol_mass)**0.5

    @property
    def internal_energy(self) -> float:
        """
        Возвращает удельную внутреннюю энергию газа Дж / кг
        """
        return self.c_v * self.temperature

    @property
    def enthalpy(self) -> float:
        """
        Возвращает удельную энтальпию газа Дж / кг
        """
        return self.c_p * self.temperature
