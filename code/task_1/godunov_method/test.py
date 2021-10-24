"""
Файл  запуском кода
"""
from solver import pressure_solve, InitialApproximationMean, solver_g
from perfect_gas_state import PerfectGasState
import numpy as np

gamma = 1.4
mol_mass = 0.029

left_pressure = 1.
left_density = 1.
left_velocity = 0.75

left_gas = PerfectGasState.build_from_pressure_density(left_pressure, left_density, gamma, mol_mass)

right_pressure = 0.1
right_density = 0.125
right_velocity = 0.0

right_gas = PerfectGasState.build_from_pressure_density(right_pressure, right_density, gamma, mol_mass)

# print(pressure_solve(InitialApproximationMean, left_gas, right_gas, left_velocity, right_velocity))

x_mesh = np.linspace(0, 1, 100)

times = np.linspace(0, 0.2, 800)

pressure, density, velocity = solver_g(
    x_mesh,
    times,
    InitialApproximationMean,
    left_gas,
    right_gas,
    left_velocity,
    right_velocity,
    x_0=0.5
)

import matplotlib.pyplot as plt

fig = plt.figure(figsize=(20, 10))

fig.suptitle('Две ударные волны', fontsize=28)
ax_1 = fig.add_subplot(311)
ax_2 = fig.add_subplot(312)
ax_3 = fig.add_subplot(313)


def plot(index):
    p_max = np.max(pressure)
    p_min = np.min(pressure)
    density_max = np.max(density)
    density_min = np.min(density)
    velocity_max = np.max(velocity)
    velocity_min = np.min(velocity)
    ax_1.clear()
    ax_2.clear()
    ax_3.clear()

    # ax_1.set_ylim(p_min, p_max)

    ax_1.set_ylabel("Давление", fontsize=24)
    # ax_2.set_ylim(density_min, density_max)
    ax_2.set_ylabel("Плотность", fontsize=24)
    # ax_3.set_ylim(velocity_min, velocity_max)
    ax_3.set_ylabel("Скорость", fontsize=24)

    ax_1.plot(x_mesh, pressure[index], linewidth=3.0)
    ax_2.plot(x_mesh, density[index], linewidth=3.0)
    ax_3.plot(x_mesh, velocity[index], linewidth=3.0)

import pickle
from matplotlib import animation

anim = animation.FuncAnimation(fig, plot, interval=200, frames=len(times) - 1)

plt.show()

fig = plt.figure(figsize=(14, 10))
plot_ax_1 = fig.add_subplot(2, 2, 1)
plot_ax_2 = fig.add_subplot(2, 2, 2)
plot_ax_3 = fig.add_subplot(2, 2, 3)

plot_ax_1.plot(x_mesh, density[-1])
plot_ax_2.plot(x_mesh, velocity[-1])
plot_ax_3.plot(x_mesh, pressure[-1])

plot_ax_1.set_title('Плотность')
plot_ax_2.set_title('Скорость')
plot_ax_3.set_title('Давление')

plot_ax_1.grid()
plot_ax_2.grid()
plot_ax_3.grid()

plt.show()
data = [density, velocity, pressure]
with open('data.pickle', 'wb') as f:
    pickle.dump(data, f)
# Writer = animation.writers['ffmpeg']
# writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)

# anim.save("Две ударные волны.mp4", writer = writer)
