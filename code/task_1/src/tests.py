# import numpy as np
# from utils import convert_conserved_to_primitive
# from euler_model import Euler
# from solver_eno import SolverEno
# from simple_solver import SimpleSolver
#
# import pickle
# import matplotlib.animation
# import matplotlib.pyplot as plt
#
# w_left = np.array([1., 0.75, 1.])
# w_right = np.array([0.125, 0.0, 0.1])
#
# x0 = 0.5
# t_final = 0.2
#
# position_nodes_number = 300
#
# euler_model = Euler(initial_conditions_=[w_left, w_right, x0], wave_speed_estimator_='pressure_based')
# solver = SolverEno(model=euler_model, cfl=0.15)
# simple_solver = SimpleSolver(model=euler_model, cfl=0.15)
#
# u, space, t_final = solver.solve_problem(position_nodes_number=position_nodes_number,
#                                          t_final=t_final)
#
# # u_simple, space_simple, t_final_simple = simple_solver.solve_problem(position_nodes_number=position_nodes_number,
# #                                                                      t_final=t_final)
#
# time_nodes_number = np.shape(u)[0]
#
# v = np.zeros((time_nodes_number, position_nodes_number, 3))
# v_simple = np.zeros((time_nodes_number, position_nodes_number, 3))
#
# for t in range(time_nodes_number):
#     for i in range(position_nodes_number):
#         v[t, i, :] = convert_conserved_to_primitive(u[t, i, :])
# #         v_simple[t, i, :] = convert_conserved_to_primitive(u_simple[t, i, :])
#
#
# def full_update(num, solution, xs, ax_list):
#     """
#     Функция для построения анимации
#     """
#     y_lim_list = [[-0.05, 1.05], [-2.1, 2.1], [-0.05, 0.45]]
#     labels = ['Точное решение', 'HLLC', 'HLLC + ENO']
#
#     for ax_num, ax in enumerate(ax_list):
#         ax.clear()
#         ax_list[0].set_ylim(y_lim_list[0])
#         ax_list[1].set_ylim(y_lim_list[1])
#         ax_list[2].set_ylim(y_lim_list[2])
#         # for u_num, u_sol in enumerate(solution):
#         #     ax.plot(xs, [u_sol[num][j][ax_num] for j in range(len(xs))], marker='o', label=labels[u_num])
#         ax.plot(xs, [solution[1][num][j][ax_num] for j in range(len(xs))], marker='o', label=labels[1], markersize=2)
#         ax.plot(xs, [solution[2][num][j][ax_num] for j in range(len(xs))], marker='o', label=labels[2], markersize=2)
#         ax.plot(xs, [solution[0][ax_num][num][j] for j in range(len(xs))], marker='o', label=labels[0], markersize=2)
#         ax.legend()
#         ax.grid()
#
#     ax_list[0].set_title('Плотность')
#     ax_list[1].set_title('Скорость')
#     ax_list[2].set_title('Давление')
#
#
# fig = plt.figure(figsize=(14, 10))
# fig.suptitle('Тест 2, T_final =' + str(0.15), fontsize=16)
# ax_1 = fig.add_subplot(2, 2, 1)
# ax_2 = fig.add_subplot(2, 2, 2)
# ax_3 = fig.add_subplot(2, 2, 3)
#
# with open('data_300_826.pickle', 'rb') as f:
#     data_new = pickle.load(f)
#
# ani = matplotlib.animation.FuncAnimation(fig, full_update, time_nodes_number,
#                                          fargs=([data_new, v_simple, v], space,
#                                                 [ax_1, ax_2, ax_3]),
#                                          interval=1, blit=False)
# # ani.save('Тест 2.gif', fps=40)
# plt.show()
#
# # with open('data_300.pickle', 'rb') as f:
# #     data_new = pickle.load(f)
#
# fig = plt.figure(figsize=(14, 10))
# plot_ax_1 = fig.add_subplot(2, 2, 1)
# plot_ax_2 = fig.add_subplot(2, 2, 2)
# plot_ax_3 = fig.add_subplot(2, 2, 3)
#
# plot_ax_1.scatter(space, v[-1][:, 0], s=5, c='blue', label='HLLC + ENO')
# # plot_ax_1.scatter(space, v_simple[-1][:, 0], s=5, c='green', label='HLLC')
# # plot_ax_1.plot(space, data_new[0][-1], color='r', label='Точное решение')
# plot_ax_1.legend()
# plot_ax_2.scatter(space, v[-1][:, 1], s=5, c='blue', label='HLLC + ENO')
# # plot_ax_2.scatter(space, v_simple[-1][:, 1], s=5, c='green', label='HLLC')
# # plot_ax_2.plot(space, data_new[1][-1], color='r', label='Точное решение')
# plot_ax_2.legend()
# plot_ax_3.scatter(space, v[-1][:, 2], s=5, c='blue', label='HLLC + ENO')
# # plot_ax_3.scatter(space, v_simple[-1][:, 2], s=5, c='green', label='HLLC')
# # plot_ax_3.plot(space, data_new[2][-1], color='r', label='Точное решение')
# plot_ax_3.legend()
#
# plot_ax_1.set_title('Плотность')
# plot_ax_2.set_title('Скорость')
# plot_ax_3.set_title('Давление')
#
# plot_ax_1.grid()
# plot_ax_2.grid()
# plot_ax_3.grid()
#
# # plt.savefig('Тест 2.png')
# plt.show()
