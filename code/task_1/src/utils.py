import numpy as np


def convert_primitive_to_conserved(w, gamma=1.4):
    u = np.zeros_like(w)
    u[0] = w[0]
    u[1] = w[0] * w[1]
    u[2] = w[0] * w[2]
    u[3] = w[3] / (gamma - 1) + 0.5 * w[0] * (w[1] ** 2 + w[2] ** 2)

    return u


def convert_conserved_to_primitive_zero(u, gamma=1.4):
    v = np.zeros_like(u)
    v[0] = u[0]
    v[1] = u[1] / u[0]
    v[2] = (u[2] - 0.5 * v[0] * v[1] ** 2) * (gamma - 1)

    # assert v[2] > 0, \
    #    'Отрицательное давление'
    # assert v[0] > 0, \
    #     'Отрицательная плотность'
    v[0] = v[0] if v[0] >= 0 else 1e-9
    v[3] = v[3] if v[3] >= 0 else 1e-9

    u_new = convert_primitive_to_conserved(v)
    # if not np.array_equal(u, u_new):
    #     a = 1

    return v, u_new


def convert_conserved_to_primitive(u, gamma=1.4):
    v = np.zeros_like(u)
    v[0] = u[0]
    v[1] = u[1] / u[0]
    v[2] = u[2] / u[0]
    v[3] = (u[2] - 0.5 * v[0] * (v[1] ** 2 + v[2] ** 2)) * (gamma - 1)

    # if v[2] <= 0:
    #     print(u, v)
    # assert v[2] > 0, \
    #     'Отрицательное давление'
    # assert v[0] > 0, \
    #     'Отрицательная плотность'
    v[0] = v[0] if v[0] >= 0 else 0
    v[3] = v[3] if v[3] >= 0 else 0
    return v


def calc_flux(u):
    f = np.zeros_like(u)
    v = convert_conserved_to_primitive(u)
    f[0] = u[1]  # ρu
    f[1] = u[1] ** 2 / u[0] + v[2]  # ρ * u^2 + p
    f[2] = u[1] * u[2] / u[0]
    f[3] = (u[2] + v[2]) * v[1]  # (E + p) * u

    return f


def calc_speed_of_sound(u, gamma=1.4):
    v = convert_conserved_to_primitive(u)
    # if np.sqrt(gamma * v[2] / v[0]) != np.sqrt(gamma * v[2] / v[0]):
    #     a = 1
    return np.sqrt(gamma * v[2] / v[0])
    # return np.sqrt(np.abs(gamma * v[2] / v[0]))


def _poly_newton_coefficient_batch(x, y):
    """
    x: list or np array contanining x data points
    y: list or np array contanining y data points
    """

    m = len(x)

    x = np.copy(x)
    a = np.copy(y)
    for k in range(1, m):
        a[k:m] = (a[k:m] - a[k - 1]) / (x[k:m] - x[k - 1])

    return a


def newton_polynomial_batch(x_data, y_data, x):
    """
    x_data: data points at x
    y_data: data points at y
    x: evaluation point(s)
    """
    a = np.zeros(np.shape(y_data))
    n = len(x_data) - 1  # Degree of polynomial
    p = np.zeros(np.shape(y_data)[1])
    for i in range(np.shape(y_data)[1]):
        a[i] = _poly_newton_coefficient(x_data, y_data[i])
        p[i] = a[i][n]

        for k in range(1, n + 1):
            p[i] = a[i][n - k] + (x - x_data[n - k]) * p[i]

    return p


def _poly_newton_coefficient(x, y):
    """
    x: list or np array contanining x data points
    y: list or np array contanining y data points
    """

    m = len(x)

    x = np.copy(x)
    a = np.copy(y)
    for k in range(1, m):
        a[k:m] = (a[k:m] - a[k - 1]) / (x[k:m] - x[k - 1])

    return a


def newton_polynomial(x_data, y_data, x):
    """
    x_data: data points at x
    y_data: data points at y
    x: evaluation point(s)
    """
    a = _poly_newton_coefficient(x_data, y_data)
    n = len(x_data) - 1  # Degree of polynomial
    p = a[n]

    for k in range(1, n + 1):
        p = a[n - k] + (x - x_data[n - k]) * p

    return p


def calc_minmod(u_i_minus_1, u_i, u_i_plus_1):
    res = np.zeros_like(u_i)

    for i in range(len(u_i)):
        a = u_i_plus_1[i] - u_i[i]
        b = u_i[i] - u_i_minus_1[i]
        mask = np.abs(a) < np.abs(b)
        res[i] = a * int(mask) + b * int(not mask) if a * b > 0. else 0.

    return res


def check_pressure_and_density(u, gamma=1.4):
    v = convert_conserved_to_primitive(u, gamma=gamma)
    if v[0] < 0 or v[2] < 0:
        res = False
    else:
        res = True

    return res


def calc_divided_difference(u, m_1, m_2, position):
    if m_2 - m_1 == 1:
        next = u[int((m_2 + m_1) / 2)]
    else:
        next = (calc_divided_difference(u, m_1 + 1, m_2, position) - calc_divided_difference(u, m_1,
                                                                                             m_2 - 1,
                                                                                             position)) / (
                       position[int(m_2 - 0.5)] - position[int(m_1 - 0.5)])
    return next


def build_polynomial_for_reconstruction(x_data, y_data, x):
    k = len(x_data)
    p = 0
    h = x_data[-1] - x_data[0] / (k - 1)
    sum = 0
    for j in range(1, k + 1):
        init_sum = 0
        for m in range(0, j):
            P = 1
            for l in range(0, j):
                if l != m:
                    P *= (x - x_data[l] - h / 2)
            init_sum += P
        sum += init_sum * calc_divided_difference(y_data, 0, j, x_data[0:j + 1])

    return sum
