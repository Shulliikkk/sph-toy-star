import numpy as np
import matplotlib.pyplot as plt
from scipy.special import gamma

def W(abs_r, h):
    '''
    Функция сглаживающего ядра
    В нашем случае гауссиан (нормальное распределение Гаусса)
    Задает распределение массы вокруг частицы
    r - расстояние от частицы
    h - длина сглаживания, та длина, на которой частица может действовать на другие частицы
    (Возвращает значаение сглаживающего ядра на расстоянии r)
    '''
    w = 1 / (h**3 + np.pi**1.5) * np.exp(-abs_r**2 / h**2)
    return w

def gradW(vec_r, h):
    '''
    Градиент сглаживающего ядра W
    vec_r - массив координат направляющего вектора r,
    направленного от частицы в какую-либо точку пространства,
    в нашем случае от i-й к j-й частице
    (Возвращает массив частных производных (wx, wy))
    '''
    # вычисление нормы вектора r
    r = np.linalg.norm(vec_r)
    # вычисление скалярной части градиента
    scalar_deriv = -2 / (h**5 * np.pi**1.5) * np.exp(-r**2 / h**2)
    return scalar_deriv * vec_r

def density(pos, m, h):
    '''
    Вычисление плотности в точке нахождения каждой из частиц
    pos = np.array([x1, y1], ..., [xn, yn]) - массив координат всех частиц
    m - масса частицы
    (Возвращает массив плотностей)
    '''
    N = pos.shape[0]
    rho = np.zeros(N)
    for i in range(N):
        # инициализация плотности i-й частицы
        rho[i] = m * W(0, h)
        for j in range(i + 1, N):
            # вычисление вектора между i-й и j-й частицами
            R_ij = pos[i, :] - pos[j, :]
            # вклад в плотность i-й и j-й частицы
            rho_ij = m * W(np.linalg.norm(R_ij), h)
            # добавление вклада
            rho[i] += rho_ij
            rho[j] += rho_ij
    return rho

def pressure(density, k, n):
    '''
    Вычисление давления в точке нахождения каждой из частиц
    pos = np.array([x1, y1], ..., [xn, yn]) - массив координат всех частиц
    k - поправочный коэффициент
    n - показатель политропы
    (Возвращает массив давлений)
    '''
    return k * density**(1 + 1 / n)

def acceleration(pos, vel, m, h, k, n, lmbda, nu):
    '''
    Вычисление ускорения каждой из частиц
    lmbda, nu - гравитационный и вязкостной коэффициенты
    (Возвращает массив ускорений
    Формат записи массива аналогичен массиву координат)
    '''
    N = pos.shape[0]
    a = np.zeros((N, 2))
    rho = density(pos, m, h)
    p = pressure(rho, k, n)
    for i in range(N):
        # добавление вкалада в ускорение из-за гравитации и вязкости
        a[i, :] += -nu * vel[i, :] - lmbda * pos[i, :]
        for j in range(i + 1, N):
            R_ij = pos[i, :] - pos[j, :]
            # добавление вкалада в ускорение из-за давления
            a_p = - m * (p[i] / rho[i]**2 + p[j] / rho[j]**2) * gradW(R_ij, h)
            a[i, :] += a_p
            a[j, :] += -a_p
    return a

def main():
    # константы
    tEnd = 12
    dt = 0.04
    Nt = int(np.ceil(tEnd/dt))
    N = 200
    M = 2
    R = 0.75
    m = M / N
    h = 0.1
    k = 0.1
    n = 1
    nu = 1
    lmbda = 4 # 2*k*(1+n)*np.pi**(-3/(2*n)) * (M*gamma(5/2+n)/R**3/gamma(1+n))**(1/n) / R**2 ~ 2.1

    # начальные параметры
    pos = np.random.standard_normal(size = (N, 2))
    vel = np.random.standard_normal(size = (N, 2))
    acc = acceleration(pos, vel, m, h, k, n, lmbda, nu)

    # формат фигуры для отрисовки
    fig = plt.figure(figsize=(4,5), dpi=80)
    grid = plt.GridSpec(3, 1, wspace=0.0, hspace=0.3)
    ax1 = plt.subplot(grid[0:2,0])
    plotRealTime = True

    # численное интегрирование методом лягушки
    for i in range(Nt):
        vel += acc * dt / 2
        pos += vel * dt
        acc = acceleration(pos, vel, m, h, k, n, lmbda, nu)
        vel += acc * dt / 2
        rho = density(pos, m, h)
        # отрисовка
        if plotRealTime or (i == Nt-1):
            plt.sca(ax1)
            plt.cla()
            cval = np.minimum((rho - 3) / 3, 1).flatten()
            plt.scatter(pos[:, 0], pos[:, 1], c = cval, cmap = plt.cm.autumn, s = 10, alpha = 0.5)
            ax1.set(xlim = (-1.4, 1.4), ylim = (-1.2, 1.2))
            ax1.set_aspect('equal', 'box')
            ax1.set_xticks([-1, 0, 1])
            ax1.set_yticks([-1, 0, 1])
            ax1.set_facecolor('black')
            ax1.set_facecolor((.1, .1, .1))
            plt.pause(0.0001)

    plt.savefig('sph.png',dpi=240)
    plt.show()
    return 0

if __name__== "__main__":
    main()
