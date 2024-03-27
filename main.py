import numpy as np
import matplotlib.pyplot as plt
import math

x_1 = [-5.0,
       -4.473684210526316,
       -3.947368421052632,
       -3.4210526315789473,
       -2.8947368421052633,
       -2.368421052631579,
       -1.8421052631578947,
       -1.3157894736842106,
       -0.7894736842105265,
       -0.2631578947368425,
       0.2631578947368416,
       0.7894736842105257,
       1.3157894736842106,
       1.8421052631578947,
       2.3684210526315788,
       2.894736842105263,
       3.421052631578947,
       3.947368421052632,
       4.473684210526315,
       5.0]
y_1 = [0.04244131815783876,
       -0.03743241795177017,
       0.013197055039641211,
       0.024917014083434722,
       -0.06447222549472982,
       0.08837633378106816,
       -0.0780213900874819,
       0.013318149842937815,
       0.14701693836654636,
       -0.7626928790095355,
       -0.7626928790095371,
       0.14701693836654556,
       0.013318149842937815,
       -0.0780213900874819,
       0.08837633378106824,
       -0.06447222549472971,
       0.024917014083434524,
       0.013197055039641211,
       -0.037432417951770394,
       0.04244131815783876]


def chebyshev(st, fin, n: int):
    roots = []
    for k in range(n):
        root = 0.5 * (st + fin) + 0.5 * (fin - st) * math.cos((2 * k + 1) * math.pi / (2 * n))
        roots.append(root)
    return roots


def draw_graph(x, y, num_dots_between):
    """array of x, array of y,
    num of new dots between first ones"""

    def divided_difference(i: int, k: int):
        sum = 0
        for m in range(i, i + k + 1):
            mul = 1
            for j in range(i, i + k + 1):
                if j != m:
                    mul *= x[m] - x[j]
            sum += 0 if mul == 0 else y[m] * (mul ** -1)

        return sum

    def newton_interpolation(x_this: float, n: int):
        sum = 0
        for k in range(n):
            mul = 1
            for m in range(k):
                mul *= x_this - x[m]
            div_dif = divided_difference(0, k)
            sum += div_dif * mul

        return sum

    def add_points(k: int):
        """add k point with same step between old ones"""
        new_x = []
        for i in range(len(x) - 1):
            new_x.extend(np.linspace(x[i], x[i + 1], num=k + 2)[:-1])
        new_x.append(x[-1])
        return new_x

    plt.figure(figsize=(10, 5))

    l = len(x)
    plt.scatter(x, y, marker="o", s=100, color='r', facecolors='none', label="Узлы функции")
    plt.plot(x, y, linestyle='--', color='r')

    x_new = add_points(num_dots_between)
    y_new = []

    for xi in x_new:
        y_new.append(newton_interpolation(xi, l))

    plt.scatter(x_new, y_new, marker="*", s=30,
                color='b', facecolors='none', label="Узлы интерполяции Ньютона")
    plt.plot(x_new, y_new, linestyle='--', color='black', label="-sinc(-1,5x)")

    plt.xlabel('x', fontsize=14)
    plt.ylabel('y', fontsize=14)
    plt.tick_params(labelsize=14)
    plt.legend(fontsize=8)

    plt.show()


# draw_graph(x_1, y_1, 1) # Данные из таблицы

st = -5
fin = 5
ln = 20
x_ch = chebyshev(st, fin, ln)
y_ch = []
for xi in x_ch:
    y_ch.append(- np.sinc(-1.5 * xi))

draw_graph(x_ch, y_ch, 1)  # Данные по точкам из многочлена Чебышева
