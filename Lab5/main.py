from random import *
from pprint import pprint
from math import sqrt
from scipy.stats import f
from scipy.stats import t as t_check
import sklearn.linear_model as lm
N, d, l = 15, 8, 1.215
print("y` = b0 + b1*x1 + b2*x2 + b3*x3 + b12*x1*x2 + b13x1*x3 + b23*x2*x3 + b123*x1*x2*x3 + b11*x1**2 + b22*x2**2 + b33*x3**2")
def create_mat(N, m):
    mat_sX = [[-4, 3], [-6, 10], [0, 3]]
    mat_1X = [[-1, -1, -1], [-1, -1, 1], [-1, 1, -1], [-1, 1, 1], [1, -1, -1], [1, -1, 1], [1, 1, -1], [1, 1, 1], [-1.215, 0, 0], [1.215, 0, 0], [0, -1.215, 0], [0, 1.215, 0], [0, 0, -1.215], [0, 0, 1.215], [0, 0, 0]]
    mat_2X = [[-1, -1, -1], [-1, -1, 1], [-1, 1, -1], [-1, 1, 1], [1, -1, -1], [1, -1, 1], [1, 1, -1], [1, 1, 1],
              [-1.215, 0, 0], [1.215, 0, 0], [0, -1.215, 0], [0, 1.215, 0], [0, 0, -1.215], [0, 0, 1.215], [0, 0, 0]]
    x_min_max = [sum(mat_sX[i][k] for i in range(3))/3 for k in range(2)]
    y_min_max = [int(200 + x_min_max[i]) for i in range(2)]
    print('Задана матриця Х:\n', mat_sX)
    print('Xср min and max:\n', x_min_max, '\nY min and max:\n', y_min_max, '\nМатриця Y:')
    mat_Y = [[randint(y_min_max[0], y_min_max[1]) for i in range(m)] for k in range(N)]
    pprint(mat_Y)
    mat_serY = [sum(mat_Y[k1])/m for k1 in range(N)]
    print('Середні значення Y:\n', mat_serY, '\nНормована матриця Х:')
    mat_X0i = [sum(mat_sX[i]) / 2 for i in range(3)]
    mat_dX = [mat_sX[i][1] - mat_X0i[i] for i in range(3)]
    mat_X = [[], [], [], [], [], [], [], [], [], [], [], [], [], [], []]
    for i in range(15):
        for j in range(3, 10):
                mat_2X[i].append(mat_2X[i][0] * mat_2X[i][1])
                mat_2X[i].append(mat_2X[i][0] * mat_2X[i][2])
                mat_2X[i].append(mat_2X[i][1] * mat_2X[i][2])
                mat_2X[i].append(mat_2X[i][0] * mat_2X[i][1] * mat_2X[i][2])
                mat_2X[i].append(mat_2X[i][0] ** 2)
                mat_2X[i].append(mat_2X[i][1] ** 2)
                mat_2X[i].append(mat_2X[i][2] ** 2)
                break
    for i in range(15):
        mat_2X[i].insert(0, 1)
    pprint(mat_2X)
    for i in range(8):
        for j in range(10):
            if j < 3:
                if mat_1X[i][j] == -1:
                    mat_X[i].append(mat_sX[j][0])
                else:
                    mat_X[i].append(mat_sX[j][1])
            if j > 3:
                mat_X[i].append(mat_X[i][0] * mat_X[i][1])
                mat_X[i].append(mat_X[i][0] * mat_X[i][2])
                mat_X[i].append(mat_X[i][1] * mat_X[i][2])
                mat_X[i].append(mat_X[i][0] * mat_X[i][1] * mat_X[i][2])
                mat_X[i].append(mat_X[i][0] ** 2)
                mat_X[i].append(mat_X[i][1] ** 2)
                mat_X[i].append(mat_X[i][2] ** 2)
                break

    for i in range(8, 15):
        for j in range(10):
            if j < 3:
                if mat_1X[i][j] == 0:
                    mat_X[i].append(mat_X0i[j])
                else:
                    mat_X[i].append(mat_1X[i][j] * mat_dX[j] + mat_X0i[j])
            else:
                mat_X[i].append(mat_X[i][0] * mat_X[i][1])
                mat_X[i].append(mat_X[i][0] * mat_X[i][2])
                mat_X[i].append(mat_X[i][1] * mat_X[i][2])
                mat_X[i].append(mat_X[i][0] * mat_X[i][1] * mat_X[i][2])
                mat_X[i].append(mat_X[i][0] ** 2)
                mat_X[i].append(mat_X[i][1] ** 2)
                mat_X[i].append(mat_X[i][2] ** 2)
                break

    return (mat_X, mat_Y, mat_serY, mat_2X)

def coef_b(x, y):
    for i in range(15):
        x[i].insert(0, 1)
    print('Нормалізована матриця Х:')
    pprint(x, width=150)
    skm = lm.LinearRegression(fit_intercept=False)
    skm.fit(x, y)
    b = skm.coef_
    b = [round(i, 3) for i in b]
    print(b)

    return b


def check(mat_serY, mat_Y, tran1, blist, matt_fullX, m):
    global b1
    d = 11
    b = 0
    b1 = 0
    d1 = 11
    mat_disY = [sum([((k1 - mat_serY[j]) ** 2) for k1 in mat_Y[j]]) / m for j in range(N)]
    print("Дисперсії в рядках:\n", mat_disY)
    print(
        '-------------------------------------------------------------------------\nПЕРЕВІРКА ОДНОРІДНОСТІ ДИСПЕРСІЇ ЗА КРИТЕРІЄМ КОХРЕНА:\n\n...\n..')
    if max(mat_disY) / sum(mat_disY) < 0.7679:
        print('Дисперсія однорідна')
    else:
        print('Дисперсія неоднорідна')
        m += 1
        main(N, m)
    print(
        '-------------------------------------------------------------------------\nПЕРЕВІРКА ЗНАЧУЩОСТІ КОЕФІЦІЄНТІВ ЗА КРИТЕРІЄМ СТЬЮДЕНТА:\n')
    S2b = sum(mat_disY) / N
    S2bs = S2b / (m * N)
    Sbs = sqrt(S2bs)
    print('Sbs:\n', Sbs)
    bb = [sum(mat_serY[k] * tran1[i][k] for k in range(N)) / N for i in range(d)]
    t = [abs(bb[i]) / Sbs for i in range(d1)]
    print('bi:\n', bb, '\nti:\n', t, '\n...\n..')
    f1, f2 = m - 1, N
    f3 = f1 * f2
    for i in range(d1):
        if t[i] < t_check.ppf(q=0.975, df=f3):
            blist[i] = 0
            d -= 1
            b += 1
            print('Виключаємо з рівняння коефіціент b', i)
    y_reg = [
        blist[0] * matt_fullX[i][0] + blist[1] * matt_fullX[i][1] + blist[2] * matt_fullX[i][2] + blist[3] * matt_fullX[i][3] + blist[4] *
        matt_fullX[i][4] + blist[5] * matt_fullX[i][5] + blist[6] * matt_fullX[i][6] + blist[7] * matt_fullX[i][7] + blist[8] * matt_fullX[i][8] + blist[9] * matt_fullX[i][9] + blist[10] * matt_fullX[i][10] for i
        in range(N)]
    print('Значення рівнянь регресій:\n', y_reg)
    print(
        '-------------------------------------------------------------------------\nПЕРЕВІРКА АДЕКВАТНОСТІ ЗА КРИТЕРІЄМ ФІШЕРА:\n')
    f4 = N - d
    Sad = (m / (N - d)) * int(sum(y_reg[i] - mat_serY[i] for i in range(N)) ** 2)
    Fp = Sad / S2b
    b1 += b
    print('Кількість значимих коефіціентів:\n', d, '\nFp:\n', Fp, '\n...\n..')
    if Fp > f.ppf(q=0.95, dfn=f4, dfd=f3):
        print('Рівняння регресії неадекватно оригіналу при рівні значимості 0.05')
    else:
        print('Рівняння регресії адекватно оригіналу при рівні значимості 0.05')


def main(N, m):
    b2 = 0
    for j in range(100):
        x, mat_Y, mat_serY, mat_2X = create_mat(N, m)
        tran1 = [list(i) for i in zip(*mat_2X)]
        b = coef_b(x, mat_serY)
        check(mat_serY, mat_Y, tran1, b, x, m)
        b2 += b1
    print("Cередня кількість незначимих коефіцієнтів на 1 ітерацію: ", b2 / 100)


main(15, 3)
