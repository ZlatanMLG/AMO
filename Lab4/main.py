from random import *
from pprint import pprint
from math import sqrt
from scipy.stats import f
from scipy.stats import t as t_check
N = 8
print("y` = b0 + b1*x1 + b2*x2 + b3*x3 + b12*x1*x2 + b13x1*x3 + b23*x2*x3 + b123*x1*x2*x3")


def create_mat(N, m):
    mat_sX = [[-20, 15], [10, 60], [15, 35]]
    mat_1X = [[1, -1, -1, -1, 1, 1, 1, -1], [1, -1, -1, 1, 1, -1, -1, 1], [1, -1, 1, -1, -1, 1, -1, 1], [1, -1, 1, 1, -1, -1, 1, -1], [1, 1, -1, -1, -1, -1, 1, 1], [1, 1, -1, 1, -1, 1, -1, -1], [1, 1, 1, -1, 1, -1, -1, -1], [1, 1, 1, 1, 1, 1, 1, 1]]
    tran1 = [list(i) for i in zip(*mat_1X)]
    x_min_max = [sum(mat_sX[i][k] for i in range(3))/3 for k in range(2)]
    y_min_max = [int(200 + x_min_max[i]) for i in range(2)]
    print('Задана матриця Х:\n', mat_sX, '\nНормована матриця Х:')
    pprint(mat_1X)
    print('Xср min and max:\n', x_min_max, '\nY min and max:\n', y_min_max, '\nМатриця Y:')
    mat_Y = [[randint(y_min_max[0], y_min_max[1]) for i in range(m)] for k in range(8)]
    pprint(mat_Y)
    mat_serY = [sum(mat_Y[k1])/m for k1 in range(8)]
    print('Середні значення Y:\n', mat_serY, '\nМатриця Х:')
    mat_X = [[-20, 10, 15], [-20, 10, 35], [-20, 60, 15], [-20, 60, 35], [15, 10, 15], [15, 10, 35], [15, 60, 15], [15, 60, 35]]
    pprint(mat_X)
    mx = [sum(mat_X[i][k] for i in range(8))/8 for k in range(3)]
    my = sum(mat_serY)/8
    print('Значення mxi:\n', mx, '\nЗначення my:\n', my)
    tran = [list(i) for i in zip(*mat_X)]
    b0 = sum(mat_serY[i] for i in range(N)) / N
    b1 = sum(mat_serY[i] * tran1[1][i] for i in range(N)) / N
    b2 = sum(mat_serY[i] * tran1[2][i] for i in range(N)) / N
    b3 = sum(mat_serY[i] * tran1[3][i] for i in range(N)) / N
    b12 = sum(mat_serY[i] * tran1[1][i] * tran1[2][i] for i in range(N)) / N
    b13 = sum(mat_serY[i] * tran1[1][i] * tran1[3][i] for i in range(N)) / N
    b23 = sum(mat_serY[i] * tran1[2][i] * tran1[3][i] for i in range(N)) / N
    b123 = sum(mat_serY[i] * tran1[1][i] * tran1[2][i] * tran1[3][i] for i in range(N)) / N
    blist = [b0, b1, b2, b3, b12, b13, b23, b123]
    matt_fullX = [[-20, 10, 15, -200, -300, 150, -3000], [-20, 10, 35, -200, -700, 350, -7000], [-20, 60, 15, -1200, -300, 900, -18000], [-20, 60, 35, -1200, -700, 2100, -42000], [15, 10, 15, 150, 225, 150, 2250], [15, 10, 35, 150, 525, 350, 5250], [15, 60, 15, 900, 225, 900, 13500], [15, 60, 35, 900, 525, 2100, 31500]]
    mat_fullX = [list(i) for i in zip(*matt_fullX)]
    pprint(mat_fullX)
    y_result = [b0 + b1 * mat_fullX[0][i] + b2 * mat_fullX[1][i] + b3 * mat_fullX[2][i] + b12 * mat_fullX[0][i] * mat_fullX[1][i] + b13 * mat_fullX[0][i] * mat_fullX[2][i] + b23 * mat_fullX[1][i] * mat_fullX[2][i] + b123 * mat_fullX[0][i] * mat_fullX[1][i] * mat_fullX[2][i] for i in range(N)]
    print(y_result)
    return (mat_serY, mat_Y, tran1, blist, matt_fullX)

def check(mat_serY, mat_Y, tran1, blist, matt_fullX, m):
    d = 8
    mat_disY = [sum([((k1 - mat_serY[j]) ** 2) for k1 in mat_Y[j]]) / m for j in range(N)]
    print("Дисперсії в рядках:\n", mat_disY)
    print('-------------------------------------------------------------------------\nПЕРЕВІРКА ОДНОРІДНОСТІ ДИСПЕРСІЇ ЗА КРИТЕРІЄМ КОХРЕНА:\n\n...\n..')
    if max(mat_disY)/sum(mat_disY) < 0.7679:
        print('Дисперсія однорідна')
    else:
        print('Дисперсія неоднорідна')
        m += 1
        main(N, m)
    print('-------------------------------------------------------------------------\nПЕРЕВІРКА ЗНАЧУЩОСТІ КОЕФІЦІЄНТІВ ЗА КРИТЕРІЄМ СТЬЮДЕНТА:\n')
    S2b = sum(mat_disY) / N
    S2bs = S2b / (m * N)
    Sbs = sqrt(S2bs)
    print('Sbs:\n', Sbs)
    bb = [sum(mat_serY[k] * tran1[i][k] for k in range(N))/N for i in range(N)]
    t = [abs(bb[i])/Sbs for i in range(N)]
    print('bi:\n', bb, '\nti:\n', t, '\n...\n..')
    f1, f2 = m - 1, N
    f3 = f1 * f2
    for i in range(N):
        if t[i] < t_check.ppf(q=0.975, df=f3):
            blist[i] = 0
            d -= 1
            print('Виключаємо з рівняння коефіціент b', i)
    y_reg = [blist[0] + blist[1] * matt_fullX[i][0] + blist[2] * matt_fullX[i][1] + blist[3] * matt_fullX[i][2] + blist[4] * matt_fullX[i][3] + blist[5] * matt_fullX[i][4] + blist[6] * matt_fullX[i][5] + blist[7] * matt_fullX[i][6] for i in range(N)]
    print('Значення рівнянь регресій:\n', y_reg)
    print('-------------------------------------------------------------------------\nПЕРЕВІРКА АДЕКВАТНОСТІ ЗА КРИТЕРІЄМ ФІШЕРА:\n')
    f4 = N - d
    Sad = (m / (N - d)) * int(sum(y_reg[i] - mat_serY[i] for i in range(N))**2)
    Fp = Sad / S2b
    print('Кількість значимих коефіціентів:\n', d, '\nFp:\n', Fp, '\n...\n..')
    if Fp > f.ppf(q=0.95, dfn=f4, dfd=f3):
        print('Рівняння регресії неадекватно оригіналу при рівні значимості 0.05')
    else:
        print('Рівняння регресії адекватно оригіналу при рівні значимості 0.05')


def main(N, m):
    mat_serY, mat_Y, tran1, blist, matt_fullX = create_mat(N, m)
    check(mat_serY, mat_Y, tran1, blist, matt_fullX, m)


main(8, 3)
