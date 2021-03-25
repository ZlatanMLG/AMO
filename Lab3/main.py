from random import *
from pprint import pprint
import numpy as np
from math import sqrt
m, N, d = 3, 4, 4
mat_sX = [[-25, -5], [-30, 45], [-5, 5]]
mat_1X = [[1, -1, -1, -1], [1, -1, 1, 1], [1, 1, -1, 1], [1, 1, 1, -1]]
tran1 = [list(i) for i in zip(*mat_1X)]
x_min_max = [sum(mat_sX[i][k] for i in range(3))/3 for k in range(2)]
y_min_max = [int(200 + x_min_max[i]) for i in range(2)]
print('Задана матриця Х:\n', mat_sX, '\nНормована матриця Х:')
pprint(mat_1X, width=17)
print('Xср min and max:\n', x_min_max, '\nY min and max:\n', y_min_max, '\nМатриця Y:')
mat_Y = [[randint(y_min_max[0], y_min_max[1]) for i in range(3)] for k in range(4)]
pprint(mat_Y, width=17)
mat_serY = [sum(mat_Y[k1])/3 for k1 in range(4)]
print('Середні значення Y:\n', mat_serY, '\nМатриця Х:')
mat_X = [[-25, -30, -5], [-25, 45, 5], [-5, -30, 5], [-5, 45, -5]]
pprint(mat_X, width=17)
mx = [sum(mat_X[i][k] for i in range(4))/4 for k in range(3)]
my = sum(mat_serY)/4
print('Значення mxi:\n', mx, '\nЗначення my:\n', my)
tran = [list(i) for i in zip(*mat_X)]
ai = [sum(tran[k][i] * mat_serY[i] for i in range(4)) / 4 for k in range(3)]
aii = [sum(tran[k][i]**2 for i in range(4)) / 4 for k in range(3)]
print('Значення ai:\n', ai, '\nЗначення aii:\n', aii)
a12 = a21 = (tran[0][0] * tran[1][0] + tran[0][1] * tran[1][1] + tran[0][2] * tran[1][2] + tran[0][3] * tran[1][3])/4
a13 = a31 = (tran[0][0] * tran[2][0] + tran[0][1] * tran[2][1] + tran[0][2] * tran[2][2] + tran[0][3] * tran[2][3])/4
a23 = a32 = (tran[1][0] * tran[2][0] + tran[1][1] * tran[2][1] + tran[1][2] * tran[2][2] + tran[1][3] * tran[2][3])/4
znamen = np.linalg.det(np.matrix([[1, mx[0], mx[1], mx[2]], [mx[0], aii[0], a12, a13], [mx[1], a12, aii[1], a32], [mx[2], a13, a23, aii[2]]]))
b0 = np.linalg.det(np.matrix([[my, mx[0], mx[1], mx[2]], [ai[0], aii[0], a12, a13], [ai[1], a12, aii[1], a32], [ai[2], a13, a23, aii[2]]]))/znamen
b1 = np.linalg.det(np.matrix([[1, my, mx[1], mx[2]], [mx[0], ai[0], a12, a13], [mx[1], ai[1], aii[1], a32], [mx[2], ai[2], a23, aii[2]]]))/znamen
b2 = np.linalg.det(np.matrix([[1, mx[0], my, mx[2]], [mx[0], aii[0], ai[0], a13], [mx[1], a12, ai[1], a32], [mx[2], a13, ai[2], aii[2]]]))/znamen
b3 = np.linalg.det(np.matrix([[1, mx[0], mx[1], my], [mx[0], aii[0], a12, ai[0]], [mx[1], a12, aii[1], ai[1]], [mx[2], a13, a23, ai[2]]]))/znamen
perevirka = [b0 + b1 * tran[0][i] + b2 * tran[1][i] + b3 * tran[2][i] for i in range(4)]
blist = [b0, b1, b2, b3]
print("Перевірка порівнянням з середніми значеннями Y:\n", perevirka)
mat_disY = [sum([((k1 - mat_serY[j]) ** 2) for k1 in mat_Y[j]]) / m for j in range(4)]
print("Дисперсії в рядках:\n", mat_disY)
print('-------------------------------------------------------------------------\nПЕРЕВІРКА ОДНОРІДНОСТІ ДИСПЕРСІЇ ЗА КРИТЕРІЄМ КОХРЕНА:\n\n...\n..')
if max(mat_disY)/sum(mat_disY) < 0.7679:
    print('Дисперсія однорідна')
else:
    print('Дисперсія неоднорідна')
print('-------------------------------------------------------------------------\nПЕРЕВІРКА ЗНАЧУЩОСТІ КОЕФІЦІЄНТІВ ЗА КРИТЕРІЄМ СТЬЮДЕНТА:\n')
S2b = sum(mat_disY) / N
S2bs = S2b / (m * N)
Sbs = sqrt(S2bs)
print('Sbs:\n', Sbs)
bb = [sum(mat_serY[k] * tran1[i][k] for k in range(N))/N for i in range(N)]
t = [abs(bb[i])/Sbs for i in range(N)]
print('bi:\n', bb, '\nti:\n', t, '\n...\n..')
for i in range(N):
    if t[i] < 2.306:
        blist[i] = 0
        d -= 1
        print('Виключаємо з рівняння коефіціент b', i)
y_reg = [blist[0] + blist[1] * mat_X[i][0] + blist[2] * mat_X[i][1] + blist[3] * mat_X[i][2] for i in range(4)]
print('Значення рівнянь регресій:\n', y_reg)
print('-------------------------------------------------------------------------\nПЕРЕВІРКА АДЕКВАТНОСТІ ЗА КРИТЕРІЄМ ФІШЕРА:\n')
Sad = (m / (N - d)) * int(sum(y_reg[i] - mat_serY[i] for i in range(N))**2)
Fp = Sad / S2b
print('Кількість значимих коефіціентів:\n', d, '\nFp:\n', Fp, '\n...\n..')
if Fp > 4.5:
    print('Рівняння регресії неадекватно оригіналу при рівні значимості 0.05')
else:
    print('Рівняння регресії адекватно оригіналу при рівні значимості 0.05')