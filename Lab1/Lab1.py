import time
start = time.time()
from random import *
from pprint import pprint
mat = [[randint(1, 20) for j in range(3)] for i in range(8)]
a_mat = [randint(1, 20) for k in range(4)]
print("--------------------------------------------------------\nЗгенерована матриця Х: ")
pprint(mat)
print("--------------------------------------------------------\nЗгенеровані значення а:\n", a_mat)
Y = [(a_mat[0]+a_mat[1]*mat[i][0]+a_mat[2]*mat[i][1]+a_mat[3]*mat[i][2]) for i in range(8)]
print("--------------------------------------------------------\nСписок значень Y:\n", Y)
x1_mat = [mat[i][0] for i in range(8)]
x01 = (max(x1_mat) + min(x1_mat)) / 2
x2_mat = [mat[i][1] for i in range(8)]
x02 = (max(x2_mat) + min(x2_mat)) / 2
x3_mat = [mat[i][2] for i in range(8)]
x03 = (max(x3_mat) + min(x3_mat)) / 2
x0_mat = [x01, x02, x03]
print("--------------------------------------------------------\nЗначення Х0:\n", x0_mat)
dx1 = x01 - min(x1_mat)
dx2 = x02 - min(x2_mat)
dx3 = x03 - min(x3_mat)
dx_mat = [dx1, dx2, dx3]
print("--------------------------------------------------------\nЗначення dX:\n", dx_mat)
Norm_mat = [[round((mat[i][j]-x0_mat[j])/dx_mat[j], 5) for j in range(3)] for i in range(8)]
print("--------------------------------------------------------\nНормовані значення:")
pprint(Norm_mat)
print("--------------------------------------------------------\nЗначення Y ет:\n", round((a_mat[0]+a_mat[1]*x01+a_mat[2]*x02+a_mat[3]*x03), 3))
print("--------------------------------------------------------\nЗначення шуканого min(Y):\n", min(Y))
print("--------------------------------------------------------\nЗначення Х для min(Y):\n", mat[Y.index(min(Y))])

print("Час роботи програми в секундах:", time.time() - start)
