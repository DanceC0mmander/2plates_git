# -*- coding: utf-8 -*-

# когда надо сбросить время, надо один блок закомментировать и одну строку раскомментировать (они отмечены комментариями)
# при дальнейших расчетах вернуть все

import numpy as np
import math
import matplotlib.pyplot as plt
import argparse
import os
import ffmpeg
import shutil
# import bigfloat
#from mpmath import *

# константы задачи
alfa = 1
gamma = 0.85
ksi = 0.1
q = 0.9
k = 0.2
beta1 = 6
beta2 = 6
sll = 1
c = 0.3#0.372748

#вспомогательные константы
eps = 1e-4 # для сравнения с нулем
eps_exp = 0.01

# шаг по пространству
dx = 0.1#0.05
x_max = 500
N = int(x_max/dx)
x = np.arange(0, x_max, dx)

print('dx=', dx)
print('N=', N)
print('len(x)', len(x))

# инициация переменных

T2_1 = []
Y2_1 = []
T1_1 = []
Y1_1 = []
T1 = [None]*N
Y1 = [None]*N
T2 = [None]*N
Y2 = [None]*N
chem_rate = [None]*N

# файлы
filename_T1 = 'T1.dat'
filename_Y1 = 'Y1.dat'
filename_T2 = 'T2.dat'
filename_Y2 = 'Y2.dat'
filename_overall = 'overall.dat'
filename_T_stat = 'T_stat.dat'
filename_parameters = 'param.dat'
filename_time = 'time.dat'
filename_x = 'x.dat'
filename_chem_rate = 'chem_rate.dat'

# считать значение текущего времени из файла. Если такого файла нет или он пустой, то текущее время t_cur = 0
subdir = ''#'1sec/'
if (os.path.isfile(subdir + filename_time) and os.stat(subdir + filename_time).st_size != 0):
	f = open(subdir + filename_time, 'r')
	t_cur = float(f.readline())
	print('Current time: %f' %t_cur)
else:
	t_cur = 0.
# t_max должно быть больше, чем t_cur
# t_cur = 0
t_max = float(input('Enter time: '))
# закомментировать, когда надо сбросить время
# while(t_max <= t_cur):
# 	print('Err: Time must be greater than current time')
# 	t_max = float(input('Enter time: '))

# шаг по времени
dt = 4e-03#5e-04#
print("dt = %f" %dt)

def exponenta(beta_var, T_var):
	return math.exp( beta_var*(T_var-1)/(1+gamma*(T_var-1)) )

# адиабатическая температура
Ub = (q + sll)/(1 + sll)


# начальные значения
if t_cur > eps:
	# считать с отдельный файлов
	with open(subdir + 'T1.dat') as file_T1:
		for line in file_T1:
			T1_1.append(float(line))
	file_T1.close()
	with open(subdir + 'Y1.dat') as file_Y1:
		for line in file_Y1:
			Y1_1.append(float(line))
	file_Y1.close()
	with open(subdir + 'T2.dat') as file_T2:
		for line in file_T2:
			T2_1.append(float(line))
	file_T2.close()
	with open(subdir + 'Y2.dat') as file_Y2:
		for line in file_Y2:
			Y2_1.append(float(line))
	file_Y2.close()
	with open(subdir + 'param.dat') as file_parameters:
		next(file_parameters)
		for line in file_parameters:
			c = float(line.split()[5])	
	file_parameters.close()
else:
	# назначить новое
	for i in range(N):
		T1_1.append((0.5-0.5*np.tanh(x[i]-230)))
		Y1_1.append((0.5+0.5*np.tanh(x[i]-230)))
		T2_1.append((0.5-0.5*np.tanh(x[i]-230)))
		Y2_1.append((0.5+0.5*np.tanh(x[i]-230)))
		# T2_1.append(0)
		# Y2_1.append(1)


# раскомментировать, когда надо сбросить время
t_cur = 0

# для вычисления скорости
sum_T1 = 0
for i in range(1, N-1):
	sum_T1 = sum_T1 + T1_1[i]


# открыть файл для записи значений beta, T1, Y1, T2
file_overall = open(filename_overall, 'w')
file_overall.write('time   T1   Y1   T2   Y2\n')

file_chem_rate = open(filename_chem_rate, 'w')

# константы для сокращений
var1 = dt/(dx*dx)

time_steps = 0

# константа демпфирующего члена
w_z = 1e-03

while t_cur <= t_max:

	c_1 = 0

	while not math.isclose(c, c_1, rel_tol=0.01, abs_tol=0.0):

		c_1 = c
		sum_T1_1 = sum_T1
		sum_T1 = 0

		var2 = 0.5*c*dt/dx

		#----- 1 уравнение -----

		# i = 1..N-2
		for i in range(1, N-1):
			
			Tm1 = max(T1_1[i], eps_exp)
			Tm2 = max(T2_1[i], eps_exp)
			var3_1 = dt*Y1_1[i]*exponenta(beta1, Tm1)
			var3_2 = dt*Y2_1[i]*exponenta(beta2, Tm2)
			var4 = dt*ksi*(T1_1[i]-T2_1[i])

			
			T1[i] = T1_1[i] + var1*(T1_1[i+1]-2*T1_1[i]+T1_1[i-1]) + var3_1 + var2*(T1_1[i+1]-T1_1[i-1]) - var4
			Y1[i] = Y1_1[i] - var3_1 + var2*(Y1_1[i+1]-Y1_1[i-1])
			T2[i] = T2_1[i] + alfa*var1*(T2_1[i+1]-2*T2_1[i]+T2_1[i-1]) + q*k*var3_2 + var2*(T2_1[i+1]-T2_1[i-1]) + sll*var4
			Y2[i] = Y2_1[i] - k*var3_2 + var2*(Y2_1[i+1]-Y2_1[i-1])


			# введение вязкости
			if i != 1 and i != N-2:
				w1 = Y1_1[i+2] - 4*Y1_1[i+1] + 6*Y1_1[i] - 4*Y1_1[i-1] + Y1_1[i-2]
				w2 = Y2_1[i+2] - 4*Y2_1[i+1] + 6*Y2_1[i] - 4*Y2_1[i-1] + Y2_1[i-2]
			elif i == 1:
				w1 = Y1_1[3] - 4*Y1_1[2] - 5*Y1_1[1] - 2*Y1_1[0]
				w2 = Y2_1[3] - 4*Y2_1[2] - 5*Y2_1[1] - 2*Y2_1[0]
			elif i == N-2:
				w1 = -2*Y1_1[N-1] + 5*Y1_1[N-2] - 4*Y1_1[N-3] + Y1_1[N-4] 
				w2 = -2*Y2_1[N-1] + 5*Y2_1[N-2] - 4*Y2_1[N-3] + Y2_1[N-4] 

			Y1[i] -= w_z*w1
			Y2[i] -= w_z*w2

			# для расчета скорости волны
			sum_T1 += T1[i]

			# скорость реакции (для выявления пульсаций)
			Tm1 = max(T1[i], eps_exp)
			chem_rate[i] = Y1[i]*exponenta(beta1, Tm1) 


		# граничные условия
		# if 2*T1[1]-T1[2] <= Ub:
		# 	T1[0] = 2*T1[1]-T1[2]
		# 	Y1[0] = 2*Y1[1]-Y1[2]
		# else:
		# 	T1[0] = Ub
		# 	Y1[0] = 0
		T1[0] = 2*T1[1]-T1[2]
		Y1[0] = 2*Y1[1]-Y1[2]
		T1[N-1] = 0
		Y1[N-1] = 1

		# if 2*T2[1]-T2[2] <= Ub:
		# 	T2[0] = 2*T2[1]-T2[2]
		# 	Y2[0] = 2*Y2[1]-Y2[2]
		# else:
		# 	T2[0] = Ub
		# 	Y2[0] = 0
		T2[0] = 2*T2[1]-T2[2]
		Y2[0] = 2*Y2[1]-Y2[2]
		T2[N-1] = 0
		Y2[N-1] = 1

		# скорость реации на границах
		Tm1 = max(T1[0], eps_exp)
		chem_rate[0] = Y1[0]*exponenta(beta1, Tm1) 
		Tm1 = max(T1[N-1], eps_exp)
		chem_rate[N-1] = Y1[N-1]*exponenta(beta1, Tm1) 

		# корректировка скорости
		dc = sum_T1 - sum_T1_1
		c = c*(1+dc)

	# перезапись значений с текущего временного шага на предыдущий
	T1_1 = T1.copy()
	Y1_1 = Y1.copy()

	T2_1 = T2.copy()
	Y2_1 = Y2.copy()

	max_chem_rate = max(chem_rate)

	# вывод текущего времени на экран каждый 1000-й шаг
	if int(t_cur/dt) % 10000 == 0:
		print('%f %f %f' %(t_cur, c, dc))

	# создаем новую директорию и записываем туда файлы через каждые ... временных шагов
	'''
	if (t_cur % 500 < dt):
		# скопировать файл overall и записать промежуточные файлы parameters
		file_overall_temp = shutil.copyfile(filename_overall, dirName_beta + '/overall_temp.dat')

		file_parameters_temp = open(dirName_beta + '/param_t='+str(int(t_cur)) + '.dat', 'w') # файл с параметрами
		file_parameters_temp.write("beta      time_steps   dx   x_max   c\n")
		file_parameters_temp.write('%f   %d   %f   %f   %f\n' %( beta, time_steps, dx, x_max, c ))
		file_parameters_temp.close()
	'''
	# записывать данные в файл через каждые n временных промежутков
	if int(t_cur/dt) % (int(1/dt)*200) == 0:
		for i in range(0, N):
			file_overall.write('%f   %f   %f   %f   %f\n' % (t_cur, T1[i], Y1[i], T2[i], Y2[i]))
		time_steps += 1

	# записывать скорость реации в файл через каждые n временных промежутков
	if int(t_cur/dt) % (int(1/dt)*10) == 0:
		file_chem_rate.write('%f   %f\n' % (t_cur, max_chem_rate))

	# увеличиваем временной шаг
	t_cur += dt


for i in range(0, N):
	file_overall.write('%f   %f   %f   %f   %f\n' % (t_cur, T1[i], Y1[i], T2[i], Y2[i]))

time_steps += 1


file_parameters = open(filename_parameters, 'w') # файл с параметрами
file_parameters.write("beta1      beta2      time_steps   dx   x_max   c\n")
file_parameters.write('%f   %f   %d   %f   %f   %f\n' %( beta1, beta2, time_steps, dx, x_max, c ))
file_parameters.close()

# открыть файлы для записи
file_T1 = open(filename_T1, 'w')
file_Y1 = open(filename_Y1, 'w')

file_T2 = open(filename_T2, 'w')
file_Y2 = open(filename_Y2, 'w')

file_time = open(filename_time, 'w')
file_x = open(filename_x, 'w')

file_time.write('%f' %t_max)
file_x.write('%f' %x_max)
# записать в файл
for i in range(0, N):
	file_T1.write('%f\n' %T1[i])
	file_Y1.write('%f\n' %Y1[i])

	file_T2.write('%f\n' %T2[i])
	file_Y2.write('%f\n' %Y2[i])


file_T1.close(); file_Y1.close()
file_T2.close(); file_Y2.close()
file_overall.close(); file_time.close(); file_x.close()
file_chem_rate.close()

print('the end')