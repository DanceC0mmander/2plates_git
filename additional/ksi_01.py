import numpy as np
import matplotlib.pyplot as plt

# стационарные значения
stat = [
		]

# нестационарные значения
nonstat = [ 
		]
# значения, которые сейчас считаются
temp = [[6, 6] 
		]

# вынимаем из массива stat все x и y
# рисуем стационарные значения
x = []
y = []
for row in stat:
	x.append(row[0])
	y.append(row[1])

fig, ax = plt.subplots()
ax.plot(x, y, 'ko')

x.clear()
y.clear()

# рисуем нестационарные значения
for row in nonstat:
	x.append(row[0])
	y.append(row[1])

ax.plot(x, y, 'ro')

x.clear()
y.clear()

# рисуем временные значения
for row in temp:
	x.append(row[0])
	y.append(row[1])
ax.plot(x, y, 'bx')

# подписываем оси
plt.xlabel('$Z_{1}$')
plt.ylabel('$Z_{2}$')

# рисуем сетку
ticks = np.arange(0, 15, 0.5)
ax.set_xticks(ticks)
ax.set_yticks(ticks)
ax.grid(which='both')

# plt.savefig(folder + 'T_t='+str(int(time[0]))+'-'+str(int(time[-1]))+'.jpg', format='jpg', dpi=1000)

plt.show()