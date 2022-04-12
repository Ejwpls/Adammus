import numpy as np
import matplotlib.pyplot as plt

theta = np.linspace(0, np.pi, 10)
radius = 1000
xall = radius * np.cos(theta)
yall = radius * np.sin(theta)

xdelta = np.diff(xall)
ydelta = np.diff(yall)

x = xall[:-1]
y = yall[:-1]

xstar = (np.random.rand(100,1)-.5)*2000
ystar = (np.random.rand(100,1))*2000
# xstar = np.atleast_2d([900, -300, 500]).T
# ystar = np.atleast_2d([400, 1000, 200]).T

det = (y * xstar - x * ystar) / (ystar * xdelta - xstar * ydelta)

mask = np.logical_and(det >= 0, det <= 1)

xint = np.sum(x * mask, axis=1) + np.sum(xdelta * mask, axis=1) * det[mask]
yint = np.sum(y * mask, axis=1) + np.sum(ydelta * mask, axis=1) * det[mask]

print(xint)
print(yint)

plt.scatter(xint, yint)
plt.plot(xall, yall)
for i, j in zip(xstar, ystar):
    plt.plot(np.append(0, i), np.append(0, j))
plt.show()
