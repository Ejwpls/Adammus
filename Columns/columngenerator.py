import math
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt

np.seterr(divide='ignore')

# points of ku (it will always include ku_ub)
steps = 100

strain_max = 0.003
stress_max = 0.0025
modulus_steel = 200000
kub = 54 / 99

# assume k_phi
k_phi = 12 / 13

fsy = 500

ku1 = np.arange(0, kub, kub / (steps / 2))
ku2 = np.arange(kub, 1, (1 - kub) / (steps / 2 - 1))
ku = np.append(ku1, ku2)
ku = np.append(ku, np.asarray([1]))

# **************************************************
size_x = np.linspace(250, 500, 6, dtype=np.uint16)
size_y = np.linspace(250,1500,26,dtype = np.uint16)
bars = np.linspace(2, 10, 9, dtype=np.uint16)
barsize_lig = np.asarray([10, 12, 16], dtype=np.uint16)
barsize_main = np.asarray([16, 20, 24, 28, 32, 36], dtype=np.uint16)
fc = np.asarray([32, 40, 50, 65, 80, 100], dtype=np.uint16)
cover = np.asarray([30, 40], dtype=np.uint16)
# **************************************************

allvalues = np.array(np.meshgrid(size_x, size_y, bars, bars, barsize_lig, barsize_main, fc, cover),
                     dtype=np.uint32).T.reshape(-1, 8)
print(f"Comparing {len(allvalues)} columns")
mask_minbarcts = np.logical_or(
    np.logical_and((allvalues[:, 0] - 2 * allvalues[:, 7]) / allvalues[:, 2] > 250, allvalues[:, 2] > 2),
    np.logical_and((allvalues[:, 1] - 2 * allvalues[:, 7]) / allvalues[:, 3] > 250, allvalues[:, 3] > 2))
allvalues = allvalues[mask_minbarcts]
print(f"Comparing {len(allvalues)} columns")
tempmask = ((allvalues[:, 2] * 2 + allvalues[:, 3] * 2 - 4) * allvalues[:, 5] ** 2 * math.pi / 4) / (
        allvalues[:, 0] * allvalues[:, 1])
mask_minandmaxreo = np.logical_and(tempmask > 0.01, tempmask < 0.04)
allvalues = allvalues[mask_minandmaxreo]
print(f"Comparing {len(allvalues)} columns")
mask_minligreo = np.logical_not(np.logical_and(allvalues[:, 4] < 12, allvalues[:, 5] > 31))
allvalues = allvalues[mask_minligreo]
print(f"Comparing {len(allvalues)} columns")
# allvalues = np.unique(allvalues,axis=0)
# print(f"Comparing {len(allvalues)} columns")

# np.savetxt('text.txt', allvalues, fmt='%i', delimiter=',')

allvalues = allvalues.T

size_x = allvalues[0][np.newaxis].T
size_y = allvalues[1][np.newaxis].T
bars_x = allvalues[2][np.newaxis].T
bars_y = allvalues[3][np.newaxis].T
barsize_lig = allvalues[4][np.newaxis].T
barsize_main = allvalues[5][np.newaxis].T
fc = allvalues[6][np.newaxis].T
cover = allvalues[7][np.newaxis].T

A_s = np.asarray((2 * bars_x + 2 * bars_y - 4) * math.pi * (barsize_main / 2) ** 2)
alpha1 = np.clip(1 - 0.003 * fc, 0.72, 0.85)
alpha2 = np.maximum(0.85 - 0.0015 * fc, 0.67)
gamma = np.maximum(0.97 - 0.0025 * fc, 0.67)
N_uo = alpha1 * fc * (size_x * size_y - A_s) + A_s * fsy
N_uot = A_s * fsy

test = np.hstack(
    (size_x, size_y, bars_x, bars_y, barsize_lig, barsize_main, fc, A_s, alpha1, alpha2, gamma, N_uo, N_uot))
maindata = pd.DataFrame(test, columns=['sizex', 'sizey', 'barsx', 'barsy', 'sizelig', 'sizebar', 'fc', 'As', 'alpha1',
                                       'alpha2', 'gamma', 'Nuo', 'Nuot'])
print(maindata.head())
del test


def axis(size_a, size_b, bars_a, bars_b):
    df = size_a - cover - barsize_lig - barsize_main / 2
    Fcc = alpha2 * fc * size_b * gamma * df * ku
    yc = gamma * df / 2 * ku

    # bar layout
    layout = np.ones([len(str(bars_a)), np.max(bars_a)]) * 2
    layout[:, 0] = bars_b.T
    idx = np.arange(0, len(bars_a)), (bars_a - 1).T
    layout[idx] = bars_b.T
    temp = (bars_a - 1)
    idx_end = temp < np.arange(layout.shape[1])
    layout[idx_end] = 0

    # dy for bar spacing
    dy = np.ones([len(str(bars_a)), np.max(bars_a)])
    dy = np.cumsum(dy, axis=1) - 1
    dy_inc = (size_a - 2 * cover - 2 * barsize_lig - barsize_main) / (bars_a - 1)
    dy = (dy_inc * dy) + (cover + barsize_lig + barsize_main / 2)
    dy[idx_end] = np.nan

    # strain per dy
    strain_calc = ku * df

    # subtract an array (dy) from each element in temp
    es = 0.003 * (1 - dy[:, :, np.newaxis] / strain_calc[:, np.newaxis])

    # strain
    os = np.clip(es * 200000, -fsy, fsy)

    # stress
    fs_calc1 = layout * barsize_main * barsize_main / 4 * math.pi
    fs_calcpos = os - (fc * alpha2)[:, np.newaxis]
    fs_calcpos = fs_calcpos * fs_calc1[:, :, np.newaxis]
    fs_calcneg = os * fs_calc1[:, :, np.newaxis]

    fs = np.where(os > 0, fs_calcpos, fs_calcneg)
    fs = np.where(np.isfinite(fs), fs, 0)

    # nu (axial load)
    nu = np.sum(np.transpose(fs, (2, 0, 1)), axis=2).T
    nu = nu + Fcc

    nu_full = np.hstack((nu, N_uo))

    # mu (moment)
    mu = fs * dy[:, :, np.newaxis] / 1000
    mu = np.where(np.isfinite(mu), mu, 0)
    mu = np.sum(np.transpose(mu, (2, 0, 1)), axis=2).T

    mu = nu * size_a / 2000 - Fcc * yc / 1000 - mu

    temp = np.zeros(len(bars_a))[np.newaxis].T

    mu_full = np.hstack((mu, temp))



    return np.array([mu_full, nu_full])

#x = np.linspace(-2, 2, 17)
#X1 = np.stack((x[::2], x[::2]**2))
#X2 = np.stack((x[1::2], 4 - x[1::2]**2))
#x_int, y_int = intersect_piecewise(X1, X2)
#plt.plot(X1[0], X1[1], 'bo-', X2[0], X2[1], 'bo-', x_int, y_int, 'rs')
#plt.show()


result = axis(600, 600, 3, 3)
mu_full = result[0]
nu_full = result[1]

print(a)

#fig, ax = plt.subplots()
#ax.plot(mu_full[150]/1000, nu_full[150]/1000, color='blue')
#ax.spines['left'].set_position('zero')
#ax.spines['bottom'].set_position('zero')
#plt.show()

