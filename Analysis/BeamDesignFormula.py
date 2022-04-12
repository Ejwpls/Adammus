import numpy as np
import matplotlib.pyplot as plt
import sympy as sp
import pandas as pd

def Beam_Simply(l,q):
    x = np.linspace(0, l, 100)  # split 'l' into 100 sections
    M = q / 2 * (l * x - x ** 2)
    V = q * (l / 2 - x)

    print("Max Shear: ", np.max(V), " kN")
    print("Max Moment: ", np.max(M), " kNm")

    plt.figure(figsize=(10, 4))
    plt.plot([0] * (l + 1), color='k')
    plt.plot(x, -M)
    plt.plot(x, V)

    plt.show()
    return

#Beam_Simply(5,20)

def solve_beam(l1, l2, q1, q2):
    l = l1 + l2  # total length
    Mx = sp.symbols('Mx')  # create symbol Mx

    # calculate Mx
    Mx = sp.solveset(Mx * l1 / 3 + q1 * l1 ** 3 / 24 + Mx * l2 / 3 + q2 * l2 ** 3 / 24, Mx).args[0]

    # sove equilibrium equations
    Va, Vb1, Vb2, Vc = sp.symbols('Va Vb1 Vb2 Vc')
    Va, Vb1 = sp.linsolve([Va + Vb1 - q1 * l1, Vb1 * l1 + Mx - (q1 * l1 ** 2) / 2],
                          (Va, Vb1)).args[0]
    Vc, Vb2 = sp.linsolve([Vb2 + Vc - q2 * l2, Vb2 * l2 + Mx - (q2 * l2 ** 2) / 2],
                          (Vc, Vb2)).args[0]
    Vb = Vb1 + Vb2

    x1 = np.arange(0, l1 + 0.1, 0.1)  # create axis x1
    x2 = np.arange(0, l2 + 0.1, 0.1)  # create axis x2

    beam1 = pd.DataFrame({"x": x1})  # create a dataframe for the first span
    beam2 = pd.DataFrame({"x": x2})  # create a dataframe for the second span

    beam1["M"] = Va * beam1.x - (q1 * beam1.x ** 2) / 2  # calculate M and store it
    beam2["M"] = Mx - (q2 * beam2.x ** 2) / 2 + Vb2 * beam2.x  # calculate M and store it

    beam1["V"] = Va - q1 * beam1.x  # calculate V and store it
    beam2["V"] = Vb2 - q2 * beam2.x  # calculate V and store it

    beam2.x = beam2.x + l1  # re-assign x for the second span

    beam = pd.concat([beam1, beam2])  # concatenate the two dataframes

    return (beam)  # return the result


header = pd.MultiIndex.from_tuples([("combo 1", "M"), ("combo 1", "V"),
                                    ("combo 2", "M"), ("combo 2", "V")])
combos = pd.DataFrame(columns=header)
combos["x"] = solve_beam(4, 5, 3.2, 4.5)["x"]

combos["combo 1"] = solve_beam(4, 5, 3.2, 4.5)
combos["combo 2"] = solve_beam(4, 5, 4.5, 3.2)
combos = combos.set_index("x")

combos = combos.astype("float")
combos.head()

fig = plt.figure(figsize=(8, 8))
ax = plt.subplot(211)
ax.invert_yaxis()

combos.loc[:, pd.IndexSlice[:, "M"]].plot(ax=ax)

ax = plt.subplot(212)
ax.invert_yaxis()
combos.loc[:, pd.IndexSlice[:, "V"]].plot(ax=ax)

plt.show()