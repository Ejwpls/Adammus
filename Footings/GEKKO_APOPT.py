from gekko import GEKKO
import numpy as np

# Create a GEKKO model
m = GEKKO(remote=False)

# Define discrete sets for each variable
L_values = list(np.arange(6, 6.85, 0.05))  # Example values for length
D_values = [1, 2, 3]        # Example values for depth
fc_values = [25, 32, 40, 50, 64]  # Example values for material cost factor

# Define variables using SOS (Special Ordered Sets)
L = m.sos1(L_values)
D = m.sos1(D_values)
fc = m.sos1(fc_values)

# Define constants
Cvr = 0.5  # Cover ratio as a constant
ColL = 0.8  # Column length as a constant
ColW = 0.5  # Column width as a constant
Pdl = 8000   # Dead load as a constant
Pll = 1500    # Live load as a constant
BP = 250    # Bearing pressure as a constant

# Define W as an intermediate variable
W = m.Intermediate(L - ColL + ColW)

# Define other intermediate variables
SWt = m.Intermediate(6 * D * L * W)
Pult = m.Intermediate(1.2 * (Pdl + SWt) + 1.5 * Pll)
BPmax = m.Intermediate((Pdl + SWt + Pll) / (L * W))
BPult = m.Intermediate(Pult / (L * W))
BPr = m.Intermediate(BPmax / BP)

# Define the objective function (e.g., minimize cost)
m.Minimize(L * W * D * fc)

# Add constraints
m.Equation(BPr <= 1)  # BPr should be less than 1
m.Equation(L >= ColL)  # Foundation length should be greater than or equal to column length
m.Equation(W >= ColW)  # Foundation width should be greater than or equal to column width

# Set solver options to use APOPT
m.options.SOLVER = 1

# Solve the optimization problem
m.solve(disp=True)

# Display results
print(f"L: {L.VALUE[0]}")
print(f"W: {W.VALUE[0]}")
print(f"D: {D.VALUE[0]}")
print(f"fc: {fc.VALUE[0]}")
print(f"BPr: {BPr.VALUE[0]}")