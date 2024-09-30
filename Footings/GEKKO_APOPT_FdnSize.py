from gekko import GEKKO
import numpy as np

# Create a GEKKO model
m = GEKKO(remote=False)

# Define discrete sets for each variable
L_values = [6350,6400,6450,6500,6550,6600]  # Example values for length
D_values = [1150,1200,1250,1300,1350,1400,1450,1500]       # Example values for depth
fc_values = [25, 32, 40, 50, 65]  # Example values for material cost factor
cts_values = [100,120,125,140,150,160,175,180,200,220,225,240,250,260,275,280,300]
bar_values = [12,16,20,24,28,32,36,40]

# Define variables using SOS (Special Ordered Sets)
L_m = m.sos1(L_values)
D_m = m.sos1(D_values)
fc = m.sos1(fc_values)
ctsL_m = m.sos1(cts_values)
ctsW_m = m.sos1(cts_values)
barL_m = m.sos1(bar_values)
# barW_m = m.sos1(bar_values)

L_m.VALUE = L_values[0]
D_m.VALUE = D_values[0]
fc.VALUE = fc_values[0]
ctsL_m.VALUE = cts_values[0]
ctsW_m.VALUE = cts_values[0]
barL_m.VALUE = bar_values[0]

#Rescale the units
L = m.Intermediate(L_m/1000)
D = m.Intermediate(D_m/1000)
ctsL = m.Intermediate(ctsL_m/1000)
ctsW = m.Intermediate(ctsW_m/1000)
barL = m.Intermediate(barL_m/1000)
barW = m.Intermediate(barL_m/1000)

# Define constants
Cvr = 0.05  # Cover ratio as a constant
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


AstL = m.Intermediate(250000 / ctsL * barL**2 * 3.141592653589793)
AstW = m.Intermediate(250000 / ctsW * barW**2 * 3.141592653589793)

dsL = m.Intermediate(D - Cvr - barL / 2)
dsW = m.Intermediate(D - Cvr - barL - barW / 2)

AstminL = m.Intermediate((228 * D**2 * m.sqrt(fc)) / dsL)
AstminW = m.Intermediate((228 * D**2 * m.sqrt(fc)) / dsW)

alpha = m.Intermediate(0.85 - 0.0015 * fc)
gamma = m.Intermediate(0.97 - 0.0025 * fc)

MultL = m.Intermediate(((7 * ColL - 10 * L) ** 2 * (-9 * SWt + 10 * BPult * L * W)) / (8000 * L * W))
MultW = m.Intermediate(((7 * ColW - 10 * W) ** 2 * (-9 * SWt + 10 * BPult * L * W)) / (8000 * L * W))

AstshrL = m.Intermediate(100 * MultL / (dsL*(50-9*gamma)))
AstshrW = m.Intermediate(100 * MultW / (dsW*(50-9*gamma)))

AstreqL = m.Intermediate(m.max2(AstminL, AstshrL))
AstreqW = m.Intermediate(m.max2(AstminW, AstshrW))

kuL = m.Intermediate(AstreqL / (2000 * alpha * dsL * fc * gamma * L))
kuW = m.Intermediate(AstreqW / (2000 * alpha * dsW * fc * gamma * W))

# phiL = phiL = m.Intermediate(0.65 + (0.85 - 0.65) * m.tanh(100 * ((1.24 - 13 * kuL / 12) - 0.65)) / 2 + 0.5)
# phiW = m.Intermediate(m.max2(0.65, m.min2(0.85, 1.24 - 13 * kuW / 12)))

phiL = m.Intermediate(0.65 + (0.85 - 0.65) * m.tanh(100 * ((1.24 - 13 * kuL / 12) - 0.65)) / 2 + 0.5)
phiW = m.Intermediate(0.65 + (0.85 - 0.65) * m.tanh(100 * ((1.24 - 13 * kuW / 12) - 0.65)) / 2 + 0.5)

fMuoL = m.Intermediate((AstL * dsL * phiL) / 2 - (AstL**2 * phiL) / (8000 * alpha * fc))
fMuoW = m.Intermediate((AstW * dsW * phiW) / 2 - (AstW**2 * phiW) / (8000 * alpha * fc))

CLR = m.Intermediate(BPult * (ColW + dsL) * (ColL + dsL))
VPult = m.Intermediate(Pult - CLR)

max_col_ratio = max(ColL / ColW, ColW / ColL, 2)
fcv = m.Intermediate(0.17 * (1 + 2 / max_col_ratio) * m.sqrt(fc))
fVP = m.Intermediate(1400 * dsW * (ColL + ColW + 2 * dsW) * fcv)

VPr = m.Intermediate(VPult / fVP)

dvL = m.Intermediate(m.max2(0.9 * dsL, 0.72 * D))
dvW = m.Intermediate(m.max2(0.9 * dsW, 0.72 * D))

VOultL = m.Intermediate(0.5 * (-ColL - 2 * dvW + L) * (BPult - (9 * SWt) / (10 * L * W)))
VOultW = m.Intermediate(0.5 * (BPult - (9 * SWt) / (10 * L * W)) * (-ColW - 2 * dvW + W))

MOultL = m.Intermediate(((ColL + 2 * dvW - L) ** 2 * (-9 * SWt + 10 * BPult * L * W)) / (80 * L * W))
MOultW = m.Intermediate(((ColW + 2 * dvW - W) ** 2 * (-9 * SWt + 10 * BPult * L * W)) / (80 * L * W))

Mur = m.Intermediate(m.max2(MultL / fMuoL, MultW / fMuoW))

ex1L = m.Intermediate(m.min2((m.max2(VOultL * dvL, MOultL) / dvL + VOultL) / (2 * 200000 * AstL) * 1000, 3 / 1000))
ex1W = m.Intermediate(m.min2((m.max2(VOultW * dvW, MOultW) / dvW + VOultW) / (2 * 200000 * AstW) * 1000, 3 / 1000))

kvL = m.Intermediate(13 / (25 * (1 + dvL) * (1 + 1500 * ex1L)))
kvW = m.Intermediate(13 / (25 * (1 + dvW) * (1 + 1500 * ex1W)))

# AngleL = m.Intermediate(29 + 7000 * ex1L)
# AngleW = m.Intermediate(29 + 7000 * ex1W)

ks = m.Intermediate(m.max2(1 / 2, (10 / 7) * (1 - D)))

fVucL = m.Intermediate(700 * dvL * m.sqrt(fc) * ks * kvL)
fVucW = m.Intermediate(700 * dvW * m.sqrt(fc) * ks * kvW)

VOr = m.Intermediate(m.max2(VOultL / fVucL, VOultW / fVucW))

Cost = m.Intermediate((AstW / 1000000 * L + AstL / 1000000 * W) * 7850 * 3.400 + L * W * D * (130.866 * m.exp(fc * 0.0111) + 45 + 130) + 2 * D * L * W * 180)
# Cost = m.Intermediate(fc*L*W*D*barL*barW/ctsL)

# Objective function (minimize cost)
m.Minimize((AstW / 1000000 * L + AstL / 1000000 * W) * 7850 * 3.400 + L * W * D * (130.866 * m.exp(fc * 0.0111) + 45 + 130) + 2 * D * L * W * 180)

# Add constraints
m.Equation(BPmax < BP)
m.Equation(VOultL < fVucL)
m.Equation(VOultW < fVucW)
m.Equation(AstL > AstreqL)
m.Equation(AstW > AstreqW)
m.Equation(m.abs(ctsL*1000 - ctsW*1000) < 50)
m.Equation(MultL < fMuoL)
m.Equation(MultW < fMuoW)
m.Equation(VPult < fVP)



# Set solver options to use APOPT
m.options.SOLVER = 1
# m.solver_options = [
#     'minlp_maximum_iterations 10000',  # Increase maximum iterations
#     'minlp_max_iter_with_int_sol 1000']  # Increase iterations with integer solutions
    # 'minlp_gap_tol 0.01']

# Solve the optimization problem
m.solve(disp=False)

# Display results
print(f"L: {L.VALUE[0]}")
print(f"W: {W.VALUE[0]}")
print(f"D: {D.VALUE[0]}")
print(f"fc: {fc.VALUE[0]}")
print(f"barL: {barL.VALUE[0]}")
print(f"ctsL: {ctsL.VALUE[0]}")
print(f"barW: {barW.VALUE[0]}")
print(f"ctsW: {ctsW.VALUE[0]}")
print("**********************")
print(f"SWt: {SWt.VALUE[0]}")
print(f"Pult: {Pult.VALUE[0]}")
print(f"BPmax: {BPmax.VALUE[0]}")
print(f"BPult: {BPult.VALUE[0]}")
print(f"BPr: {BPr.VALUE[0]}")
print(f"AstL: {AstL.VALUE[0]}")
print(f"AstW: {AstW.VALUE[0]}")
print(f"dsL: {dsL.VALUE[0]}")
print(f"dsW: {dsW.VALUE[0]}")
print(f"AstminL: {AstminL.VALUE[0]}")
print(f"AstminW: {AstminW.VALUE[0]}")
print(f"alpha: {alpha.VALUE[0]}")
print(f"gamma: {gamma.VALUE[0]}")
print(f"MultL: {MultL.VALUE[0]}")
print(f"MultW: {MultW.VALUE[0]}")
print(f"AstshrL: {AstshrL.VALUE[0]}")
print(f"AstshrW: {AstshrW.VALUE[0]}")
print(f"AstreqL: {AstreqL.VALUE[0]}")
print(f"AstreqW: {AstreqW.VALUE[0]}")
print(f"kuL: {kuL.VALUE[0]}")
print(f"kuW: {kuW.VALUE[0]}")
print(f"phiL: {phiL.VALUE[0]}")
print(f"phiW: {phiW.VALUE[0]}")
print(f"fMuoL: {fMuoL.VALUE[0]}")
print(f"fMuoW: {fMuoW.VALUE[0]}")
print(f"CLR: {CLR.VALUE[0]}")
print(f"VPult: {VPult.VALUE[0]}")
print(f"fcv: {fcv.VALUE[0]}")
print(f"fVP: {fVP.VALUE[0]}")
print(f"dvL: {dvL.VALUE[0]}")
print(f"dvW: {dvW.VALUE[0]}")
print(f"VOultL: {VOultL.VALUE[0]}")
print(f"VOultW: {VOultW.VALUE[0]}")
print(f"MOultL: {MOultL.VALUE[0]}")
print(f"MOultW: {MOultW.VALUE[0]}")
print(f"ex1L: {ex1L.VALUE[0]}")
print(f"ex1W: {ex1W.VALUE[0]}")
print(f"kvL: {kvL.VALUE[0]}")
print(f"kvW: {kvW.VALUE[0]}")
print(f"ks: {ks.VALUE[0]}")
print(f"fVucL: {fVucL.VALUE[0]}")
print(f"fVucW: {fVucW.VALUE[0]}")
print(f"Mur: {Mur.VALUE[0]}")
print(f"VPr: {VPr.VALUE[0]}")
print(f"VOr: {VOr.VALUE[0]}")
print(f"Cost: {Cost.VALUE[0]}")
