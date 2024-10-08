import math
from dataclasses import dataclass
from scipy.optimize import minimize

@dataclass
class FoundationParams:
    ftg_cover: float
    col_length: float
    col_width: float
    load_dead: float
    load_live: float
    bp: float
    fc: float = 25
    barL: float = 0.012
    ctsL: float = 0.1
    barW: float = 0.012
    ctsW: float = 0.1

def calculate_SWt(x, params):
    l, w, d = x
    return 6 * d * l * w

def calculate_Pult(x, params):
    return 1.2 * (params.load_dead + calculate_SWt(x, params)) + 1.5 * params.load_live

def calculate_BPmax(x, params):
    l, w, d = x
    return (params.load_dead + calculate_SWt(x, params) + params.load_live) / (l * w)

def calculate_BPult(x, params):
    l, w, d = x
    return calculate_Pult(x, params) / (l * w)

def calculate_AstL(x, params):
    return 250000 / params.ctsL * params.barL**2 * math.pi

def calculate_AstW(x, params):
    return 250000 / params.ctsW * params.barW**2 * math.pi

def calculate_dsL(x, params):
    l, w, d = x
    return d - params.ftg_cover - params.barL / 2

def calculate_dsW(x, params):
    l, w, d = x
    return d - params.ftg_cover - params.barL - params.barW / 2

def calculate_AstminL(x, params):
    l, w, d = x
    return (228 * d**2 * math.sqrt(params.fc)) / calculate_dsL(x, params)

def calculate_AstminW(x, params):
    l, w, d = x
    return (228 * d**2 * math.sqrt(params.fc)) / calculate_dsW(x, params)

def calculate_alpha(params):
    return 0.85 - 0.0015 * params.fc

def calculate_gamma(params):
    return 0.97 - 0.0025 * params.fc

def calculate_MultL(x, params):
    l, w, d = x
    return ((7 * params.col_length - 10 * l) ** 2 * (-9 * calculate_SWt(x, params) + 10 * calculate_BPult(x, params) * l * w)) / (8000 * l * w)

def calculate_MultW(x, params):
    l, w, d = x
    return ((7 * params.col_width - 10 * w) ** 2 * (-9 * calculate_SWt(x, params) + 10 * calculate_BPult(x, params) * l * w)) / (8000 * l * w)

def calculate_AstshrL(x, params):
    return 5 * calculate_MultL(x, params) / (2 * calculate_dsL(x, params))

def calculate_AstshrW(x, params):
    return 5 * calculate_MultW(x, params) / (2 * calculate_dsW(x, params))

def calculate_AstreqL(x, params):
    return max(calculate_AstminL(x, params), calculate_AstshrL(x, params))

def calculate_AstreqW(x, params):
    return max(calculate_AstminW(x, params), calculate_AstshrW(x, params))

def calculate_kuL(x, params):
    l, w, d = x
    return calculate_AstreqL(x, params) / (2000 * calculate_alpha(params) * calculate_dsL(x, params) * params.fc * calculate_gamma(params) * l)

def calculate_kuW(x, params):
    l, w, d = x
    return calculate_AstreqW(x, params) / (2000 * calculate_alpha(params) * calculate_dsW(x, params) * params.fc * calculate_gamma(params) * w)

def calculate_phiL(x, params):
    return min(0.85, max(0.65, 1.24 - 13 * calculate_kuL(x, params) / 12))

def calculate_phiW(x, params):
    return min(0.85, max(0.65, 1.24 - 13 * calculate_kuW(x, params) / 12))

def calculate_fMuoL(x, params):
    return (calculate_AstL(x, params) * calculate_dsL(x, params) * calculate_phiL(x, params)) / 2 - (calculate_AstL(x, params)**2 * calculate_phiL(x, params)) / (8000 * calculate_alpha(params) * params.fc)

def calculate_fMuoW(x, params):
    return (calculate_AstW(x, params) * calculate_dsW(x, params) * calculate_phiW(x, params)) / 2 - (calculate_AstW(x, params)**2 * calculate_phiW(x, params)) / (8000 * calculate_alpha(params) * params.fc)

def calculate_CLR(x, params):
    return calculate_BPult(x, params) * (params.col_width + calculate_dsL(x, params)) * (params.col_length + calculate_dsL(x, params))

def calculate_VPult(x, params):
    return calculate_Pult(x, params) - calculate_CLR(x, params)

def calculate_fcv(params):
    return 0.17 * (1 + 2 / max(params.col_length / params.col_width, params.col_width / params.col_length, 2)) * math.sqrt(params.fc)

def calculate_fVP(x, params):
    return 1400 * calculate_dsW(x, params) * (params.col_length + params.col_width + 2 * calculate_dsW(x, params)) * calculate_fcv(params)

def calculate_dvL(x, params):
    l, w, d = x
    return max(0.9 * calculate_dsL(x, params), 0.72 * d)

def calculate_dvW(x, params):
    l, w, d = x
    return max(0.9 * calculate_dsW(x, params), 0.72 * d)

def calculate_VOultL(x, params):
    l, w, d = x
    return 0.5 * (-params.col_length - 2 * calculate_dvW(x, params) + l) * (calculate_BPult(x, params) - (9 * calculate_SWt(x, params)) / (10 * l * w))

def calculate_VOultW(x, params):
    l, w, d = x
    return 0.5 * (calculate_BPult(x, params) - (9 * calculate_SWt(x, params)) / (10 * l * w)) * (-params.col_width - 2 * calculate_dvW(x, params) + w)

def calculate_MOultL(x, params):
    l, w, d = x
    return ((params.col_length + 2 * calculate_dvW(x, params) - l) ** 2 * (-9 * calculate_SWt(x, params) + 10 * calculate_BPult(x, params) * l * w)) / (80 * l * w)

def calculate_MOultW(x, params):
    l, w, d = x
    return ((params.col_width + 2 * calculate_dvW(x, params) - w) ** 2 * (-9 * calculate_SWt(x, params) + 10 * calculate_BPult(x, params) * l * w)) / (80 * l * w)

def calculate_ex1L(x, params):
    return min((max(calculate_VOultL(x, params) * calculate_dvL(x, params), calculate_MOultL(x, params)) / calculate_dvL(x, params) + calculate_VOultL(x, params)) / (2 * 200000 * calculate_AstL(x, params)) * 1000, 3 / 1000)

def calculate_ex1W(x, params):
    return min((max(calculate_VOultW(x, params) * calculate_dvW(x, params), calculate_MOultW(x, params)) / calculate_dvW(x, params) + calculate_VOultW(x, params)) / (2 * 200000 * calculate_AstW(x, params)) * 1000, 3 / 1000)

def calculate_kvL(x, params):
    return 13 / (25 * (1 + calculate_dvL(x, params)) * (1 + 1500 * calculate_ex1L(x, params)))

def calculate_kvW(x, params):
    return 13 / (25 * (1 + calculate_dvW(x, params)) * (1 + 1500 * calculate_ex1W(x, params)))

def calculate_AngleL(x, params):
    return 29 + 7000 * calculate_ex1L(x, params)

def calculate_AngleW(x, params):
    return 29 + 7000 * calculate_ex1W(x, params)

def calculate_ks(x, params):
    l, w, d = x
    return max(1 / 2, (10 / 7) * (1 - d))

def calculate_fVucL(x, params):
    return 700 * calculate_dvL(x, params) * math.sqrt(params.fc) * calculate_ks(x, params) * calculate_kvL(x, params)

def calculate_fVucW(x, params):
    return 700 * calculate_dvW(x, params) * math.sqrt(params.fc) * calculate_ks(x, params) * calculate_kvW(x, params)

def calculate_Bpr(x, params):
    return calculate_BPmax(x, params) / params.bp

def calculate_Mur(x, params):
    return max(calculate_MultL(x, params) / calculate_fMuoL(x, params), calculate_MultW(x, params) / calculate_fMuoW(x, params))

def calculate_VPr(x, params):
    return calculate_VPult(x, params) / calculate_fVP(x, params)

def calculate_VOr(x, params):
    return max(calculate_VOultL(x, params) / calculate_fVucL(x, params), calculate_VOultW(x, params) / calculate_fVucW(x, params))

def calculate_Cost(x, params):
    l, w, d = x
    return (calculate_AstW(x, params) / 1000000 * l + calculate_AstL(x, params) / 1000000 * w) * 7850 * 3.400 + l * w * d * (130.866 * math.exp(params.fc * 0.0111) + 45 + 130) + 2 * d * l * w * 180

def objective_function(x, params):
    return calculate_Cost(x, params)

def optimize_foundation(params):
    initial_guess = [2, 2, 0.5]  # Initial guess for L, W, D
    bounds = [(max(params.col_length, params.col_width), 10),  # L and W between max(col_length, col_width) and 10m
              (max(params.col_length, params.col_width), 10),
              (0.3, 3)]  # D between 0.3m and 3m

    constraints = [
        {'type': 'ineq', 'fun': lambda x: 1 - calculate_Bpr(x, params)},
        {'type': 'ineq', 'fun': lambda x: 1 - calculate_Mur(x, params)},
        {'type': 'ineq', 'fun': lambda x: 1 - calculate_VPr(x, params)},
        {'type': 'ineq', 'fun': lambda x: 1 - calculate_VOr(x, params)},
        {'type': 'ineq', 'fun': lambda x: calculate_AstreqL(x, params) - calculate_AstL(x, params)},
        {'type': 'ineq', 'fun': lambda x: calculate_AstreqW(x, params) - calculate_AstW(x, params)}
    ]

    result = minimize(objective_function, initial_guess, args=(params,), method='SLSQP', bounds=bounds, constraints=constraints)

    return result.x

# Example usage:
params = FoundationParams(
    ftg_cover=0.060,
    col_length=0.600,
    col_width=0.500,
    load_dead=4000,
    load_live=1500,
    bp=250,
)

optimized_dimensions = optimize_foundation(params)
l, w, d = optimized_dimensions

print(f"Optimized Foundation L: {l:.2f}")
print(f"Optimized Foundation W: {w:.2f}")
print(f"Optimized Foundation D: {d:.2f}")
print(f"Optimized Cost: {calculate_Cost(optimized_dimensions, params):.2f}")
print(f"Final ratios:")
print(f"Bpr: {calculate_Bpr(optimized_dimensions, params):.2f}")
print(f"Mur: {calculate_Mur(optimized_dimensions, params):.2f}")
print(f"VPr: {calculate_VPr(optimized_dimensions, params):.2f}")
print(f"VOr: {calculate_VOr(optimized_dimensions, params):.2f}")
print(f"AstL/AstreqL: {calculate_AstL(optimized_dimensions, params)/calculate_AstreqL(optimized_dimensions, params):.2f}")
print(f"AstW/AstreqW: {calculate_AstW(optimized_dimensions, params)/calculate_AstreqW(optimized_dimensions, params):.2f}")