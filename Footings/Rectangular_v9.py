import math
from dataclasses import dataclass
from scipy.optimize import differential_evolution
import dataclasses

@dataclass
class FoundationParams:
    ftg_cover: float
    col_length: float
    col_width: float
    load_dead: float
    load_live: float
    bp: float
    fc: float = 25
    barL: float = 0.016
    ctsL: float = 0.2
    barW: float = 0.016
    ctsW: float = 0.2

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
    if len(x) == 3:
        l, w, d = x
        barL, ctsL = params.barL, params.ctsL
    else:
        l, w, d, barL, ctsL, barW, ctsW = x
    return 250000 / ctsL * barL**2 * math.pi

def calculate_AstW(x, params):
    if len(x) == 3:
        l, w, d = x
        barW, ctsW = params.barW, params.ctsW
    else:
        l, w, d, barL, ctsL, barW, ctsW = x
    return 250000 / ctsW * barW**2 * math.pi

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
    return 100 * calculate_MultL(x, params) / (calculate_dsL(x, params)*(50-9*calculate_gamma(params)))

def calculate_AstshrW(x, params):
    return 100 * calculate_MultW(x, params) / (calculate_dsW(x, params)*(50-9*calculate_gamma(params)))

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

def calculate_ReoRatio(x, params):
    return max((calculate_AstL(x, params) / calculate_AstreqL(x, params), calculate_AstW(x, params) / calculate_AstreqW(x, params)))

def objective_function(x, params):
    return calculate_Cost(x, params)

def optimize_foundation(params):
    bounds = [
        (max(params.col_length, params.col_width) * 1.2, 10),  # L and W
        (max(params.col_length, params.col_width) * 1.2, 10),
        (0.4, 1.5),  # D between 0.4m and 1.5m
        (0.016, 0.025),  # barL and barW between 16mm and 25mm
        (0.15, 0.3),  # ctsL and ctsW between 150mm and 300mm
        (0.016, 0.025),
        (0.15, 0.3)
    ]

    def objective_wrapper(x):
        l, w, d, barL, ctsL, barW, ctsW = x
        modified_params = dataclasses.replace(params, barL=barL, ctsL=ctsL, barW=barW, ctsW=ctsW)
        cost = calculate_Cost((l, w, d), modified_params)
        
        # Calculate ratios
        bpr = calculate_Bpr((l, w, d), modified_params)
        mur = calculate_Mur((l, w, d), modified_params)
        vpr = calculate_VPr((l, w, d), modified_params)
        vor = calculate_VOr((l, w, d), modified_params)
        reo_ratio = calculate_ReoRatio((l, w, d), modified_params)
        
        penalty = 0
        # Stricter penalties for ratios that should be between 0 and 1
        for ratio in [bpr, mur, vpr, vor, reo_ratio]:
            if ratio < 0 or ratio > 1:
                penalty += 1e6 * (max(0, -ratio) + max(0, ratio - 1))
        
        return cost + penalty

    result = differential_evolution(objective_wrapper, bounds, strategy='best1bin', 
                                    maxiter=10000, popsize=15, tol=0.01, 
                                    mutation=(0.5, 1), recombination=0.7, 
                                    seed=None, callback=None, disp=True, polish=True, 
                                    init='latinhypercube', atol=0, updating='immediate', 
                                    workers=1)

    print("Optimization successful:", result.success)
    print("Optimization message:", result.message)

    return result.x

# Example usage:
# params = FoundationParams(
#     ftg_cover=0.060,
#     col_length=0.600,
#     col_width=0.500,
#     load_dead=4000,
#     load_live=1500,
#     bp=250,
# )

# optimized_values = optimize_foundation(params)
# l, w, d, barL, ctsL, barW, ctsW = optimized_values

# print(f"Optimized Foundation L: {l:.3f}")
# print(f"Optimized Foundation W: {w:.3f}")
# print(f"Optimized Foundation D: {d:.3f}")
# print(f"Optimized barL: {barL:.3f}")
# print(f"Optimized ctsL: {ctsL:.3f}")
# print(f"Optimized barW: {barW:.3f}")
# print(f"Optimized ctsW: {ctsW:.3f}")
# print(f"Optimized Cost: {calculate_Cost(optimized_values[:3], dataclasses.replace(params, barL=barL, ctsL=ctsL, barW=barW, ctsW=ctsW)):.2f}")
# print(f"Final ratios:")
# print(f"Bpr: {calculate_Bpr(optimized_values[:3], params):.2f}")
# print(f"Mur: {calculate_Mur(optimized_values[:3], params):.2f}")
# print(f"VPr: {calculate_VPr(optimized_values[:3], params):.2f}")
# print(f"VOr: {calculate_VOr(optimized_values[:3], params):.2f}")
# print(f"ReoRatio: {calculate_ReoRatio(optimized_values[:3], params):.2f}")

# # Print constraint values
# print("\nConstraint values:")
# print(f"Bpr constraint: {0.999 - calculate_Bpr(optimized_values[:3], params):.6f}")
# print(f"Mur constraint: {0.999 - calculate_Mur(optimized_values[:3], params):.6f}")
# print(f"VPr constraint: {0.999 - calculate_VPr(optimized_values[:3], params):.6f}")
# print(f"VOr constraint: {0.999 - calculate_VOr(optimized_values[:3], params):.6f}")
# print(f"ReoRatio constraint: {1.5 - calculate_ReoRatio(optimized_values[:3], params):.6f}")


def print_all_results(params, l, w, d, barL, ctsL, barW, ctsW):
    x = (l, w, d)
    modified_params = dataclasses.replace(params, barL=barL, ctsL=ctsL, barW=barW, ctsW=ctsW)
    
    print(f"Results for L={l:.3f}, W={w:.3f}, D={d:.3f}, barL={barL:.3f}, ctsL={ctsL:.3f}, barW={barW:.3f}, ctsW={ctsW:.3f}")
    print(f"SWt: {calculate_SWt(x, modified_params):.2f}")
    print(f"Pult: {calculate_Pult(x, modified_params):.2f}")
    print(f"BPmax: {calculate_BPmax(x, modified_params):.2f}")
    print(f"BPult: {calculate_BPult(x, modified_params):.2f}")
    print(f"AstL: {calculate_AstL(x, modified_params):.2f}")
    print(f"AstW: {calculate_AstW(x, modified_params):.2f}")
    print(f"dsL: {calculate_dsL(x, modified_params):.2f}")
    print(f"dsW: {calculate_dsW(x, modified_params):.2f}")
    print(f"AstminL: {calculate_AstminL(x, modified_params):.2f}")
    print(f"AstminW: {calculate_AstminW(x, modified_params):.2f}")
    print(f"alpha: {calculate_alpha(modified_params):.2f}")
    print(f"gamma: {calculate_gamma(modified_params):.2f}")
    print(f"MultL: {calculate_MultL(x, modified_params):.2f}")
    print(f"MultW: {calculate_MultW(x, modified_params):.2f}")
    print(f"AstshrL: {calculate_AstshrL(x, modified_params):.2f}")
    print(f"AstshrW: {calculate_AstshrW(x, modified_params):.2f}")
    print(f"AstreqL: {calculate_AstreqL(x, modified_params):.2f}")
    print(f"AstreqW: {calculate_AstreqW(x, modified_params):.2f}")
    print(f"kuL: {calculate_kuL(x, modified_params):.2f}")
    print(f"kuW: {calculate_kuW(x, modified_params):.2f}")
    print(f"phiL: {calculate_phiL(x, modified_params):.2f}")
    print(f"phiW: {calculate_phiW(x, modified_params):.2f}")
    print(f"fMuoL: {calculate_fMuoL(x, modified_params):.2f}")
    print(f"fMuoW: {calculate_fMuoW(x, modified_params):.2f}")
    print(f"CLR: {calculate_CLR(x, modified_params):.2f}")
    print(f"VPult: {calculate_VPult(x, modified_params):.2f}")
    print(f"fcv: {calculate_fcv(modified_params):.2f}")
    print(f"fVP: {calculate_fVP(x, modified_params):.2f}")
    print(f"dvL: {calculate_dvL(x, modified_params):.2f}")
    print(f"dvW: {calculate_dvW(x, modified_params):.2f}")
    print(f"VOultL: {calculate_VOultL(x, modified_params):.2f}")
    print(f"VOultW: {calculate_VOultW(x, modified_params):.2f}")
    print(f"MOultL: {calculate_MOultL(x, modified_params):.2f}")
    print(f"MOultW: {calculate_MOultW(x, modified_params):.2f}")
    print(f"ex1L: {calculate_ex1L(x, modified_params):.6f}")
    print(f"ex1W: {calculate_ex1W(x, modified_params):.6f}")
    print(f"kvL: {calculate_kvL(x, modified_params):.2f}")
    print(f"kvW: {calculate_kvW(x, modified_params):.2f}")
    print(f"AngleL: {calculate_AngleL(x, modified_params):.2f}")
    print(f"AngleW: {calculate_AngleW(x, modified_params):.2f}")
    print(f"ks: {calculate_ks(x, modified_params):.2f}")
    print(f"fVucL: {calculate_fVucL(x, modified_params):.2f}")
    print(f"fVucW: {calculate_fVucW(x, modified_params):.2f}")
    print(f"Bpr: {calculate_Bpr(x, modified_params):.2f}")
    print(f"Mur: {calculate_Mur(x, modified_params):.2f}")
    print(f"VPr: {calculate_VPr(x, modified_params):.2f}")
    print(f"VOr: {calculate_VOr(x, modified_params):.2f}")
    print(f"Cost: {calculate_Cost(x, modified_params):.2f}")

# Example usage:
params = FoundationParams(
    ftg_cover=0.050,
    col_length=0.800,
    col_width=0.500,
    load_dead=8000,
    load_live=1500,
    bp=250,
)

# Input specific values
l = 6.45
w = 6.15
d = 1.7
barL = 0.040
ctsL = 0.1
barW = 0.040
ctsW = 0.1

print_all_results(params, l, w, d, barL, ctsL, barW, ctsW)