import numpy as np
import numexpr as ne
import math
import dataclasses as dc

def create_concrete_combos(start_L, end_L, start_D, end_D, step):
    # Create array A (L values)
    L = np.arange(start_L, end_L + 0.05, 0.05)
    
    # Create array B (D values)
    D = np.arange(start_D, end_D + 0.05, 0.05)
    
    # Create array C (fc values)
    fc_values = np.array([25, 32, 40, 50, 65])
    
    # Create all combinations using numpy's meshgrid
    L_mesh, D_mesh, fc_mesh = np.meshgrid(L, D, fc_values, indexing='ij')
        
    # Create W_mesh by adding step to L_mesh
    W_mesh = L_mesh + step
    
    # Combine the meshes into the final array
    combinations = np.column_stack((L_mesh.ravel(), W_mesh.ravel(), D_mesh.ravel(), fc_mesh.ravel()))
    
    # Calculate A*B*C*D using numexpr
    L_vals = combinations[:, 0]
    W_vals = combinations[:, 1]
    D_vals = combinations[:, 2]
    fc_vals = combinations[:, 3]
    product_ABCD = ne.evaluate('L_vals * W_vals * D_vals * fc_vals')
    
    # Sort combinations based on A*B*C*D
    sorted_indices = np.argsort(product_ABCD)
    sorted_combinations = combinations[sorted_indices]
    
    return sorted_combinations

def create_steel_combos():
    # Create array A
    A = np.array([0.012,0.016,0.020,0.024,0.028,0.032,0.036,0.040])
    
    # Create array B
    B = np.array([0.1, 0.120, 0.125, 0.14, 0.15, 0.16, 0.175, 0.18, 0.2, 0.22, 0.225, 0.24, 0.25, 0.26, 0.275, 0.28, 0.3])
    

    # Create all combinations using numpy's meshgrid
    A_mesh1, B_mesh1, A_mesh2, B_mesh2 = np.meshgrid(A, B, A, B, indexing='ij')
    combinations = np.stack([A_mesh1.ravel(), B_mesh1.ravel(), A_mesh2.ravel(), B_mesh2.ravel()], axis=-1)
    
    # Calculate A*B*A*B using numexpr
    A_vals1 = combinations[:, 0]
    B_vals1 = combinations[:, 1]
    A_vals2 = combinations[:, 2]
    B_vals2 = combinations[:, 3]
    product_ABAB = ne.evaluate('A_vals1 * B_vals1 * A_vals2 * B_vals2')
    
    # Sort combinations based on A*B*A*B
    sorted_indices = np.argsort(product_ABAB)
    sorted_combinations = combinations[sorted_indices]
    
    return sorted_combinations

def concrete_filter(conc_arr):
    # Extract individual columns
    A = conc_arr[:, 0]
    B = conc_arr[:, 1]
    C = conc_arr[:, 2]
    
    # Create a boolean mask for the condition using numexpr
    mask = ne.evaluate('(0 < (8000 + 1500 + 6 * A * B * C) / (250 * A * B)) & ((8000 + 1500 + 6 * A * B * C) / (250 * A * B) < 1)')
    
    # Apply the mask to filter the array
    return conc_arr[mask]

def steel_filter(conc_arr, steel_arr):
    # Get the smallest value in the last column of array1
    depth = np.min(conc_arr[:, -2])
    
    # Prepare the dictionary for numexpr
    context = {
        'barL': steel_arr[:, 0],
        'ctsL': steel_arr[:, 1],
        'barW': steel_arr[:, 2],
        'ctsW': steel_arr[:, 3],
        'D': depth,
        'pi': math.pi,
        'sqrt_25': math.sqrt(25)
    }

    # Calculate the conditions using numexpr
    condition1 = ne.evaluate('31250 * barL**2 * (-barL - 2*0.05 + 2*D) * pi / (57 * ctsL * D**2 * sqrt_25) > 1', context)
    condition2 = ne.evaluate('31250 * barW**2 * (-barW - 2*0.05 + 2*D) * pi / (57 * ctsW * D**2 * sqrt_25) > 1', context)

    # Combine conditions using logical AND
    mask = np.logical_and(condition1, condition2)
    
    # Filter array2 based on the mask
    filtered_rows = steel_arr[mask]
    
    return filtered_rows  # Return the filtered array2

@dc.dataclass
class Foundation:
    L: float = 0.0
    W: float = 0.0
    D: float = 0.0
    Cvr: float = 0.0
    ColL: float = 0.0
    ColW: float = 0.0
    Pdl: float = 0.0
    Pll: float = 0.0
    BP: float = 0.0
    fc: float = 0.0
    barL: float = 0.0
    ctsL: float = 0.0
    barW: float = 0.0
    ctsW: float = 0.0

    def __post_init__(self):
        for field in dc.fields(self):
            setattr(self, field.name, np.atleast_1d(getattr(self, field.name)))

def calculate_SWt(fdn):
    return 6 * fdn.D * fdn.L * fdn.W

def calculate_Pult(fdn):
    return 1.2 * (fdn.Pdl + calculate_SWt(fdn)) + 1.5 * fdn.Pll

def calculate_BPmax(fdn):
    return (fdn.Pdl + calculate_SWt(fdn) + fdn.Pll) / (fdn.L * fdn.W)

def calculate_BPult(fdn):
    return calculate_Pult(fdn) / (fdn.L * fdn.W)

def calculate_AstL(fdn):
    return 250000 / fdn.ctsL * fdn.barL**2 * np.pi

def calculate_AstW(fdn):
    return 250000 / fdn.ctsW * fdn.barW**2 * np.pi

def calculate_dsL(fdn):
    return fdn.D - fdn.Cvr - fdn.barL / 2

def calculate_dsW(fdn):
    return fdn.D - fdn.Cvr - fdn.barL - fdn.barW / 2

def calculate_AstminL(fdn):
    return (228 * fdn.D**2 * np.sqrt(fdn.fc)) / calculate_dsL(fdn)

def calculate_AstminW(fdn):
    return (228 * fdn.D**2 * np.sqrt(fdn.fc)) / calculate_dsW(fdn)

def calculate_alpha(fdn):
    return 0.85 - 0.0015 * fdn.fc

def calculate_gamma(fdn):
    return 0.97 - 0.0025 * fdn.fc

def calculate_MultL(fdn):
    return ((7 * fdn.ColL - 10 * fdn.L) ** 2 * (-9 * calculate_SWt(fdn) + 10 * calculate_BPult(fdn) * fdn.L * fdn.W)) / (8000 * fdn.L * fdn.W)

def calculate_MultW(fdn):
    return ((7 * fdn.ColW - 10 * fdn.W) ** 2 * (-9 * calculate_SWt(fdn) + 10 * calculate_BPult(fdn) * fdn.L * fdn.W)) / (8000 * fdn.L * fdn.W)

def calculate_AstshrL(fdn):
    return 100 * calculate_MultL(fdn) / (calculate_dsL(fdn)*(50-9*calculate_gamma(fdn)))

def calculate_AstshrW(fdn):
    return 100 * calculate_MultW(fdn) / (calculate_dsW(fdn)*(50-9*calculate_gamma(fdn)))

def calculate_AstreqL(fdn):
    return np.maximum(calculate_AstminL(fdn), calculate_AstshrL(fdn))

def calculate_AstreqW(fdn):
    return np.maximum(calculate_AstminW(fdn), calculate_AstshrW(fdn))

def calculate_kuL(fdn):
    return calculate_AstreqL(fdn) / (2000 * calculate_alpha(fdn) * calculate_dsL(fdn) * fdn.fc * calculate_gamma(fdn) * fdn.L)

def calculate_kuW(fdn):
    return calculate_AstreqW(fdn) / (2000 * calculate_alpha(fdn) * calculate_dsW(fdn) * fdn.fc * calculate_gamma(fdn) * fdn.W)

def calculate_phiL(fdn):
    return np.minimum(0.85, np.maximum(0.65, 1.24 - 13 * calculate_kuL(fdn) / 12))

def calculate_phiW(fdn):
    return np.minimum(0.85, np.maximum(0.65, 1.24 - 13 * calculate_kuW(fdn) / 12))

def calculate_fMuoL(fdn):
    return (calculate_AstL(fdn) * calculate_dsL(fdn) * calculate_phiL(fdn)) / 2 - (calculate_AstL(fdn)**2 * calculate_phiL(fdn)) / (8000 * calculate_alpha(fdn) * fdn.fc)

def calculate_fMuoW(fdn):
    return (calculate_AstW(fdn) * calculate_dsW(fdn) * calculate_phiW(fdn)) / 2 - (calculate_AstW(fdn)**2 * calculate_phiW(fdn)) / (8000 * calculate_alpha(fdn) * fdn.fc)

def calculate_CLR(fdn):
    return calculate_BPult(fdn) * (fdn.ColW + calculate_dsL(fdn)) * (fdn.ColL + calculate_dsL(fdn))

def calculate_VPult(fdn):
    return calculate_Pult(fdn) - calculate_CLR(fdn)

def calculate_fcv(fdn):
    return 0.17 * (1 + 2 / np.maximum.reduce([fdn.ColL / fdn.ColW, fdn.ColW / fdn.ColL, np.full(len(fdn.L), 2.0)])) * np.sqrt(fdn.fc)

def calculate_fVP(fdn):
    return 1400 * calculate_dsW(fdn) * (fdn.ColL + fdn.ColW + 2 * calculate_dsW(fdn)) * calculate_fcv(fdn)

def calculate_dvL(fdn):
    return np.maximum(0.9 * calculate_dsL(fdn), 0.72 * fdn.D)

def calculate_dvW(fdn):
    return np.maximum(0.9 * calculate_dsW(fdn), 0.72 * fdn.D)

def calculate_VOultL(fdn):
    return 0.5 * (-fdn.ColL - 2 * calculate_dvW(fdn) + fdn.L) * (calculate_BPult(fdn) - (9 * calculate_SWt(fdn)) / (10 * fdn.L * fdn.W))

def calculate_VOultW(fdn):
    return 0.5 * (calculate_BPult(fdn) - (9 * calculate_SWt(fdn)) / (10 * fdn.L * fdn.W)) * (-fdn.ColW - 2 * calculate_dvW(fdn) + fdn.W)

def calculate_MOultL(fdn):
    return ((fdn.ColL + 2 * calculate_dvW(fdn) - fdn.L) ** 2 * (-9 * calculate_SWt(fdn) + 10 * calculate_BPult(fdn) * fdn.L * fdn.W)) / (80 * fdn.L * fdn.W)

def calculate_MOultW(fdn):
    return ((fdn.ColW + 2 * calculate_dvW(fdn) - fdn.W) ** 2 * (-9 * calculate_SWt(fdn) + 10 * calculate_BPult(fdn) * fdn.L * fdn.W)) / (80 * fdn.L * fdn.W)

def calculate_ex1L(fdn):
    return np.minimum((np.maximum(calculate_VOultL(fdn) * calculate_dvL(fdn), calculate_MOultL(fdn)) / calculate_dvL(fdn) + calculate_VOultL(fdn)) / (2 * 200000 * calculate_AstL(fdn)) * 1000, 3 / 1000)

def calculate_ex1W(fdn):
    return np.minimum((np.maximum(calculate_VOultW(fdn) * calculate_dvW(fdn), calculate_MOultW(fdn)) / calculate_dvW(fdn) + calculate_VOultW(fdn)) / (2 * 200000 * calculate_AstW(fdn)) * 1000, 3 / 1000)

def calculate_kvL(fdn):
    return 13 / (25 * (1 + calculate_dvL(fdn)) * (1 + 1500 * calculate_ex1L(fdn)))

def calculate_kvW(fdn):
    return 13 / (25 * (1 + calculate_dvW(fdn)) * (1 + 1500 * calculate_ex1W(fdn)))

def calculate_AngleL(fdn):
    return 29 + 7000 * calculate_ex1L(fdn)

def calculate_AngleW(fdn):
    return 29 + 7000 * calculate_ex1W(fdn)

def calculate_ks(fdn):
    return np.maximum(1 / 2, (10 / 7) * (1 - fdn.D))

def calculate_fVucL(fdn):
    return 700 * calculate_dvL(fdn) * np.sqrt(fdn.fc) * calculate_ks(fdn) * calculate_kvL(fdn)

def calculate_fVucW(fdn):
    return 700 * calculate_dvW(fdn) * np.sqrt(fdn.fc) * calculate_ks(fdn) * calculate_kvW(fdn)

def calculate_Bpr(fdn):
    return calculate_BPmax(fdn) / fdn.BP

def calculate_Mur(fdn):
    return np.maximum(calculate_MultL(fdn) / calculate_fMuoL(fdn), calculate_MultW(fdn) / calculate_fMuoW(fdn))

def calculate_VPr(fdn):
    return calculate_VPult(fdn) / calculate_fVP(fdn)

def calculate_VOr(fdn):
    return np.maximum(calculate_VOultL(fdn) / calculate_fVucL(fdn), calculate_VOultW(fdn) / calculate_fVucW(fdn))

def calculate_Cost(fdn):
    return (calculate_AstW(fdn) / 1000000 * fdn.L + calculate_AstL(fdn) / 1000000 * fdn.W) * 7850 * 3.400 + fdn.L * fdn.W * fdn.D * (130.866 * np.exp(fdn.fc * 0.0111) + 45 + 130) + 2 * fdn.D * fdn.L * fdn.W * 180

def calculate_ReoRatio(fdn):
    return np.maximum(calculate_AstL(fdn) / calculate_AstreqL(fdn), calculate_AstW(fdn) / calculate_AstreqW(fdn))

def calc_starting_ranges(fdn):
    extent = np.max([fdn.L.max(), fdn.W.max()])
    D_start = np.round(extent/6 / 0.05) * 0.05-0.1
    D_end = np.round(extent/2 / 0.05) * 0.05+0.1
    D = np.arange(D_start, D_end + 0.05, 0.05)
    _length = len(D)
    L = np.full(_length, fdn.L)
    W = np.full(_length, fdn.W)
    ColL = np.full(_length, fdn.ColL)
    ColW = np.full(_length, fdn.ColW)
    BP = np.full(_length, fdn.BP)
    fc = np.full(_length, fdn.fc)
    barL = np.full(_length, fdn.barL)
    ctsL = np.full(_length, fdn.ctsL)
    barW = np.full(_length, fdn.barW)
    ctsW = np.full(_length, fdn.ctsW)
    Pdl = np.full(_length, fdn.Pdl)
    Pll = np.full(_length, fdn.Pll)
    Cvr = np.full(_length, fdn.Cvr)
    return ([L, W, D, Cvr, ColL, ColW, Pdl, Pll, BP, fc, barL, ctsL, barW, ctsW])

def calc_min_L(fdn):
    L = (fdn.BP *(fdn.ColL-fdn.ColW)+np.sqrt(fdn.BP*(fdn.BP*(fdn.ColL-fdn.ColW)**2+4 *(fdn.Pdl+fdn.Pll))))/(2* fdn.BP) 
    return (np.ceil(L/0.05)*0.05)[0]

def calc_min_W(fdn):
    return (calc_min_L(fdn)-fdn.ColL+fdn.ColW)[0]

def calc_min_D(fdn):
    params = fdn
    params.fc = 65
    params.barL=0.040
    params.ctsL=0.100
    params.barW=0.040
    params.ctsW=0.100
    params.L = calc_min_L(params)
    params.W = params.L+params.ColW-params.ColL
    temparray = calc_starting_ranges(params)
    startingfdn = Foundation(L=temparray[0],W=temparray[1],D=temparray[2],Cvr=temparray[3],ColL=temparray[4],ColW=temparray[5],Pdl=temparray[6],Pll=temparray[7],BP=temparray[8],fc=temparray[9],barL=temparray[10],ctsL=temparray[11],barW=temparray[12],ctsW=temparray[13])
    mask = calculate_fVucL(startingfdn)/calculate_VOultL(startingfdn) > 1
    return np.min(startingfdn.D[mask])

def calc_max_D(fdn):
    params = fdn
    params.fc = 25
    params.barL=0.012
    params.ctsL=0.300
    params.barW=0.012
    params.ctsW=0.300
    params.L = calc_min_L(params)
    params.W = params.L+params.ColW-params.ColL
    temparray = calc_starting_ranges(params)
    startingfdn = Foundation(L=temparray[0],W=temparray[1],D=temparray[2],Cvr=temparray[3],ColL=temparray[4],ColW=temparray[5],Pdl=temparray[6],Pll=temparray[7],BP=temparray[8],fc=temparray[9],barL=temparray[10],ctsL=temparray[11],barW=temparray[12],ctsW=temparray[13])
    mask = calculate_fVucL(startingfdn)/calculate_VOultL(startingfdn) > 1
    return np.min(startingfdn.D[mask])

def calc_max_L(fdn):
    _diff = fdn.ColW-fdn.ColL
    _temp_a = fdn.BP-6*calc_max_D(fdn)
    max_L = (_diff*_temp_a-np.sqrt(_temp_a)*np.sqrt(4*(fdn.Pdl+fdn.Pll)+_diff**2*_temp_a))/(-2*_temp_a)
    return (np.ceil(max_L/0.05)*0.05)[0]

def calc_max_W(fdn):
    return (calc_max_L(fdn)-fdn.ColL+fdn.ColW)[0]

def print_foundation_results(fdn,idx=0):
    print(f"Results for L={fdn.L[idx]:.3f}, W={fdn.W[idx]:.3f}, D={fdn.D[idx]:.3f}, barL={fdn.barL[idx]:.3f}, ctsL={fdn.ctsL[idx]:.3f}, barW={fdn.barW[idx]:.3f}, ctsW={fdn.ctsW[idx]:.3f},Cvr={ fdn.Cvr[idx]:.3f}")
    print(f"SWt: {calculate_SWt(fdn)[idx]:.2f}")
    print(f"Pult: {calculate_Pult(fdn)[idx]:.2f}")
    print(f"BPmax: {calculate_BPmax(fdn)[idx]:.2f}")
    print(f"BPult: {calculate_BPult(fdn)[idx]:.2f}")
    print(f"AstL: {calculate_AstL(fdn)[idx]:.0f}")
    print(f"AstW: {calculate_AstW(fdn)[idx]:.0f}")
    print(f"dsL: {calculate_dsL(fdn)[idx]:.3f}")
    print(f"dsW: {calculate_dsW(fdn)[idx]:.3f}")
    print(f"AstminL: {calculate_AstminL(fdn)[idx]:.0f}")
    print(f"AstminW: {calculate_AstminW(fdn)[idx]:.0f}")
    print(f"alpha: {calculate_alpha(fdn)[idx]:.4f}")
    print(f"gamma: {calculate_gamma(fdn)[idx]:.4f}")
    print(f"MultL: {calculate_MultL(fdn)[idx]:.2f}")
    print(f"MultW: {calculate_MultW(fdn)[idx]:.2f}")
    print(f"AstshrL: {calculate_AstshrL(fdn)[idx]:.0f}")
    print(f"AstshrW: {calculate_AstshrW(fdn)[idx]:.0f}")
    print(f"AstreqL: {calculate_AstreqL(fdn)[idx]:.0f}")
    print(f"AstreqW: {calculate_AstreqW(fdn)[idx]:.0f}")
    print(f"kuL: {calculate_kuL(fdn)[idx]:.2f}")
    print(f"kuW: {calculate_kuW(fdn)[idx]:.2f}")
    print(f"phiL: {calculate_phiL(fdn)[idx]:.2f}")
    print(f"phiW: {calculate_phiW(fdn)[idx]:.2f}")
    print(f"fMuoL: {calculate_fMuoL(fdn)[idx]:.2f}")
    print(f"fMuoW: {calculate_fMuoW(fdn)[idx]:.2f}")
    print(f"CLR: {calculate_CLR(fdn)[idx]:.2f}")
    print(f"VPult: {calculate_VPult(fdn)[idx]:.2f}")
    print(f"fcv: {calculate_fcv(fdn)[idx]:.2f}")
    print(f"fVP: {calculate_fVP(fdn)[idx]:.2f}")
    print(f"dvL: {calculate_dvL(fdn)[idx]:.3f}")
    print(f"dvW: {calculate_dvW(fdn)[idx]:.3f}")
    print(f"VOultL: {calculate_VOultL(fdn)[idx]:.1f}")
    print(f"VOultW: {calculate_VOultW(fdn)[idx]:.1f}")
    print(f"MOultL: {calculate_MOultL(fdn)[idx]:.1f}")
    print(f"MOultW: {calculate_MOultW(fdn)[idx]:.1f}")
    print(f"ex1L: {calculate_ex1L(fdn)[idx]:.6f}")
    print(f"ex1W: {calculate_ex1W(fdn)[idx]:.6f}")
    print(f"kvL: {calculate_kvL(fdn)[idx]:.2f}")
    print(f"kvW: {calculate_kvW(fdn)[idx]:.2f}")
    print(f"AngleL: {calculate_AngleL(fdn)[idx]:.2f}")
    print(f"AngleW: {calculate_AngleW(fdn)[idx]:.2f}")
    print(f"ks: {calculate_ks(fdn)[idx]:.2f}")
    print(f"fVucL: {calculate_fVucL(fdn)[idx]:.2f}")
    print(f"fVucW: {calculate_fVucW(fdn)[idx]:.2f}")
    print(f"Bpr: {calculate_Bpr(fdn)[idx]:.2f}")
    print(f"Mur: {calculate_Mur(fdn)[idx]:.2f}")
    print(f"VPr: {calculate_VPr(fdn)[idx]:.2f}")
    print(f"VOr: {calculate_VOr(fdn)[idx]:.2f}")
    print(f"Cost: {calculate_Cost(fdn)[idx]:.2f}")

# test = Foundation(L=6.5,W=6.2,D=2.05,ColL=0.8,ColW=0.5,Pdl=8000,Pll=1500,BP=250,fc=25,barL=0.024,ctsL=0.15,barW=0.024,ctsW=0.15,Cvr=0.05)
test = Foundation(L=6.3,W=6.0,D=2.55,ColL=0.8,ColW=0.5,Pdl=8000,Pll=1500,BP=250,fc=25,barL=0.012,ctsL=0.3,barW=0.012,ctsW=0.3,Cvr=0.05)


conc = create_concrete_combos(calc_min_L(test),calc_max_L(test),calc_min_D(test),calc_max_D(test),test.ColW-test.ColL)
conc = concrete_filter(conc)
reo = create_steel_combos()
reo = steel_filter(conc, reo)

conc_shape = conc.shape
reo_shape = reo.shape
pi = math.pi

# Reshape conc and reo for broadcasting
conc_reshaped = conc[:, np.newaxis, :]
reo_reshaped = reo[np.newaxis, :, :]

# Create arrays for each parameter
L = np.broadcast_to(conc_reshaped[:, :, 0], (conc_shape[0], reo_shape[0]))
W = np.broadcast_to(conc_reshaped[:, :, 1], (conc_shape[0], reo_shape[0]))
D = np.broadcast_to(conc_reshaped[:, :, 2], (conc_shape[0], reo_shape[0]))
fc = np.broadcast_to(conc_reshaped[:, :, 3], (conc_shape[0], reo_shape[0]))
ColL = np.full_like(L, test.ColL)
ColW = np.full_like(L, test.ColW)
Pdl = np.full_like(L, test.Pdl)
Pll = np.full_like(L, test.Pll)
BP = np.full_like(L, test.BP)
barL = np.broadcast_to(reo_reshaped[:, :, 0], (conc_shape[0], reo_shape[0]))
ctsL = np.broadcast_to(reo_reshaped[:, :, 1], (conc_shape[0], reo_shape[0]))
barW = np.broadcast_to(reo_reshaped[:, :, 2], (conc_shape[0], reo_shape[0]))
ctsW = np.broadcast_to(reo_reshaped[:, :, 3], (conc_shape[0], reo_shape[0]))
Cvr = np.full_like(L, test.Cvr)

# Use numexpr for complex calculations
SWt = ne.evaluate('6 * D * L * W')

Pult = ne.evaluate('1.2 * (Pdl + SWt) + 1.5 * Pll')
BPmax = ne.evaluate('(Pdl + SWt + Pll) / (L * W)')
BPult = ne.evaluate('Pult / (L * W)')

AstL = ne.evaluate('250000 / ctsL * barL**2 * pi')
AstW = ne.evaluate('250000 / ctsW * barW**2 * pi')

dsL = ne.evaluate('D - Cvr - barL / 2')
dsW = ne.evaluate('D - Cvr - barL - barW / 2')

AstminL = ne.evaluate('(228 * D**2 * sqrt(fc)) / dsL')
AstminW = ne.evaluate('(228 * D**2 * sqrt(fc)) / dsW')

alpha = ne.evaluate('0.85 - 0.0015 * fc')
gamma = ne.evaluate('0.97 - 0.0025 * fc')

MultL = ne.evaluate('((7 * ColL - 10 * L) ** 2 * (-9 * SWt + 10 * BPult * L * W)) / (8000 * L * W)')
MultW = ne.evaluate('((7 * ColW - 10 * W) ** 2 * (-9 * SWt + 10 * BPult * L * W)) / (8000 * L * W)')

AstshrL = ne.evaluate('100 * MultL / (dsL * (50 - 9 * gamma))')
AstshrW = ne.evaluate('100 * MultW / (dsW * (50 - 9 * gamma))')

# Use NumPy's maximum function instead of numexpr
AstreqL = np.maximum(AstminL, AstshrL)
AstreqW = np.maximum(AstminW, AstshrW)

kuL = ne.evaluate('AstreqL / (2000 * alpha * dsL * fc * gamma * L)')
kuW = ne.evaluate('AstreqW / (2000 * alpha * dsW * fc * gamma * W)')

# Use NumPy's minimum and maximum functions
phiL = np.minimum(0.85, np.maximum(0.65, ne.evaluate('1.24 - 13 * kuL / 12')))
phiW = np.minimum(0.85, np.maximum(0.65, ne.evaluate('1.24 - 13 * kuW / 12')))

fMuoL = ne.evaluate('(AstL * dsL * phiL) / 2 - (AstL**2 * phiL) / (8000 * alpha * fc)')
fMuoW = ne.evaluate('(AstW * dsW * phiW) / 2 - (AstW**2 * phiW) / (8000 * alpha * fc)')

CLR = ne.evaluate('BPult * (ColW + dsL) * (ColL + dsL)')
VPult = ne.evaluate('Pult - CLR')

# Use NumPy's maximum function for fcv calculation
fcv = ne.evaluate('0.17 * (1 + 2 / col_ratio) * sqrt(fc)', 
                  {'col_ratio': np.maximum(ColL / ColW, np.maximum(ColW / ColL, 2.0))})
fVP = ne.evaluate('1400 * dsW * (ColL + ColW + 2 * dsW) * fcv')

dvL = np.maximum(ne.evaluate('0.9 * dsL'), ne.evaluate('0.72 * D'))
dvW = np.maximum(ne.evaluate('0.9 * dsW'), ne.evaluate('0.72 * D'))

VOultL = ne.evaluate('0.5 * (-ColL - 2 * dvW + L) * (BPult - (9 * SWt) / (10 * L * W))')
VOultW = ne.evaluate('0.5 * (BPult - (9 * SWt) / (10 * L * W)) * (-ColW - 2 * dvW + W)')

MOultL = ne.evaluate('((ColL + 2 * dvW - L) ** 2 * (-9 * SWt + 10 * BPult * L * W)) / (80 * L * W)')
MOultW = ne.evaluate('((ColW + 2 * dvW - W) ** 2 * (-9 * SWt + 10 * BPult * L * W)) / (80 * L * W)')

# Use NumPy's maximum and minimum functions
ex1L = np.minimum(ne.evaluate('(max_val / dvL + VOultL) / (2 * 200000 * AstL) * 1000', 
                              {'max_val': np.maximum(VOultL * dvL, MOultL)}), 3 / 1000)
ex1W = np.minimum(ne.evaluate('(max_val / dvW + VOultW) / (2 * 200000 * AstW) * 1000', 
                              {'max_val': np.maximum(VOultW * dvW, MOultW)}), 3 / 1000)

kvL = ne.evaluate('13 / (25 * (1 + dvL) * (1 + 1500 * ex1L))')
kvW = ne.evaluate('13 / (25 * (1 + dvW) * (1 + 1500 * ex1W))')

AngleL = ne.evaluate('29 + 7000 * ex1L')
AngleW = ne.evaluate('29 + 7000 * ex1W')

ks = np.maximum(1 / 2, ne.evaluate('(10 / 7) * (1 - D)'))

fVucL = ne.evaluate('700 * dvL * where(sqrt(fc) < 8, sqrt(fc), 8) * ks * kvL')
fVucW = ne.evaluate('700 * dvW * where(sqrt(fc) < 8, sqrt(fc), 8) * ks * kvW')

Bpr = ne.evaluate('BPmax / BP')
Mur = np.maximum(ne.evaluate('MultL / fMuoL'), ne.evaluate('MultW / fMuoW'))
VPr = ne.evaluate('VPult / fVP')
VOr = np.maximum(ne.evaluate('VOultL / fVucL'), ne.evaluate('VOultW / fVucW'))

Cost = ne.evaluate('(AstW / 1000000 * L + AstL / 1000000 * W) * 7850 * 3.400 + L * W * D * (130.866 * exp(fc * 0.0111) + 45 + 130) + 2 * D * L * W * 180')

# Create a mask for valid footings
valid_mask = ne.evaluate('(Bpr <= 1) & (Mur <= 1) & (VPr <= 1) & (VOr <= 1)')

# Flatten the arrays
L_flat = L[valid_mask]
W_flat = W[valid_mask]
D_flat = D[valid_mask]
barL_flat = barL[valid_mask]
ctsL_flat = ctsL[valid_mask]
barW_flat = barW[valid_mask]
ctsW_flat = ctsW[valid_mask]
Cost_flat = Cost[valid_mask]
fc_flat = fc[valid_mask]

# Sort by cost
sort_indices = np.argsort(Cost_flat)

# Print the results
for i in range(min(1, len(sort_indices))):
    idx = sort_indices[i]
    print(f"Footing {i+1}:")
    print(f"L={L_flat[idx]:.2f}, W={W_flat[idx]:.2f}, D={D_flat[idx]:.2f},barL={barL_flat[idx]:.3f}, ctsL={ctsL_flat[idx]:.3f}, barW={barW_flat[idx]:.3f}, ctsW={ctsW_flat[idx]:.3f},fc={fc_flat[idx]:.2f}")
    print(f"Cost: {Cost_flat[idx]:.2f}")
    print(f"SWt: {SWt[valid_mask][idx]:.2f}")
    print(f"Pult: {Pult[valid_mask][idx]:.2f}")
    print(f"BPmax: {BPmax[valid_mask][idx]:.2f}")
    print(f"BPult: {BPult[valid_mask][idx]:.2f}")
    print(f"AstL: {AstL[valid_mask][idx]:.2f}")
    print(f"AstW: {AstW[valid_mask][idx]:.2f}")
    print(f"dsL: {dsL[valid_mask][idx]:.4f}")
    print(f"dsW: {dsW[valid_mask][idx]:.4f}")
    print(f"AstminL: {AstminL[valid_mask][idx]:.2f}")
    print(f"AstminW: {AstminW[valid_mask][idx]:.2f}")
    print(f"alpha: {alpha[valid_mask][idx]:.4f}")
    print(f"gamma: {gamma[valid_mask][idx]:.4f}")
    print(f"MultL: {MultL[valid_mask][idx]:.2f}")
    print(f"MultW: {MultW[valid_mask][idx]:.2f}")
    print(f"AstshrL: {AstshrL[valid_mask][idx]:.2f}")
    print(f"AstshrW: {AstshrW[valid_mask][idx]:.2f}")
    print(f"AstreqL: {AstreqL[valid_mask][idx]:.2f}")
    print(f"AstreqW: {AstreqW[valid_mask][idx]:.2f}")
    print(f"kuL: {kuL[valid_mask][idx]:.4f}")
    print(f"kuW: {kuW[valid_mask][idx]:.4f}")
    print(f"phiL: {phiL[valid_mask][idx]:.4f}")
    print(f"phiW: {phiW[valid_mask][idx]:.4f}")
    print(f"fMuoL: {fMuoL[valid_mask][idx]:.2f}")
    print(f"fMuoW: {fMuoW[valid_mask][idx]:.2f}")
    print(f"CLR: {CLR[valid_mask][idx]:.2f}")
    print(f"VPult: {VPult[valid_mask][idx]:.2f}")
    print(f"fcv: {fcv[valid_mask][idx]:.4f}")
    print(f"fVP: {fVP[valid_mask][idx]:.2f}")
    print(f"dvL: {dvL[valid_mask][idx]:.4f}")
    print(f"dvW: {dvW[valid_mask][idx]:.4f}")
    print(f"VOultL: {VOultL[valid_mask][idx]:.2f}")
    print(f"VOultW: {VOultW[valid_mask][idx]:.2f}")
    print(f"MOultL: {MOultL[valid_mask][idx]:.2f}")
    print(f"MOultW: {MOultW[valid_mask][idx]:.2f}")
    print(f"ex1L: {ex1L[valid_mask][idx]:.6f}")
    print(f"ex1W: {ex1W[valid_mask][idx]:.6f}")
    print(f"kvL: {kvL[valid_mask][idx]:.4f}")
    print(f"kvW: {kvW[valid_mask][idx]:.4f}")
    print(f"AngleL: {AngleL[valid_mask][idx]:.2f}")
    print(f"AngleW: {AngleW[valid_mask][idx]:.2f}")
    print(f"ks: {ks[valid_mask][idx]:.4f}")
    print(f"fVucL: {fVucL[valid_mask][idx]:.2f}")
    print(f"fVucW: {fVucW[valid_mask][idx]:.2f}")
    print(f"Bpr: {Bpr[valid_mask][idx]:.2f}")
    print(f"Mur: {Mur[valid_mask][idx]:.2f}")
    print(f"VPr: {VPr[valid_mask][idx]:.2f}")
    print(f"VOr: {VOr[valid_mask][idx]:.2f}")
    print()
