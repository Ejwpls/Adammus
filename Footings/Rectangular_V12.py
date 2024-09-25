import numpy as np
import numexpr as ne
import math
import dataclasses as dc
from foundation_export import create_mem_file

GC_PI = math.pi

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

def create_concrete_combos(start_L, end_L, start_D, end_D, step):
    L = np.arange(start_L, end_L + 0.05, 0.05)
    D = np.arange(start_D, end_D + 0.05, 0.05)
    fc_values = np.array([25, 32, 40, 50, 65])
    
    L_mesh, D_mesh, fc_mesh = np.meshgrid(L, D, fc_values, indexing='ij')
    W_mesh = L_mesh + step
    
    combinations = np.column_stack((L_mesh.ravel(), W_mesh.ravel(), D_mesh.ravel(), fc_mesh.ravel()))
    
    product_ABCD = ne.evaluate('L_vals * W_vals * D_vals * fc_vals', {
        'L_vals': combinations[:, 0],
        'W_vals': combinations[:, 1],
        'D_vals': combinations[:, 2],
        'fc_vals': combinations[:, 3]
    })
    
    sorted_indices = np.argsort(product_ABCD)
    sorted_combinations = combinations[sorted_indices]
    
    # Filter the combinations
    A, B, C = sorted_combinations[:, 0], sorted_combinations[:, 1], sorted_combinations[:, 2]
    mask = ne.evaluate('(0 < (8000 + 1500 + 6 * A * B * C) / (250 * A * B)) & ((8000 + 1500 + 6 * A * B * C) / (250 * A * B) < 1)')
    
    return sorted_combinations[mask]

def create_steel_combos(conc_arr):
    A = np.array([0.012,0.016,0.020,0.024,0.028,0.032,0.036,0.040])
    B = np.array([0.1, 0.120, 0.125, 0.14, 0.15, 0.16, 0.175, 0.18, 0.2, 0.22, 0.225, 0.24, 0.25, 0.26, 0.275, 0.28, 0.3])
    
    A_mesh1, B_mesh1, A_mesh2, B_mesh2 = np.meshgrid(A, B, A, B, indexing='ij')
    combinations = np.stack([A_mesh1.ravel(), B_mesh1.ravel(), A_mesh2.ravel(), B_mesh2.ravel()], axis=-1)
    
    product_ABAB = ne.evaluate('A_vals1 * B_vals1 * A_vals2 * B_vals2', {
        'A_vals1': combinations[:, 0],
        'B_vals1': combinations[:, 1],
        'A_vals2': combinations[:, 2],
        'B_vals2': combinations[:, 3]
    })
    
    sorted_indices = np.argsort(product_ABAB)
    sorted_combinations = combinations[sorted_indices]
    
    # Filter the combinations
    depth = np.min(conc_arr[:, -2])
    
    context = {
        'barL': sorted_combinations[:, 0],
        'ctsL': sorted_combinations[:, 1],
        'barW': sorted_combinations[:, 2],
        'ctsW': sorted_combinations[:, 3],
        'D': depth,
        'pi': math.pi,
        'sqrt_25': math.sqrt(25)
    }

    condition1 = ne.evaluate('31250 * barL**2 * (-barL - 2*0.05 + 2*D) * pi / (57 * ctsL * D**2 * sqrt_25) > 1', context)
    condition2 = ne.evaluate('31250 * barW**2 * (-barW - 2*0.05 + 2*D) * pi / (57 * ctsW * D**2 * sqrt_25) > 1', context)

    mask = np.logical_and(condition1, condition2)
    
    return sorted_combinations[mask]

def calculate_foundation_values(fdn):
    vars = {
        'L': fdn.L, 'W': fdn.W, 'D': fdn.D, 'Cvr': fdn.Cvr,
        'ColL': fdn.ColL, 'ColW': fdn.ColW, 'Pdl': fdn.Pdl,
        'Pll': fdn.Pll, 'BP': fdn.BP, 'fc': fdn.fc,
        'barL': fdn.barL, 'ctsL': fdn.ctsL,
        'barW': fdn.barW, 'ctsW': fdn.ctsW
    }

    vars['SWt'] = ne.evaluate('6 * D * L * W', vars)
    vars['Pult'] = ne.evaluate('1.2 * (Pdl + SWt) + 1.5 * Pll', vars)
    vars['BPmax'] = ne.evaluate('(Pdl + SWt + Pll) / (L * W)', vars)
    vars['BPult'] = ne.evaluate('Pult / (L * W)', vars)
    vars['AstL'] = ne.evaluate('250000 / ctsL * barL**2 * GC_PI', vars)
    vars['AstW'] = ne.evaluate('250000 / ctsW * barW**2 * GC_PI', vars)
    vars['dsL'] = ne.evaluate('D - Cvr - barL / 2', vars)
    vars['dsW'] = ne.evaluate('D - Cvr - barL - barW / 2', vars)
    vars['AstminL'] = ne.evaluate('(228 * D**2 * sqrt(fc)) / dsL', vars)
    vars['AstminW'] = ne.evaluate('(228 * D**2 * sqrt(fc)) / dsW', vars)
    vars['alpha'] = ne.evaluate('0.85 - 0.0015 * fc', vars)
    vars['gamma'] = ne.evaluate('0.97 - 0.0025 * fc', vars)
    vars['MultL'] = ne.evaluate('((7 * ColL - 10 * L) ** 2 * (-9 * SWt + 10 * BPult * L * W)) / (8000 * L * W)', vars)
    vars['MultW'] = ne.evaluate('((7 * ColW - 10 * W) ** 2 * (-9 * SWt + 10 * BPult * L * W)) / (8000 * L * W)', vars)
    vars['AstshrL'] = ne.evaluate('100 * MultL / (dsL * (50 - 9 * gamma))', vars)
    vars['AstshrW'] = ne.evaluate('100 * MultW / (dsW * (50 - 9 * gamma))', vars)
    vars['AstreqL'] = ne.evaluate('where(AstminL > AstshrL, AstminL, AstshrL)', vars)
    vars['AstreqW'] = ne.evaluate('where(AstminW > AstshrW, AstminW, AstshrW)', vars)
    vars['kuL'] = ne.evaluate('AstreqL / (2000 * alpha * dsL * fc * gamma * L)', vars)
    vars['kuW'] = ne.evaluate('AstreqW / (2000 * alpha * dsW * fc * gamma * W)', vars)
    vars['phiL'] = ne.evaluate('where(0.85 < 1.24 - 13 * kuL / 12, 0.85, where(0.65 > 1.24 - 13 * kuL / 12, 0.65, 1.24 - 13 * kuL / 12))', vars)
    vars['phiW'] = ne.evaluate('where(0.85 < 1.24 - 13 * kuW / 12, 0.85, where(0.65 > 1.24 - 13 * kuW / 12, 0.65, 1.24 - 13 * kuW / 12))', vars)
    vars['fMuoL'] = ne.evaluate('(AstL * dsL * phiL) / 2 - (AstL**2 * phiL) / (8000 * alpha * fc)', vars)
    vars['fMuoW'] = ne.evaluate('(AstW * dsW * phiW) / 2 - (AstW**2 * phiW) / (8000 * alpha * fc)', vars)
    vars['CLR'] = ne.evaluate('BPult * (ColW + dsL) * (ColL + dsL)', vars)
    vars['VPult'] = ne.evaluate('Pult - CLR', vars)
    vars['fcv'] = ne.evaluate('0.17 * (1 + 2 / where(ColL / ColW > ColW / ColL, where(ColL / ColW > 2.0, ColL / ColW, 2.0), where(ColW / ColL > 2.0, ColW / ColL, 2.0))) * sqrt(fc)', vars)
    vars['fVP'] = ne.evaluate('1400 * dsW * (ColL + ColW + 2 * dsW) * fcv', vars)
    vars['dvL'] = ne.evaluate('where(0.9 * dsL > 0.72 * D, 0.9 * dsL, 0.72 * D)', vars)
    vars['dvW'] = ne.evaluate('where(0.9 * dsW > 0.72 * D, 0.9 * dsW, 0.72 * D)', vars)
    vars['VOultL'] = ne.evaluate('0.5 * (-ColL - 2 * dvW + L) * (BPult - (9 * SWt) / (10 * L * W))', vars)
    vars['VOultW'] = ne.evaluate('0.5 * (BPult - (9 * SWt) / (10 * L * W)) * (-ColW - 2 * dvW + W)', vars)
    vars['MOultL'] = ne.evaluate('((ColL + 2 * dvW - L) ** 2 * (-9 * SWt + 10 * BPult * L * W)) / (80 * L * W)', vars)
    vars['MOultW'] = ne.evaluate('((ColW + 2 * dvW - W) ** 2 * (-9 * SWt + 10 * BPult * L * W)) / (80 * L * W)', vars)
    vars['ex1L'] = ne.evaluate('where((where(VOultL * dvL > MOultL, VOultL * dvL, MOultL) / dvL + VOultL) / (2 * 200000 * AstL) * 1000 < 3 / 1000, (where(VOultL * dvL > MOultL, VOultL * dvL, MOultL) / dvL + VOultL) / (2 * 200000 * AstL) * 1000, 3 / 1000)', vars)
    vars['ex1W'] = ne.evaluate('where((where(VOultW * dvW > MOultW, VOultW * dvW, MOultW) / dvW + VOultW) / (2 * 200000 * AstW) * 1000 < 3 / 1000, (where(VOultW * dvW > MOultW, VOultW * dvW, MOultW) / dvW + VOultW) / (2 * 200000 * AstW) * 1000, 3 / 1000)', vars)
    vars['kvL'] = ne.evaluate('13 / (25 * (1 + dvL) * (1 + 1500 * ex1L))', vars)
    vars['kvW'] = ne.evaluate('13 / (25 * (1 + dvW) * (1 + 1500 * ex1W))', vars)
    vars['AngleL'] = ne.evaluate('29 + 7000 * ex1L', vars)
    vars['AngleW'] = ne.evaluate('29 + 7000 * ex1W', vars)
    vars['ks'] = ne.evaluate('where(1 / 2 > (10 / 7) * (1 - D), 1 / 2, (10 / 7) * (1 - D))', vars)
    vars['fVucL'] = ne.evaluate('700 * dvL * where(sqrt(fc) < 8, sqrt(fc), 8) * ks * kvL', vars)
    vars['fVucW'] = ne.evaluate('700 * dvW * where(sqrt(fc) < 8, sqrt(fc), 8) * ks * kvW', vars)
    vars['Bpr'] = ne.evaluate('BPmax / BP', vars)
    vars['Mur'] = ne.evaluate('where(MultL / fMuoL > MultW / fMuoW, MultL / fMuoL, MultW / fMuoW)', vars)
    vars['VPr'] = ne.evaluate('VPult / fVP', vars)
    vars['VOr'] = ne.evaluate('where(VOultL / fVucL > VOultW / fVucW, VOultL / fVucL, VOultW / fVucW)', vars)
    vars['Cost'] = ne.evaluate('(AstW / 1000000 * (L + barL*0.01165+0.0202) + AstL / 1000000 * (W + barW*0.01165+0.0202)) * 7850 * 3.400 + L * W * D * (130.866 * exp(fc * 0.0111) + 45 + 130) + 2 * D * L * W * 180', vars)

    return vars

def print_foundation_results(results, idx):
    print("-----------------------------------------")
    
    def print_row(label, value_l, value_w=None, unit='', decimals=2):
        format_str = f"{{:>{decimals}.{decimals}f}}"
        if value_w is None:
            print(f"{label}\t\t\t{format_str.format(value_l)} {unit}")
        else:
            print(f"{label}\t\tL: {format_str.format(value_l)} {unit}\tW: {format_str.format(value_w)} {unit}")

    print_row("Pdl", results['Pdl'][idx], unit="kN", decimals=0)
    print_row("Pll", results['Pll'][idx], unit="kN", decimals=0)
    print_row("Col L", results['ColL'][idx]*1000, unit="mm", decimals=0)
    print_row("Col W", results['ColW'][idx]*1000, unit="mm", decimals=0)
    print_row("BP", results['BP'][idx], unit="kPa", decimals=0)
    print_row("Cvr", results['Cvr'][idx]*1000, unit="mm", decimals=0)
    print("-----------------------------------------")
    print_row("L", results['L'][idx]*1000, unit="mm", decimals=0)
    print_row("W", results['W'][idx]*1000, unit="mm", decimals=0)
    print_row("D", results['D'][idx]*1000, unit="mm", decimals=0)
    print_row("fc", results['fc'][idx], unit="MPa", decimals=0)
    print_row("bar", results['barL'][idx]*1000, results['barW'][idx]*1000, "mm", 0)
    print_row("cts", results['ctsL'][idx]*1000, results['ctsW'][idx]*1000, "mm", 0)
    print("-----------------------------------------")
    print_row("SWt", results['SWt'][idx], unit="kN", decimals=0)
    print_row("Pult", results['Pult'][idx], unit="kN", decimals=0)
    print_row("BPmax", results['BPmax'][idx], unit="kPa", decimals=1)
    print_row("BPult", results['BPult'][idx], unit="kPa", decimals=1)
    print_row("Ast", results['AstL'][idx], results['AstW'][idx], "mm²", 0)
    print_row("ds", results['dsL'][idx]*1000, results['dsW'][idx]*1000, "mm", 0)
    print_row("Astmin", results['AstminL'][idx], results['AstminW'][idx], "mm²", 0)
    print_row("alpha", results['alpha'][idx], decimals=4)
    print_row("gamma", results['gamma'][idx], decimals=4)
    print_row("Mult", results['MultL'][idx], results['MultW'][idx], "kNm", 1)
    print_row("Astshr", results['AstshrL'][idx], results['AstshrW'][idx], "mm²", 0)
    print_row("Astreq", results['AstreqL'][idx], results['AstreqW'][idx], "mm²", 0)
    print_row("ku", results['kuL'][idx], results['kuW'][idx], decimals=3)
    print_row("phi", results['phiL'][idx], results['phiW'][idx], decimals=3)
    print_row("fMuo", results['fMuoL'][idx], results['fMuoW'][idx], "kNm", 1)
    print_row("CLR", results['CLR'][idx], unit="kN", decimals=0)
    print_row("VPult", results['VPult'][idx], unit="kN", decimals=1)
    print_row("fcv", results['fcv'][idx], unit="kN", decimals=1)
    print_row("fVP", results['fVP'][idx], unit="kN", decimals=1)
    print_row("dv", results['dvL'][idx]*1000, results['dvW'][idx]*1000, "mm", 0)
    print_row("VOult", results['VOultL'][idx], results['VOultW'][idx], "kN", 1)
    print_row("MOult", results['MOultL'][idx], results['MOultW'][idx], "kNm", 1)
    print_row("ex1",results['ex1L'][idx], results['ex1W'][idx], "", 6)
    print_row("kv", results['kvL'][idx], results['kvW'][idx], decimals=3)
    print_row("Angle", results['AngleL'][idx], results['AngleW'][idx], "°", 1)
    print_row("ks", results['ks'][idx], decimals=3)
    print_row("fVuc", results['fVucL'][idx], results['fVucW'][idx], "kN", 1)
    print("-----------------------------------------")
    print_row("Bpr", results['Bpr'][idx], decimals=3)
    print_row("Mur", results['Mur'][idx], decimals=3)
    print_row("VPr", results['VPr'][idx], decimals=3)
    print_row("VOr", results['VOr'][idx], decimals=3)
    print_row("Cost", results['Cost'][idx], unit="$", decimals=2)
    print("-----------------------------------------")

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
    
    # Calculate all values
    results = calculate_foundation_values(startingfdn)
    
    # Use the results directly
    mask = results['fVucL'] / results['VOultL'] > 1
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
    
    # Calculate all values
    results = calculate_foundation_values(startingfdn)
    
    # Use the results directly
    mask = results['fVucL'] / results['VOultL'] > 1
    return np.min(startingfdn.D[mask])

def calc_max_L(fdn):
    _diff = fdn.ColW-fdn.ColL
    _temp_a = fdn.BP-6*calc_max_D(fdn)
    max_L = (_diff*_temp_a-np.sqrt(_temp_a)*np.sqrt(4*(fdn.Pdl+fdn.Pll)+_diff**2*_temp_a))/(-2*_temp_a)
    return (np.ceil(max_L/0.05)*0.05)[0]

def calc_max_W(fdn):
    return (calc_max_L(fdn)-fdn.ColL+fdn.ColW)[0]

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

def main():
    test = Foundation(L=6.3, W=6.0, D=2.55, ColL=0.8, ColW=0.5, Pdl=8000, Pll=1500, BP=250, fc=25, barL=0.012, ctsL=0.3, barW=0.012, ctsW=0.3, Cvr=0.05)

    conc = create_concrete_combos(calc_min_L(test), calc_max_L(test), calc_min_D(test), calc_max_D(test), test.ColW - test.ColL)
    reo = create_steel_combos(conc)

    combined_fdn = Foundation(
        L=conc[:, np.newaxis, 0].repeat(reo.shape[0], axis=1),
        W=conc[:, np.newaxis, 1].repeat(reo.shape[0], axis=1),
        D=conc[:, np.newaxis, 2].repeat(reo.shape[0], axis=1),
        fc=conc[:, np.newaxis, 3].repeat(reo.shape[0], axis=1),
        ColL=np.full_like(conc[:, np.newaxis, 0].repeat(reo.shape[0], axis=1), test.ColL),
        ColW=np.full_like(conc[:, np.newaxis, 0].repeat(reo.shape[0], axis=1), test.ColW),
        Pdl=np.full_like(conc[:, np.newaxis, 0].repeat(reo.shape[0], axis=1), test.Pdl),
        Pll=np.full_like(conc[:, np.newaxis, 0].repeat(reo.shape[0], axis=1), test.Pll),
        BP=np.full_like(conc[:, np.newaxis, 0].repeat(reo.shape[0], axis=1), test.BP),
        barL=reo[np.newaxis, :, 0].repeat(conc.shape[0], axis=0),
        ctsL=reo[np.newaxis, :, 1].repeat(conc.shape[0], axis=0),
        barW=reo[np.newaxis, :, 2].repeat(conc.shape[0], axis=0),
        ctsW=reo[np.newaxis, :, 3].repeat(conc.shape[0], axis=0),
        Cvr=np.full_like(conc[:, np.newaxis, 0].repeat(reo.shape[0], axis=1), test.Cvr)
    )

    results = calculate_foundation_values(combined_fdn)

    valid_mask = ne.evaluate('(Bpr <= 1) & (Mur <= 1) & (VPr <= 1) & (VOr <= 1)', results)

    flat_results = {k: v.flatten()[valid_mask.flatten()] for k, v in results.items() if isinstance(v, np.ndarray)}
    sort_indices = np.argsort(flat_results['Cost'])


    # After finding the optimal foundation
    foundation_export_parameters = {
        'foundation_ID': f'PF_TEST',
        'author': 'Bobda Builder',
        'date': '[date]', 
        'fc': flat_results['fc'][sort_indices[0]],
        'L': flat_results['L'][sort_indices[0]],
        'W': flat_results['W'][sort_indices[0]],
        'D': flat_results['D'][sort_indices[0]],
        'ColL': flat_results['ColL'][sort_indices[0]],
        'ColW': flat_results['ColW'][sort_indices[0]],
        'Pdl': flat_results['Pdl'][sort_indices[0]],
        'Pll': flat_results['Pll'][sort_indices[0]],
        'BP': flat_results['BP'][sort_indices[0]],
        'Cvr': flat_results['Cvr'][sort_indices[0]],
        'barL': flat_results['barL'][sort_indices[0]],
        'CTSL': flat_results['ctsL'][sort_indices[0]],
        'barW': flat_results['barW'][sort_indices[0]],
        'CTSW': flat_results['ctsW'][sort_indices[0]]
    }

    create_mem_file(foundation_export_parameters)

    num_results_to_print = min(1, len(sort_indices))
    for i in range(num_results_to_print):
        print_foundation_results(flat_results, sort_indices[i])

if __name__ == "__main__":
    main()