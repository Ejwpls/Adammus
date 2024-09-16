import itertools
import math
from typing import Tuple, List
import time
from dataclasses import dataclass

REO_DENSITY: float = 7850  # kg/m3

def decode_hash(hash_input: List[str]) -> List[List[float]]:
    def decode_single_hash(hash_string: str) -> List[float]:
        return [
            round(int(hash_string[0:3]) * 0.05, 2),
            round(int(hash_string[3:6]) * 0.05, 2),
            round(int(hash_string[6:9]) * 0.05, 2),
            int(hash_string[9:11]),
            round((int(hash_string[11]) + 3) * 0.004, 3),
            round(int(hash_string[12:14]) * 0.005, 3),
            round((int(hash_string[14]) + 3) * 0.004, 3),
            round(int(hash_string[15:17]) * 0.005, 3)
        ]

    return [decode_single_hash(h) for h in hash_input]

def create_hash_for_combination(l, w, d, cs, rl, rcl, rw, rcw):
    return (f"{int(round(l/0.05)):03d}"
            f"{int(round(w/0.05)):03d}"
            f"{int(round(d/0.05)):03d}"
            f"{int(cs):02d}"
            f"{int(round(rl/0.004))-3:1d}"
            f"{int(round(rcl/0.005)):02d}"
            f"{int(round(rw/0.004))-3:1d}"
            f"{int(round(rcw/0.005)):02d}")

@dataclass
class FoundationData:
    l: float
    w: float
    d: float
    fc: float
    barL: float
    ctsL: float
    barW: float
    ctsW: float
    col_length: float
    col_width: float
    load_dead: float
    load_live: float
    ftg_cover: float
    bp: float

def calculate_foundation_values(data: FoundationData):
    SWt = 6 * data.d * data.l * data.w
    Pult = 1.2 * (data.load_dead + SWt) + 1.5 * data.load_live
    BPmax = (data.load_dead + SWt + data.load_live) / (data.l * data.w)
    BPult = Pult / (data.l * data.w)
    AstL = 250000 / data.ctsL * data.barL**2 * math.pi
    AstW = 250000 / data.ctsW * data.barW**2 * math.pi
    dsL = data.d - data.ftg_cover - data.barL / 2
    dsW = data.d - data.ftg_cover - data.barL - data.barW / 2
    AstminL = (228 * data.d**2 * math.sqrt(data.fc)) / dsL
    AstminW = (228 * data.d**2 * math.sqrt(data.fc)) / dsW
    alpha = 0.85 - 0.0015 * data.fc
    gamma = 0.97 - 0.0025 * data.fc
    MultL = ((7 * data.col_length - 10 * data.l) ** 2 * (-9 * SWt + 10 * BPult * data.l * data.w)) / (8000 * data.l * data.w)
    MultW = ((7 * data.col_width - 10 * data.w) ** 2 * (-9 * SWt + 10 * BPult * data.l * data.w)) / (8000 * data.l * data.w)
    AstshrL = 5 * MultL / (2 * dsL)
    AstshrW = 5 * MultW / (2 * dsW)
    AstreqL = max(AstminL, AstshrL)
    AstreqW = max(AstminW, AstshrW)
    
    # Add this check
    if AstreqL > AstL or AstreqW > AstW:
        return None  # Invalid combination

    kuL = AstreqL / (2000 * alpha * dsL * data.fc * gamma * data.l)
    kuW = AstreqW / (2000 * alpha * dsW * data.fc * gamma * data.w)
    phiL = min(0.85, max(0.65, 1.24 - 13 * kuL / 12))
    phiW = min(0.85, max(0.65, 1.24 - 13 * kuW / 12))
    fMuoL = (AstL * dsL * phiL) / 2 - (AstL**2 * phiL) / (8000 * alpha * data.fc)
    fMuoW = (AstW * dsW * phiW) / 2 - (AstW**2 * phiW) / (8000 * alpha * data.fc)
    CLR = BPult * (data.col_width + dsL) * (data.col_length + dsL)
    VPult = Pult - CLR
    fcv = 0.17 * (1 + 2 / max(data.col_length / data.col_width, data.col_width / data.col_length, 2)) * math.sqrt(data.fc)
    fVP = 1400 * dsW * (data.col_length + data.col_width + 2 * dsW) * fcv
    dvL = max(0.9 * dsL, 0.72 * data.d)
    dvW = max(0.9 * dsW, 0.72 * data.d)
    VOultL = 0.5 * (-data.col_length - 2 * dvW + data.l) * (BPult - (9 * SWt) / (10 * data.l * data.w))
    VOultW = 0.5 * (BPult - (9 * SWt) / (10 * data.l * data.w)) * (-data.col_width - 2 * dvW + data.w)
    MOultL = ((data.col_length + 2 * dvW - data.l) ** 2 * (-9 * SWt + 10 * BPult * data.l * data.w)) / (80 * data.l * data.w)
    MOultW = ((data.col_width + 2 * dvW - data.w) ** 2 * (-9 * SWt + 10 * BPult * data.l * data.w)) / (80 * data.l * data.w)
    ex1L = min((max(VOultL * dvL, MOultL) / dvL + VOultL) / (2 * 200000 * AstL) * 1000, 3 / 1000)
    ex1W = min((max(VOultW * dvW, MOultW) / dvW + VOultW) / (2 * 200000 * AstW) * 1000, 3 / 1000)
    kvL = 13 / (25 * (1 + dvL) * (1 + 1500 * ex1L))
    kvW = 13 / (25 * (1 + dvW) * (1 + 1500 * ex1W))
    AngleL = 29 + 7000 * ex1L
    AngleW = 29 + 7000 * ex1W
    ks = max(1 / 2, (10 / 7) * (1 - data.d))
    fVucL = 700 * dvL * math.sqrt(data.fc) * ks * kvL
    fVucW = 700 * dvW * math.sqrt(data.fc) * ks * kvW
    Bpr = BPmax / data.bp
    Mur = max(MultL / fMuoL, MultW / fMuoW)
    VPr = VPult / fVP
    VOr = max(VOultL / fVucL, VOultW / fVucW)

    cost = (AstW / 1000000 * data.l + AstL / 1000000 * data.w) * 7850 * 3.400 + data.l * data.w * data.d * (130.866 * math.exp(data.fc * 0.0111) + 45 + 130) + 2 * data.d * data.l * data.w * 180

    return cost, Bpr, Mur, VPr, VOr, AstreqL, AstreqW, AstL, AstW

def generate_valid_combinations(
    COL_LENGTH, COL_WIDTH, FTG_COVER, 
    ftg_lengths, ftg_widths, ftg_depths, ftg_conc_strengths,
    ftg_reosizes_L, ftg_reocts_L, ftg_reosizes_W, ftg_reocts_W,
    LOAD_DEAD, LOAD_LIVE, BEARINGPRESSURE
):
    valid_combinations = []
    count = 0
    max_count = 10

    # Pre-compute some constant values
    factor = 0.00129826

    for l, w in itertools.product(ftg_lengths, ftg_widths):
        if ((COL_LENGTH == COL_WIDTH and l == w) or (COL_LENGTH > COL_WIDTH and l > w) or (COL_WIDTH > COL_LENGTH and w > l)) and abs(l - w) <= 0.05 + abs(COL_LENGTH - COL_WIDTH):
            print(l, w)
            for d, cs in itertools.product(ftg_depths, ftg_conc_strengths):
                d_factor = factor * d**2 / (d - FTG_COVER)
                reo_combinations = list(itertools.product(ftg_reosizes_L, ftg_reocts_L, ftg_reosizes_W, ftg_reocts_W))
                reo_combinations.sort(key=lambda x: x[0]**2 / x[1] + x[2]**2 / x[3])
                for rl, rcl, rw, rcw in reo_combinations:
                    if rl**2 / rcl < d_factor and rw**2 / rcw < d_factor:
                        foundation_data = FoundationData(l, w, d, cs, rl, rcl, rw, rcw, COL_LENGTH, COL_WIDTH, LOAD_DEAD, LOAD_LIVE, FTG_COVER, BEARINGPRESSURE)
                        result = calculate_foundation_values(foundation_data)
                        if result is not None:
                            cost, Bpr, Mur, VPr, VOr, AstreqL, AstreqW, AstL, AstW = result
                            if Bpr <= 1 and Mur <= 1 and VPr <= 1 and VOr <= 1 and AstreqL <= AstL and AstreqW <= AstW:
                                valid_combinations.append((cost, l, w, d, cs, rl, rcl, rw, rcw, Bpr, Mur, VPr, VOr))
                                count += 1
                                if count == max_count:
                                    return valid_combinations

    return valid_combinations

def generate_foundation_size_ranges(
    COLUMN_LENGTH: float,
    COLUMN_WIDTH: float,
    BEARINGPRESSURE: float,
    LOAD_DEAD: float,
    LOAD_LIVE: float,
) -> Tuple[
    List[float],
    List[float],
    List[float],
    List[float],
    List[float],
    List[float],
    List[float],
    List[float],
]:
    # Calculate common terms once
    area = (LOAD_DEAD + LOAD_LIVE)/BEARINGPRESSURE

    # Calculate START_WIDTH and START_LENGTH
    START_WIDTH = 0.5*(COLUMN_WIDTH - COLUMN_LENGTH + math.sqrt(4*area+(COLUMN_LENGTH-COLUMN_WIDTH)**2))
    START_LENGTH = COLUMN_LENGTH - COLUMN_WIDTH + START_WIDTH

    # Round up to nearest 0.05 and calculate max values
    FTG_LENGTH_MIN = math.ceil(START_LENGTH / 0.05) * 0.05
    FTG_WIDTH_MIN = math.ceil(START_WIDTH / 0.05) * 0.05
    FTG_LENGTH_MAX = FTG_LENGTH_MIN + 0.5
    FTG_WIDTH_MAX = FTG_WIDTH_MIN + 0.5

    # Generate lists more efficiently
    FTG_LENGTHS = [round(x, 2) for x in range(int(FTG_LENGTH_MIN * 100), int(FTG_LENGTH_MAX * 100) + 6, 5)]
    FTG_LENGTHS = [x / 100 for x in FTG_LENGTHS]

    FTG_WIDTHS = [round(x, 2) for x in range(int(FTG_WIDTH_MIN * 100), int(FTG_WIDTH_MAX * 100) + 6, 5)]
    FTG_WIDTHS = [x / 100 for x in FTG_WIDTHS]

    # Calculate FTG_DPTH_MIN more efficiently
    column_sum = COLUMN_LENGTH + COLUMN_WIDTH
    FTG_DPTH_MIN = (
        math.ceil(
            (
                -column_sum / 4
                + 0.05
                * math.sqrt(
                    25 * column_sum**2
                    + 3 / 238 * math.sqrt(5) * (4 * LOAD_DEAD + 5 * LOAD_LIVE)
                )
            )
            / 0.05
        )
        * 0.05
    )

    FTG_DPTH_MAX = FTG_DPTH_MIN + 0.5
    FTG_DPTHS = [round(x, 2) for x in range(int(FTG_DPTH_MIN * 100), int(FTG_DPTH_MAX * 100) + 6, 5)]
    FTG_DPTHS = [x / 100 for x in FTG_DPTHS]

    # Use predefined lists for constant values
    FTG_CONC_STRENGTHS = [25, 32, 40, 50, 65]
    FTG_REO_SIZES = [0.012, 0.016, 0.02, 0.024, 0.028, 0.032, 0.036, 0.04]
    FTG_REO_CTS = sorted(set([round(x, 3) for x in range(100, 301, 25)] + [round(x, 3) for x in range(100, 301, 20)]))
    FTG_REO_CTS = [x / 1000 for x in FTG_REO_CTS]
    

    return (
        FTG_LENGTHS,
        FTG_WIDTHS,
        FTG_DPTHS,
        FTG_CONC_STRENGTHS,
        FTG_REO_SIZES,
        FTG_REO_CTS,
        FTG_REO_SIZES,
        FTG_REO_CTS,
    )

def print_foundation_sizes(combinations, num_rows):
    print("\nFoundation Sizes:")
    print("Length\tWidth\tDepth\tfc\tL Bar\tL CTS\tW Bar\tW CTS\tBpr\tMur\tVPr\tVOr\tCost")
    for combination in combinations[:num_rows]:
        cost, l, w, d, cs, rl, rcl, rw, rcw, Bpr, Mur, VPr, VOr = combination
        print(f"{l:.3f}\t{w:.3f}\t{d:.3f}\t{cs}\t{rl:.3f}\t{rcl:.3f}\t{rw:.3f}\t{rcw:.3f}\t{Bpr:.3f}\t{Mur:.3f}\t{VPr:.3f}\t{VOr:.3f}\t{cost:.2f}")

def filter_combinations(combinations, **kwargs):
    def matches_criteria(combo, criteria):
        return all(
            getattr(combo, attr, combo[index]) == value
            for (attr, index), value in criteria.items()
        )

    criteria = {
        ('length', 1): kwargs.get('length'),
        ('width', 2): kwargs.get('width'),
        ('depth', 3): kwargs.get('depth'),
        ('fc', 4): kwargs.get('fc'),
        ('reo_size_L', 5): kwargs.get('reo_size_L'),
        ('reo_cts_L', 6): kwargs.get('reo_cts_L'),
        ('reo_size_W', 7): kwargs.get('reo_size_W'),
        ('reo_cts_W', 8): kwargs.get('reo_cts_W')
    }
    
    # Remove None values from criteria
    criteria = {k: v for k, v in criteria.items() if v is not None}

    return [combo for combo in combinations if matches_criteria(combo, criteria)]

# Example usage
if __name__ == "__main__":
    start_time = time.time()

    COL_LENGTH: float = 0.5
    COL_WIDTH: float = 0.5
    FTG_COVER: float = 0.060
    BEARINGPRESSURE: float = 500
    LOAD_DEAD: float = 4000
    LOAD_LIVE: float = 1500

    # Generate arrays
    size_ranges = generate_foundation_size_ranges(
        COL_LENGTH, COL_WIDTH, BEARINGPRESSURE, LOAD_DEAD, LOAD_LIVE
    )

    # Generate valid combinations
    combinations = generate_valid_combinations(COL_LENGTH, COL_WIDTH, FTG_COVER, *size_ranges, LOAD_DEAD, LOAD_LIVE, BEARINGPRESSURE)
    
    # Filter combinations by specific criteria
    # filtered_combinations = filter_combinations(combinations, length=4.6, fc=25)
    
    # Print the first 8 rows of filtered sizes
    print_foundation_sizes(combinations,10)

    end_time = time.time()
    execution_time = end_time - start_time
    print(f"\nExecution time: {execution_time:.2f} seconds")
