import itertools
import math
from typing import Tuple, List
import time

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

def generate_valid_combinations(
    COL_LENGTH, COL_WIDTH, FTG_COVER, 
    ftg_lengths, ftg_widths, ftg_depths, ftg_conc_strengths,
    ftg_reosizes_L, ftg_reocts_L, ftg_reosizes_W, ftg_reocts_W
):
    valid_combinations = []
    count = 0
    max_count = 100000

    # Pre-compute some constant values
    col_length_less_equal = COL_LENGTH <= COL_WIDTH
    factor = 0.00129826

    for l, w in itertools.product(ftg_lengths, ftg_widths):
        if ((col_length_less_equal and l <= w) or (not col_length_less_equal and l >= w)) and not (1 < l/w <= 1.1 or 0.9 <= l/w < 1):
            for d, cs in itertools.product(ftg_depths, ftg_conc_strengths):
                d_factor = factor * d**2 / (d - FTG_COVER)
                for rl, rcl, rw, rcw in itertools.product(ftg_reosizes_L, ftg_reocts_L, ftg_reosizes_W, ftg_reocts_W):
                    if rl**2 / rcl < d_factor and rw**2 / rcw < d_factor:
                        valid_combinations.append((None, l, w, d, cs, rl, rcl, rw, rcw))
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
    load_sum = LOAD_DEAD + LOAD_LIVE
    column_diff = COLUMN_LENGTH - COLUMN_WIDTH
    common_term = math.sqrt(
        BEARINGPRESSURE * (BEARINGPRESSURE * column_diff**2 + 4 * load_sum)
    )

    # Calculate START_WIDTH and START_LENGTH
    START_WIDTH = (COLUMN_WIDTH - COLUMN_LENGTH) / 2 + common_term / (
        2 * BEARINGPRESSURE
    )
    START_LENGTH = (COLUMN_LENGTH - COLUMN_WIDTH) / 2 + common_term / (
        2 * BEARINGPRESSURE
    )

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
    print("Hash\t\t\tLength\tWidth\tDepth\tfc\tL Bar\tL CTS\tW Bar\tW CTS")
    for combination in combinations[:num_rows]:
        _, l, w, d, cs, rl, rcl, rw, rcw = combination
        hash_value = create_hash_for_combination(l, w, d, cs, rl, rcl, rw, rcw)
        print(f"{hash_value}\t{l:.3f}\t{w:.3f}\t{d:.3f}\t{cs}\t{rl:.3f}\t{rcl:.3f}\t{rw:.3f}\t{rcw:.3f}")

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
    LOAD_DEAD: float = 8000
    LOAD_LIVE: float = 2500

    # Generate arrays
    size_ranges = generate_foundation_size_ranges(
        COL_LENGTH, COL_WIDTH, BEARINGPRESSURE, LOAD_DEAD, LOAD_LIVE
    )

    # Generate valid combinations
    combinations = generate_valid_combinations(COL_LENGTH, COL_WIDTH, FTG_COVER, *size_ranges)
    
    # Filter combinations by specific criteria
    filtered_combinations = filter_combinations(combinations, length=4.6, fc=25)
    
    # Print the first 8 rows of filtered sizes
    print_foundation_sizes(filtered_combinations, 8)

    end_time = time.time()
    execution_time = end_time - start_time
    print(f"\nExecution time: {execution_time:.2f} seconds")
