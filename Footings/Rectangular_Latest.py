import numpy as np
from numba import jit, vectorize
import math

REO_DENSITY: float = 7850  # kg/m3

@vectorize(['boolean(float64, float64, float64, float64)'])
def check_reinforcement_constraint(bar_size, spacing, depth, FTG_cover):
    return 250000 * bar_size**2 * np.pi / spacing < 1019.65 * depth**2 / (depth - FTG_cover)

@jit(nopython=True)
def is_valid_lw(l, w, col_length, col_width):
    return ((col_length <= col_width and l <= w) or (col_length > col_width and l >= w)) and not (1 < l/w <= 1.1 or 0.9 <= l/w < 1)

@jit(nopython=True)
def generate_valid_combinations(col_length, col_width, FTG_cover, ftg_lengths, ftg_widths, ftg_depths, 
                                ftg_conc_strengths, ftg_reosizes_L, ftg_reocts_L, ftg_reosizes_W, ftg_reocts_W):
    valid_combinations = []
    for l in ftg_lengths:
        for w in ftg_widths:
            if is_valid_lw(l, w, col_length, col_width):
                for d in ftg_depths:
                    for cs in ftg_conc_strengths:
                        for rl in ftg_reosizes_L:
                            for rcl in ftg_reocts_L:
                                for rw in ftg_reosizes_W:
                                    for rcw in ftg_reocts_W:
                                        if (check_reinforcement_constraint(rl, rcl, d, FTG_cover) and
                                            check_reinforcement_constraint(rw, rcw, d, FTG_cover)):
                                            valid_combinations.append((l, w, d, cs, rl, rcl, rw, rcw))
    return valid_combinations

def generate_foundation_size_ranges(
    COLUMN_LENGTH: float,
    COLUMN_WIDTH: float,
    BEARINGPRESSURE: float,
    LOAD_DEAD: float,
    LOAD_LIVE: float
):
    # Calculate common terms once
    load_sum = LOAD_DEAD + LOAD_LIVE
    column_diff = COLUMN_LENGTH - COLUMN_WIDTH
    common_term = math.sqrt(BEARINGPRESSURE * (BEARINGPRESSURE * column_diff**2 + 4 * load_sum))

    # Calculate START_WIDTH and START_LENGTH
    START_WIDTH = (COLUMN_WIDTH - COLUMN_LENGTH) / 2 + common_term / (2 * BEARINGPRESSURE)
    START_LENGTH = (COLUMN_LENGTH - COLUMN_WIDTH) / 2 + common_term / (2 * BEARINGPRESSURE)

    # Round up to nearest 0.05 and calculate max values
    FTG_LENGTH_MIN = math.ceil(START_LENGTH / 0.05) * 0.05
    FTG_WIDTH_MIN = math.ceil(START_WIDTH / 0.05) * 0.05
    FTG_LENGTH_MAX = FTG_LENGTH_MIN + 0.5
    FTG_WIDTH_MAX = FTG_WIDTH_MIN + 0.5

    # Generate arrays more efficiently
    FTG_LENGTHS = np.arange(FTG_LENGTH_MIN, FTG_LENGTH_MAX + 0.051, 0.05).round(2)
    FTG_WIDTHS = np.arange(FTG_WIDTH_MIN, FTG_WIDTH_MAX + 0.051, 0.05).round(2)

    # Calculate FTG_DPTH_MIN more efficiently
    column_sum = COLUMN_LENGTH + COLUMN_WIDTH
    FTG_DPTH_MIN = math.ceil(
        (
            -column_sum / 4
            + 0.05 * math.sqrt(
                25 * column_sum**2
                + 3 / 238 * math.sqrt(5) * (4 * LOAD_DEAD + 5 * LOAD_LIVE)
            )
        )
        / 0.05
    ) * 0.05

    FTG_DPTH_MAX = FTG_DPTH_MIN + 0.5
    FTG_DPTHS = np.arange(FTG_DPTH_MIN, FTG_DPTH_MAX + 0.051, 0.05).round(2)

    # Use predefined arrays for constant values
    FTG_CONC_STRENGTHS = np.array([25, 32, 40, 50, 65], dtype=np.float32)
    FTG_REO_SIZES = np.array([0.012, 0.016, 0.02, 0.024, 0.028, 0.032, 0.036, 0.04], dtype=np.float32)
    FTG_REO_CTS = np.unique(np.concatenate([
        np.arange(0.1, 0.301, 0.025),
        np.arange(0.1, 0.301, 0.02)
    ])).round(3)

    return (FTG_LENGTHS, FTG_WIDTHS, FTG_DPTHS, FTG_CONC_STRENGTHS, FTG_REO_SIZES, FTG_REO_CTS, FTG_REO_SIZES, FTG_REO_CTS)

# Example usage
if __name__ == "__main__":
    COL_LENGTH = 0.5
    COL_WIDTH = 0.5
    FTG_COVER = 0.060
    BEARINGPRESSURE = 500  # Example value, adjust as needed
    LOAD_DEAD = 8000  # Example value, adjust as needed
    LOAD_LIVE = 2500  # Example value, adjust as needed
    
    # Generate arrays and pass them directly to generate_valid_combinations
    results = generate_valid_combinations(
        COL_LENGTH, COL_WIDTH, FTG_COVER, 
        *generate_foundation_size_ranges(COL_LENGTH, COL_WIDTH, BEARINGPRESSURE, LOAD_DEAD, LOAD_LIVE)
    )

    print(f"Number of valid combinations: {len(results)}")