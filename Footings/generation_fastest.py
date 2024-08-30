# import itertools
# import numpy as np

# COL_LENGTH = 0.5
# COL_WIDTH = 0.4

# # Define your arrays
# ftg_lengths = np.array([4.0, 4.05, 4.1, 4.15, 4.2, 4.25, 4.3, 4.35, 4.4, 4.45, 4.5])
# ftg_widths = np.array([4.0, 4.05, 4.1, 4.15, 4.2, 4.25, 4.3, 4.35, 4.4, 4.45, 4.5])
# ftg_depths = np.array([1.5, 1.55, 1.6, 1.65, 1.7, 1.75, 1.8, 1.85, 1.9, 1.95, 2.0, 2.05, 2.1, 2.15, 2.2, 2.25, 2.3, 2.35, 2.4, 2.45, 2.5])
# ftg_conc_strengths = np.array([25.0, 32.0, 40.0, 50.0, 65.0])
# ftg_reosizes_L = np.array([0.012, 0.016, 0.02, 0.024, 0.028, 0.032, 0.036, 0.04])
# ftg_reocts_L = np.array([0.1, 0.12, 0.125, 0.14, 0.15, 0.16, 0.175, 0.18, 0.2, 0.22, 0.225, 0.24, 0.25, 0.26, 0.275, 0.28, 0.3])
# ftg_reosizes_W = np.array([0.012, 0.016, 0.02, 0.024, 0.028, 0.032, 0.036, 0.04])
# ftg_reocts_W = np.array([0.1, 0.12, 0.125, 0.14, 0.15, 0.16, 0.175, 0.18, 0.2, 0.22, 0.225, 0.24, 0.25, 0.26, 0.275, 0.28, 0.3])

# def generate_valid_combinations():
#     count = 0
#     for l, w in itertools.product(ftg_lengths, ftg_widths):
#         if (COL_LENGTH <= COL_WIDTH and l <= w) or (COL_LENGTH > COL_WIDTH and l >= w):
#             for d, cs, rl, rcl, rw, rcw in itertools.product(ftg_depths, ftg_conc_strengths, ftg_reosizes_L, ftg_reocts_L, ftg_reosizes_W, ftg_reocts_W):
#                 count += 1
#                 if count == 3000000:
#                     return tuple(float(x) for x in (l, w, d, cs, rl, rcl, rw, rcw))
#     return None

# result = generate_valid_combinations()
# print(result)


import itertools
import numpy as np

COL_LENGTH = 0.5
COL_WIDTH = 0.4
FTG_DPTH_MIN = 1.5
FTG_COVER = 0.075

# Define your arrays
ftg_lengths = np.array([4.0, 4.05, 4.1, 4.15, 4.2, 4.25, 4.3, 4.35, 4.4, 4.45, 4.5])
ftg_widths = np.array([4.0, 4.05, 4.1, 4.15, 4.2, 4.25, 4.3, 4.35, 4.4, 4.45, 4.5])
ftg_depths = np.array([1.5, 1.55, 1.6, 1.65, 1.7, 1.75, 1.8, 1.85, 1.9, 1.95, 2.0, 2.05, 2.1, 2.15, 2.2, 2.25, 2.3, 2.35, 2.4, 2.45, 2.5])
ftg_conc_strengths = np.array([25.0, 32.0, 40.0, 50.0, 65.0])
ftg_reosizes_L = np.array([0.012, 0.016, 0.02, 0.024, 0.028, 0.032, 0.036, 0.04])
ftg_reocts_L = np.array([0.1, 0.12, 0.125, 0.14, 0.15, 0.16, 0.175, 0.18, 0.2, 0.22, 0.225, 0.24, 0.25, 0.26, 0.275, 0.28, 0.3])
ftg_reosizes_W = np.array([0.012, 0.016, 0.02, 0.024, 0.028, 0.032, 0.036, 0.04])
ftg_reocts_W = np.array([0.1, 0.12, 0.125, 0.14, 0.15, 0.16, 0.175, 0.18, 0.2, 0.22, 0.225, 0.24, 0.25, 0.26, 0.275, 0.28, 0.3])

def check_reinforcement_constraint(bar_sizes, spacings, depths):
    return 250000 * bar_sizes[:, np.newaxis, np.newaxis]**2 * np.pi / spacings[np.newaxis, :, np.newaxis] < 1019.65 * depths[np.newaxis, np.newaxis, :]**2 / (depths[np.newaxis, np.newaxis, :] - FTG_COVER)

def generate_valid_combinations():
    # Calculate total number of combinations
    total_combinations = (
        len(ftg_lengths) * len(ftg_widths) * len(ftg_depths) * len(ftg_conc_strengths) *
        len(ftg_reosizes_L) * len(ftg_reocts_L) * len(ftg_reosizes_W) * len(ftg_reocts_W)
    )
    print(f"Starting number of combinations: {total_combinations}")

    # Pre-compute valid length-width combinations
    valid_lw = [(l, w) for l, w in itertools.product(ftg_lengths, ftg_widths)
                if ((COL_LENGTH <= COL_WIDTH and l <= w) or (COL_LENGTH > COL_WIDTH and l >= w)) and
                   not (1 < l/w <= 1.1 or 0.9 <= l/w < 1)]

    # Pre-compute valid reinforcement combinations
    valid_reo_L = check_reinforcement_constraint(ftg_reosizes_L, ftg_reocts_L, ftg_depths)
    valid_reo_W = check_reinforcement_constraint(ftg_reosizes_W, ftg_reocts_W, ftg_depths)

    valid_combinations = []
    for l, w in valid_lw:
        for d_idx, d in enumerate(ftg_depths):
            valid_rl_rcl = np.argwhere(valid_reo_L[:, :, d_idx])
            valid_rw_rcw = np.argwhere(valid_reo_W[:, :, d_idx])
            for cs in ftg_conc_strengths:
                for (rl_idx, rcl_idx), (rw_idx, rcw_idx) in itertools.product(valid_rl_rcl, valid_rw_rcw):
                    valid_combinations.append((
                        float(l), float(w), float(d), float(cs),
                        float(ftg_reosizes_L[rl_idx]), float(ftg_reocts_L[rcl_idx]),
                        float(ftg_reosizes_W[rw_idx]), float(ftg_reocts_W[rcw_idx])
                    ))

    print(f"Final number of valid combinations: {len(valid_combinations)}")
    return valid_combinations

result = generate_valid_combinations()
print(f"Number of valid combinations: {len(result)}")