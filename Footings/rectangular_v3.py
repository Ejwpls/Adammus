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

def create_hash(l, w, d, cs):
    # Create a simple hash using the first two decimal places of each parameter
    return f"{int(l*100):04d}{int(w*100):04d}{int(d*100):04d}{int(cs):02d}"

def check_reinforcement_constraint(bar_sizes, spacings, depths):
    return 250000 * bar_sizes[:, np.newaxis, np.newaxis]**2 * np.pi / spacings[np.newaxis, :, np.newaxis] < 1019.65 * depths[np.newaxis, np.newaxis, :]**2 / (depths[np.newaxis, np.newaxis, :] - FTG_COVER)

def generate_valid_combinations():
    valid_combinations = []
    count = 0

    for l, w in itertools.product(ftg_lengths, ftg_widths):
        if ((COL_LENGTH <= COL_WIDTH and l <= w) or (COL_LENGTH > COL_WIDTH and l >= w)) and not (1 < l/w <= 1.1 or 0.9 <= l/w < 1) and l*w>=20:
            for d, cs, rl, rcl, rw, rcw in itertools.product(ftg_depths, ftg_conc_strengths, ftg_reosizes_L, ftg_reocts_L, ftg_reosizes_W, ftg_reocts_W):
                if (rl**2 / rcl < 0.00129826 * d**2 / (d - FTG_COVER) and
                    rw**2 / rcw < 0.00129826 * d**2 / (d - FTG_COVER)):
                    # Create a hash for the combination
                    combo_hash = create_hash(l, w, d, cs)
                    valid_combinations.append((combo_hash, float(l), float(w), float(d), float(cs), float(rl), float(rcl), float(rw), float(rcw)))
                    count += 1
                    if count == 100000:
                        return valid_combinations

    print(f"Found {count} valid combinations")
    return valid_combinations

result = generate_valid_combinations()
print(f"Number of valid combinations: {len(result)}")

# Print the first 10 rows of working foundations
print("\nFirst 10 rows of working foundations:")
print("Hash     | Length | Width | Depth | fc | Reo Size L | Reo Spacing L | Reo Size W | Reo Spacing W")
print("-" * 100)
for i, combo in enumerate(result[:10]):
    print(f"{combo[0]} | {combo[1]:.2f}   | {combo[2]:.2f}  | {combo[3]:.2f}  | {combo[4]:.0f} | {combo[5]:.3f}      | {combo[6]:.3f}         | {combo[7]:.3f}      | {combo[8]:.3f}")