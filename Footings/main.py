import time
import numpy as np
from SizesGeneration import (
    generate_foundation_size_ranges,
    gen_FDN_sizes,
    print_foundation_sizes
)

def find_suitable_foundation(combinations, LOAD_DEAD, LOAD_LIVE, BEARINGPRESSURE):
    L, W, D = combinations[:, 0], combinations[:, 1], combinations[:, 2]
    mask = (LOAD_DEAD + LOAD_LIVE + 6 * D * L * W) < BEARINGPRESSURE * L * W
    suitable_indices = np.where(mask)[0]
    return combinations[suitable_indices[0]] if suitable_indices.size > 0 else None

def main():
    start_time = time.time()

    COL_LENGTH: float = 0.5
    COL_WIDTH: float = 0.5
    FTG_COVER: float = 0.060
    BEARINGPRESSURE: float = 250
    LOAD_DEAD: float = 8000
    LOAD_LIVE: float = 1000

    # Generate arrays
    size_ranges = generate_foundation_size_ranges(
        COL_LENGTH, COL_WIDTH, BEARINGPRESSURE, LOAD_DEAD, LOAD_LIVE
    )

    # Generate valid combinations
    combinations = gen_FDN_sizes(COL_LENGTH, COL_WIDTH, size_ranges[0], size_ranges[1], size_ranges[2])
    
    # Find suitable foundation
    suitable_foundation = find_suitable_foundation(combinations, LOAD_DEAD, LOAD_LIVE, BEARINGPRESSURE)
    
    print(suitable_foundation)

    end_time = time.time()
    execution_time = end_time - start_time
    print(f"\nExecution time: {execution_time:.2f} seconds")

if __name__ == "__main__":
    main()