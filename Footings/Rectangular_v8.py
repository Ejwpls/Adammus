import math
from SizesGeneration import gen_ranges, gen_starting_values

class Foundation:
    def __init__(self, ftg_cover, col_length, col_width, load_dead, load_live, bearing_pressure):
        self.ftg_cover = ftg_cover
        self.col_length = col_length
        self.col_width = col_width
        self.load_dead = load_dead
        self.load_live = load_live
        self.bp = bearing_pressure
        self.fc = 25
        self.barL = 0.012
        self.ctsL = 0.1
        self.barW = 0.012
        self.ctsW = 0.1

        self.l, self.w, self.d = gen_starting_values(self.col_length, self.col_width, self.bp, self.load_dead, self.load_live)
        
        self.optimize()

    def optimize(self):
        size_increment = 0.05
        while True:
            # Check bearing capacity
            total_load = (
                self.load_dead
                + self.load_live
                + 6 * self.d * self.l * self.w
            )
            bearing_capacity = self.bp * self.l * self.w
            bearing_ratio = total_load / bearing_capacity

            # Calculate ex1L (new equation)
            Q = 4 * self.load_dead + 5 * self.load_live + 6 * self.d * self.l * self.w
            T = max((18 * self.d) / 25, -(9/20) * (2 * self.barL + self.barW + 2 * self.ftg_cover - 2 * self.d))
            P = max(-(9/20) * (self.barL + 2 * self.ftg_cover - 2 * self.d), (18 * self.d) / 25)
            
            ex1L = min(
                3/1000,
                (self.ctsL * ((3 * Q * (-self.col_length + self.l - 2 * T)) / (self.l * self.w) + 
                            (20 * max((3 * P * Q * (-self.col_length + self.l - 2 * T)) / (20 * self.l * self.w),
                                        (3 * Q * (self.col_length - self.l + 2 * T)**2) / (80 * self.l * self.w))) / P)) / 
                (2000000000 * self.barL**2 * math.pi)
            )

            ex1W = min(
                3/1000,
                (self.ctsW * ((3 * Q * (-self.col_width - 2 * T + self.w)) / (self.l * self.w) + 
                            (20 * max((3 * Q * (self.col_width + 2 * T - self.w)**2) / (80 * self.l * self.w),
                                        (3 * Q * T * (-self.col_width - 2 * T + self.w)) / (20 * self.l * self.w))) / T)) / 
                (2000000000 * self.barW**2 * math.pi)
            )

            # Check punching shear
            fVucL = ((1 + 1500 * ex1L) * 
                    (-((27 * self.d) / 5) + ((3 * self.load_live) / 2 + 6/5 * (self.load_dead + 6 * self.d * self.l * self.w)) / (self.l * self.w)) * 
                    (1 + max((18 * self.d) / 25, 9/10 * (-(self.barL/2) - self.ftg_cover + self.d))) * 
                    (-self.col_length + self.l - 2 * max((18 * self.d) / 25, 9/10 * (-self.barL - self.barW/2 - self.ftg_cover + self.d)))) / \
                    (728 * math.sqrt(self.fc) * max(1/2, (10 * (1 - self.d)) / 7) * max((18 * self.d) / 25, 9/10 * (-(self.barL/2) - self.ftg_cover + self.d)))

            fVucW = ((1 + 1500 * ex1W) * 
                    (-((27 * self.d) / 5) + ((3 * self.load_live) / 2 + 6/5 * (self.load_dead + 6 * self.d * self.l * self.w)) / (self.l * self.w)) * 
                    (-self.col_width + self.w - 2 * max((18 * self.d) / 25, 9/10 * (-self.barL - self.barW/2 - self.ftg_cover + self.d))) * 
                    (1 + max((18 * self.d) / 25, 9/10 * (-self.barL - self.barW/2 - self.ftg_cover + self.d)))) / \
                    (728 * math.sqrt(self.fc) * max(1/2, (10 * (1 - self.d)) / 7) * max((18 * self.d) / 25, 9/10 * (-self.barL - self.barW/2 - self.ftg_cover + self.d)))

            onewaysheer_ratio = max(fVucL, fVucW)

            # Check all conditions and adjust dimensions if necessary
            if bearing_ratio > 1:
                if self.col_length == self.col_width:
                    self.l += size_increment
                    self.w += size_increment
                elif (self.l - self.col_length) / 2 > (self.w - self.col_width) / 2:
                    self.w += size_increment
                else:
                    self.l += size_increment
            elif onewaysheer_ratio > 1:
                self.d += size_increment
            else:
                break  # All conditions are satisfied, exit the loop

        print(f"Optimized dimensions: L={self.l}, W={self.w}, D={self.d}")
        print(f"Final ratios: Bearing={bearing_ratio:.2f}, fVucL={fVucL:.2f}, fVucW={fVucW:.2f}")

# Example usage:
foundation = Foundation(
    ftg_cover=0.060,
    col_length=0.600,
    col_width=0.500,
    load_dead=4000,
    load_live=1500,
    bearing_pressure=250,
)

print(f"Optimized Foundation L: {foundation.l}")
print(f"Optimized Foundation W: {foundation.w}")
print(f"Optimized Foundation D: {foundation.d}")