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
        self.barL = 0.012
        self.barW = 0.012

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
            ex1L = self.calculate_ex1L()  # You need to implement this method

            # Check punching shear
            fVucL = ((1 + 1500 * ex1L) * 
                    (-((27 * self.d) / 5) + ((3 * self.load_live) / 2 + 6/5 * (self.load_dead + 6 * self.d * self.l * self.w)) / (self.l * self.w)) * 
                    (1 + max((18 * self.d) / 25, 9/10 * (-(self.barL/2) - self.ftg_cover + self.d))) * 
                    (-self.col_length + self.l - 2 * max((18 * self.d) / 25, 9/10 * (-self.barL - self.barW/2 - self.ftg_cover + self.d)))) / \
                    (728 * math.sqrt(self.fc) * max(1/2, (10 * (1 - self.d)) / 7) * max((18 * self.d) / 25, 9/10 * (-(self.barL/2) - self.ftg_cover + self.d)))

            fVucW = ((1 + 1500 * self.ex1W) * 
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
