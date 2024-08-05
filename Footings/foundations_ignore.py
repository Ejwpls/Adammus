import numpy as np
import pandas as pd
import math
from dataclasses import dataclass
from typing import Callable, Dict, Any

np.seterr(all="raise")  # Error on overflow

@dataclass
class FoundationParameters:
    REO_DENSITY: float = 7850  # kg/m3
    LOAD_DEAD: float = 8000  # kN
    LOAD_LIVE: float = 2500  # kN
    BEARINGPRESSURE: float = 500  # kPa
    FTG_COVER: float = 0.060  # m
    COLUMN_WIDTH: float = 0.500  # m

@dataclass
class Calculation:
    func: Callable
    short_name: str = "Unnamed"
    long_name: str = "Unnamed"
    definition: str = "No definition"
    units: str = "-"

class FoundationCalculator:
    def __init__(self, params: FoundationParameters):
        self.params = params
        self.calculations: Dict[str, Calculation] = {}

    def add_calculation(self, name: str, calc: Calculation):
        self.calculations[name] = calc

    def calculate(self, sizes: np.ndarray) -> np.ndarray:
        results = {name: sizes[name] for name in sizes.dtype.names}
        
        while len(results) < len(self.calculations) + len(sizes.dtype.names):
            for name, calc in self.calculations.items():
                if name not in results:
                    args = {arg: results[arg] for arg in calc.func.__code__.co_varnames if arg in results}
                    if len(args) == len(calc.func.__code__.co_varnames):
                        results[name] = calc.func(**args)
        
        new_dtype = sizes.dtype.descr + [(name, 'f8') for name in results if name not in sizes.dtype.names]
        new_sizes = np.empty(sizes.shape, dtype=new_dtype)
        for name in sizes.dtype.names:
            new_sizes[name] = sizes[name]
        for name, result in results.items():
            if name not in sizes.dtype.names:
                new_sizes[name] = result
        
        return new_sizes

class FoundationSizer:
    def __init__(self, params: FoundationParameters, calculator: FoundationCalculator):
        self.params = params
        self.calculator = calculator
        self.sizes = self.generate_foundation_sizes()
        self.calculate_all_properties()
        self.column_names = list(self.sizes.dtype.names)

    def generate_foundation_sizes(self) -> np.ndarray:
        FTG_LEN_MIN = (
            math.ceil(
                math.sqrt(
                    (self.params.LOAD_DEAD + self.params.LOAD_LIVE)
                    / self.params.BEARINGPRESSURE
                )
                / 0.05
            )
            * 0.05
        )
        FTG_LEN_MAX = FTG_LEN_MIN + 1
        FTG_LENS = np.round(
            np.arange(FTG_LEN_MIN, FTG_LEN_MAX + 0.001, 0.05, dtype=np.float32), 2
        )

        FTG_DPTH_MIN = (
            math.ceil(
                (
                    4
                    * math.sqrt(3570)
                    * math.sqrt(
                        (4760 * self.params.COLUMN_WIDTH**2) / 3
                        + (1 + (3 * self.params.BEARINGPRESSURE) / 19040)
                        * (self.params.LOAD_DEAD + self.params.LOAD_LIVE)
                    )
                    - (9520 + 3 * self.params.BEARINGPRESSURE)
                    * self.params.COLUMN_WIDTH
                )
                / (19040 + 3 * self.params.BEARINGPRESSURE)
                / 0.05
            )
            * 0.05
        )
        FTG_DPTH_MAX = FTG_DPTH_MIN * 2
        FTG_DPTHS = np.round(
            np.arange(FTG_DPTH_MIN, FTG_DPTH_MAX + 0.001, 0.05, dtype=np.float32), 2
        )

        FTG_CONC_STRENGTHS = np.array([20, 25, 32, 40, 50, 65], dtype=np.float32)

        FTG_REO_SIZES = np.round(
            np.array(
                [0.012, 0.016, 0.02, 0.024, 0.028, 0.032, 0.036, 0.04],
                dtype=np.float32,
            ),
            3,
        )

        FTG_REO_CTS = np.unique(
            np.round(
                np.concatenate(
                    [
                        np.arange(0.1, 0.301, 0.025, dtype=np.float32),
                        np.arange(0.08, 0.301, 0.02, dtype=np.float32),
                    ]
                ),
                3,
            )
        )

        foundation_sizes = np.array(
            np.meshgrid(
                FTG_LENS, FTG_DPTHS, FTG_CONC_STRENGTHS, FTG_REO_SIZES, FTG_REO_CTS
            )
        ).T.reshape(-1, 5)

        return np.rec.array(
            foundation_sizes,
            dtype=[
                ("FtgLength", "f4"),
                ("FtgDepth", "f4"),
                ("fc", "f4"),
                ("ReoSize", "f4"),
                ("ReoCts", "f4"),
            ],
        )

    def calculate_all_properties(self) -> None:
        self.sizes = self.calculator.calculate(self.sizes)

    def remove_fails(self) -> None:
        valid_mask = (
            (self.sizes["Mrat"] <= 1)
            & (self.sizes["VLrat"] <= 1)
            & (self.sizes["VPrat"] <= 1)
            & (self.sizes["Astact"] >= self.sizes["Ast"])
            & (self.sizes["BPrat"] <= 1)
        )
        self.sizes = self.sizes[valid_mask]

    def filter_match(self, **kwargs) -> None:
        mask = np.ones(len(self.sizes), dtype=bool)
        for key, value in kwargs.items():
            if key in self.column_names:
                mask &= self.sizes[key] == value
        self.sizes = self.sizes[mask]

    def sort_by_cost(self) -> None:
        self.sizes = np.sort(self.sizes, order="Cost")
    
    def print_foundation_details(self, row_index: int) -> None:
        if row_index < 0 or row_index >= len(self.sizes):
            print(f"Row index {row_index} is out of range.")
            return

        row = self.sizes[row_index]
        max_name_length = max(len(f"{calc.long_name} ({calc.short_name})") for calc in self.calculator.calculations.values())
        max_value_length = max(len(f"{row[name]:.3f}" if isinstance(row[name], float) else str(row[name])) for name in self.calculator.calculations)

        print(f"Details for foundation size at row {row_index}:")
        for name, calc in self.calculator.calculations.items():
            value = row[name]
            formatted_value = f"{value:.3f}" if isinstance(value, float) else str(value)
            print(f"{calc.long_name} ({calc.short_name}):{' ' * (max_name_length - len(f'{calc.long_name} ({calc.short_name})')+5)}{formatted_value:>{max_value_length}} {calc.units}")

# Usage
params = FoundationParameters()
calculator = FoundationCalculator(params)

# Add calculations
calculator.add_calculation("Pult", Calculation(
    func=lambda FtgDepth, FtgLength: (
        1.2 * (6 * FtgDepth * FtgLength**2 + params.LOAD_DEAD) + 1.5 * params.LOAD_LIVE
    ),
    short_name="Pu",
    long_name="Ultimate Load",
    units="kN"
))

calculator.add_calculation("BPmax", Calculation(
    func=lambda FtgDepth, FtgLength: (
        6 * FtgDepth * FtgLength**2 + params.LOAD_LIVE + params.LOAD_DEAD
    ) / (FtgLength**2),
    short_name="BPm",
    long_name="Maximum Bearing Pressure",
    units="kPa"
))

calculator.add_calculation("BPrat", Calculation(
    func=lambda BPmax: BPmax / params.BEARINGPRESSURE,
    short_name="BPr",
    long_name="Bearing Pressure Ratio",
    units=""
))

calculator.add_calculation("Dom", Calculation(
    func=lambda FtgDepth, ReoSize: FtgDepth - params.FTG_COVER - ReoSize / 2,
    short_name="d",
    long_name="Effective Depth",
    units="m"
))

calculator.add_calculation("BPult", Calculation(
    func=lambda Pult, FtgLength: Pult / FtgLength**2,
    short_name="BPu",
    long_name="Ultimate Bearing Pressure",
    units="kPa"
))

calculator.add_calculation("CLR", Calculation(
    func=lambda BPult, Dom: BPult * (params.COLUMN_WIDTH + Dom) ** 2,
    short_name="CLR",
    long_name="Column Load Reaction",
    units="kN"
))

calculator.add_calculation("VPult", Calculation(
    func=lambda Pult, CLR: Pult - CLR,
    short_name="VPu",
    long_name="Ultimate Punching Shear",
    units="kN"
))

calculator.add_calculation("fVP", Calculation(
    func=lambda ReoSize, Dom, fc: (
        952 * (ReoSize - Dom) * (ReoSize - Dom - params.COLUMN_WIDTH) * np.sqrt(fc)
    ),
    short_name="fVP",
    long_name="Punching Shear Capacity",
    units="kN"
))

calculator.add_calculation("VPrat", Calculation(
    func=lambda VPult, fVP: VPult / fVP,
    short_name="VPr",
    long_name="Punching Shear Ratio",
    units=""
))

calculator.add_calculation("dv", Calculation(
    func=lambda Dom, FtgDepth: np.maximum(0.9 * Dom, 0.72 * FtgDepth),
    short_name="dv",
    long_name="Effective Shear Depth",
    units="m"
))

calculator.add_calculation("VLult", Calculation(
    func=lambda BPult, dv, FtgLength: -0.5
    * BPult
    * (params.COLUMN_WIDTH + 2 * dv - FtgLength),
    short_name="VLu",
    long_name="Ultimate Linear Shear",
    units="kN"
))

calculator.add_calculation("kvo", Calculation(
    func=lambda dv: 2 / (10 + 13 * dv),
    short_name="kvo",
    long_name="Initial Shear Factor",
    units=""
))

calculator.add_calculation("kv", Calculation(
    func=lambda kvo: np.minimum(kvo, 0.15),
    short_name="kv",
    long_name="Shear Factor",
    units=""
))

calculator.add_calculation("ks", Calculation(
    func=lambda FtgDepth: np.maximum(0.5, (10 / 7) * (1 - FtgDepth)),
    short_name="ks",
    long_name="Size Factor",
    units=""
))

calculator.add_calculation("fVuc", Calculation(
    func=lambda dv, fc, ks, kv: 700 * dv * np.sqrt(fc) * ks * kv,
    short_name="fVuc",
    long_name="Concrete Shear Strength",
    units="kN"
))

calculator.add_calculation("VLrat", Calculation(
    func=lambda VLult, fVuc: VLult / fVuc,
    short_name="VLr",
    long_name="Linear Shear Ratio",
    units=""
))

calculator.add_calculation("Mult", Calculation(
    func=lambda BPult, FtgLength: (BPult * (7 * params.COLUMN_WIDTH - 10 * FtgLength) ** 2) / 800,
    short_name="Mu",
    long_name="Ultimate Moment",
    units="kN·m"
))

calculator.add_calculation("Astshr", Calculation(
    func=lambda Mult, Dom: (5 * Mult) / (2 * Dom),
    short_name="Ass",
    long_name="Steel Area for Shear",
    units="mm²"
))

calculator.add_calculation("Astmin", Calculation(
    func=lambda FtgDepth, fc, Dom: (228 * FtgDepth**2 * np.sqrt(fc)) / Dom,
    short_name="Asm",
    long_name="Minimum Steel Area",
    units="mm²"
))

calculator.add_calculation("Ast", Calculation(
    func=lambda Astmin, Astshr: np.maximum(Astmin, Astshr),
    short_name="As",
    long_name="Required Steel Area",
    units="mm²"
))

calculator.add_calculation("alpha", Calculation(
    func=lambda fc: 0.85 - 0.0015 * fc,
    short_name="α",
    long_name="Stress Block Factor",
    units=""
))

calculator.add_calculation("gamma", Calculation(
    func=lambda fc: 0.97 - 0.0025 * fc,
    short_name="γ",
    long_name="Stress Block Depth Factor",
    units=""
))

calculator.add_calculation("ku", Calculation(
    func=lambda Ast, alpha, Dom, fc, gamma, FtgLength: Ast
    / (2000 * alpha * Dom * fc * gamma * FtgLength),
    short_name="ku",
    long_name="Neutral Axis Parameter",
    units=""
))

calculator.add_calculation("phi", Calculation(
    func=lambda ku: np.minimum(0.85, np.maximum(0.65, 1.24 - 1.08333 * ku)),
    short_name="φ",
    long_name="Capacity Reduction Factor",
    units=""
))

calculator.add_calculation("Astact", Calculation(
    func=lambda ReoSize, ReoCts: (250000 * ReoSize**2 * np.pi) / ReoCts,
    short_name="Asa",
    long_name="Actual Steel Area",
    units="mm²"
))

calculator.add_calculation("fMuo", Calculation(
    func=lambda Astact, Dom, phi, alpha, fc, FtgLength: (Astact * Dom * phi) / 2
    - (Astact**2 * phi) / (8000 * alpha * fc * FtgLength),
    short_name="fMuo",
    long_name="Moment Capacity",
    units="kN·m"
))

calculator.add_calculation("Mrat", Calculation(
    func=lambda Mult, fMuo: Mult / fMuo,
    short_name="Mr",
    long_name="Moment Ratio",
    units=""
))

calculator.add_calculation("Cost", Calculation(
    func=lambda Astact, FtgLength, FtgDepth, fc: (
        Astact / 1000000 * FtgLength * params.REO_DENSITY * 2 * 3400
        + FtgLength**2 * FtgDepth * (130.866 * np.exp(fc * 0.0111) + 45 + 130)
        + 4 * FtgDepth * FtgLength * 180
    ),
    short_name="Cost",
    long_name="Total Cost",
    units="$"
))

# Initialize FoundationSizer and perform calculations
foundation = FoundationSizer(params, calculator)

# Remove fails, sort by cost, and print results
foundation.remove_fails()
foundation.sort_by_cost()
foundation.print_foundation_details(0)