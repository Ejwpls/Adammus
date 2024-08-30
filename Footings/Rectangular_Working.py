import numpy as np
import pandas as pd
import math
from dataclasses import dataclass
from itertools import product
from typing import Callable, Dict, Optional

np.seterr(all="raise")  # Error on overflow


@dataclass
class FoundationParameters:
    REO_DENSITY: float = 7850  # kg/m3
    LOAD_DEAD: float = 8000  # kN
    LOAD_LIVE: float = 2500  # kN
    BEARINGPRESSURE: float = 500  # kPa
    FTG_COVER: float = 0.060  # m
    COLUMN_LENGTH: float = 0.500  # m
    COLUMN_WIDTH: float = 0.500  # m


@dataclass
class Calculation:
    func: Callable
    short_name: str
    long_name: str
    definition: str = ""
    format_func: Optional[Callable[[float], str]] = None

    def format_value(self, value: float) -> str:
        if self.format_func:
            return self.format_func(value)
        elif isinstance(value, float):
            return f"{value:.3f}"
        else:
            return str(value)


class FoundationCalculator:
    def __init__(self, params: FoundationParameters):
        self.params = params
        self.calculations: Dict[str, Calculation] = {}

    def add_calculation(self, name: str, calc: Calculation):
        self.calculations[name] = calc

    def calculate(self, sizes: np.ndarray) -> np.ndarray:
        results = {name: sizes[name] for name in sizes.dtype.names}
        calculated = set(sizes.dtype.names)

        for name, calc in self.calculations.items():
            args = calc.func.__code__.co_varnames
            missing_args = [arg for arg in args if arg not in results]

            if not missing_args:
                results[name] = calc.func(**{arg: results[arg] for arg in args})
                calculated.add(name)
            else:
                print(f"Skipping {name}. Missing arguments: {', '.join(missing_args)}")

        uncalculated = set(self.calculations.keys()) - calculated
        if len(uncalculated) > 0:
            raise ValueError(
                f"Unable to calculate all fields. Uncalculated fields: {', '.join(uncalculated)}"
            )

        new_dtype = sizes.dtype.descr + [
            (name, "f8") for name in results if name not in sizes.dtype.names
        ]
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
        self.properties = self.define_properties()

    def define_properties(self):
        return [
            ("FtgLength", "L", "Foundation Length", lambda x: f"{x*1000:.0f} mm"),
            ("FtgWidth", "W", "Foundation Width", lambda x: f"{x*1000:.0f} mm"),
            ("FtgDepth", "D", "Foundation Depth", lambda x: f"{x*1000:.0f} mm"),
            ("fc", "fc", "Concrete Strength", lambda x: f"{x:.0f} MPa"),
            ("ReoSizeL", "Bar L", "Reinforcement Size (L)", lambda x: f"N{x*1000:.0f}"),
            (
                "ReoCtsL",
                "CTS L",
                "Reinforcement Spacing (L)",
                lambda x: f"{x*1000:.0f} mm",
            ),
            ("ReoSizeW", "Bar W", "Reinforcement Size (W)", lambda x: f"N{x*1000:.0f}"),
            (
                "ReoCtsW",
                "CTS W",
                "Reinforcement Spacing (W)",
                lambda x: f"{x*1000:.0f} mm",
            ),
        ]

    def generate_foundation_sizes(self) -> np.ndarray:
        # Calculate starting length and width
        START_WIDTH = (
            self.params.COLUMN_WIDTH - self.params.COLUMN_LENGTH
        ) / 2 + math.sqrt(
            self.params.BEARINGPRESSURE
            * (
                self.params.BEARINGPRESSURE
                * (self.params.COLUMN_LENGTH - self.params.COLUMN_WIDTH) ** 2
                + 4 * (self.params.LOAD_DEAD + self.params.LOAD_LIVE)
            )
        ) / (2 * self.params.BEARINGPRESSURE)
        START_LENGTH = (
            self.params.COLUMN_LENGTH - self.params.COLUMN_WIDTH
        ) / 2 + math.sqrt(
            self.params.BEARINGPRESSURE
            * (
                self.params.BEARINGPRESSURE
                * (self.params.COLUMN_LENGTH - self.params.COLUMN_WIDTH) ** 2
                + 4 * (self.params.LOAD_DEAD + self.params.LOAD_LIVE)
            )
        ) / (2 * self.params.BEARINGPRESSURE)

        # Round up to nearest 0.05
        FTG_LENGTH_MIN = math.ceil(START_LENGTH / 0.05) * 0.05
        FTG_WIDTH_MIN = math.ceil(START_WIDTH / 0.05) * 0.05

        FTG_LENGTH_MAX = FTG_LENGTH_MIN + 0.5
        FTG_WIDTH_MAX = FTG_WIDTH_MIN + 0.5
        FTG_LENGTHS = np.round(
            np.arange(FTG_LENGTH_MIN, FTG_LENGTH_MAX + 0.001, 0.05, dtype=np.float32), 2
        )
        FTG_WIDTHS = np.round(
            np.arange(FTG_WIDTH_MIN, FTG_WIDTH_MAX + 0.001, 0.05, dtype=np.float32), 2
        )

        FTG_DPTH_MIN = (
            math.ceil(
                (
                    -(self.params.COLUMN_LENGTH / 4)
                    - self.params.COLUMN_WIDTH / 4
                    + 1
                    / 20
                    * math.sqrt(
                        25 * (self.params.COLUMN_LENGTH + self.params.COLUMN_WIDTH) ** 2
                        + 3
                        / 238
                        * math.sqrt(5)
                        * (4 * self.params.LOAD_DEAD + 5 * self.params.LOAD_LIVE)
                    )
                )
                / 0.05
            )
            * 0.05
        )

        FTG_DPTH_MAX = FTG_DPTH_MIN + 0.5
        FTG_DPTHS = np.round(
            np.arange(FTG_DPTH_MIN, FTG_DPTH_MAX + 0.001, 0.05, dtype=np.float32), 2
        )

        FTG_CONC_STRENGTHS = np.array([25, 32, 40, 50, 65], dtype=np.float32)

        FTG_REO_SIZES = np.round(
            np.array(
                [0.012, 0.016, 0.02, 0.024, 0.028, 0.032, 0.036, 0.04], dtype=np.float32
            ),
            3,
        )

        FTG_REO_CTS = np.unique(
            np.round(
                np.concatenate(
                    [
                        np.arange(0.1, 0.301, 0.025, dtype=np.float32),
                        np.arange(0.1, 0.301, 0.02, dtype=np.float32),
                    ]
                ),
                3,
            )
        )

        # Calculate the area of steel for each combination of reo size and spacing
        def calculate_steel_area(bar_size, spacing):
            return 250000 * bar_size**2 * np.pi / spacing

        # Generate combinations for length and width directions
        reo_combinations = np.array(
            list(product(FTG_REO_SIZES, FTG_REO_CTS)), dtype=np.float32
        )

        # Vectorized calculation of steel areas
        steel_areas = calculate_steel_area(
            reo_combinations[:, 0], reo_combinations[:, 1]
        )

        # Filter out combinations where area of steel is less than the minimum required
        valid_reo_combinations = reo_combinations[
            steel_areas
            > 1019.65 * FTG_DPTH_MIN**2 / (FTG_DPTH_MIN - self.params.FTG_COVER)
        ]

        FOUNDATION_SIZE_DTYPE = np.dtype(
            [
                ("FtgLength", "f4"),
                ("FtgWidth", "f4"),
                ("FtgDepth", "f4"),
                ("fc", "f4"),
                ("ReoSizeL", "f4"),
                ("ReoCtsL", "f4"),
                ("ReoSizeW", "f4"),
                ("ReoCtsW", "f4"),
            ]
        )
        n_reo = len(valid_reo_combinations) ** 2
        dummy = np.arange(n_reo)
        # Create base combinations
        base_combinations = (
            np.array(
                np.meshgrid(
                    FTG_LENGTHS,
                    FTG_WIDTHS,
                    FTG_DPTHS,
                    FTG_CONC_STRENGTHS,
                    dummy,
                    indexing="ij",
                )
            )
            .reshape(5, -1)
            .T
        )
        total_len = (
            len(FTG_LENGTHS)
            * len(FTG_WIDTHS)
            * len(FTG_DPTHS)
            * len(FTG_CONC_STRENGTHS)
        )
        replacement = np.c_[
            valid_reo_combinations.repeat(len(valid_reo_combinations), 0),
            np.tile(valid_reo_combinations, (len(valid_reo_combinations), 1)),
        ]
        replacement = np.tile(replacement, (total_len, 1))
        result = np.column_stack((base_combinations[:, :-1], replacement))

        structured_array = np.empty(result.shape[0], dtype=FOUNDATION_SIZE_DTYPE)
        structured_array["FtgLength"] = result[:, 0]
        structured_array["FtgWidth"] = result[:, 1]
        structured_array["FtgDepth"] = result[:, 2]
        structured_array["fc"] = result[:, 3]
        structured_array["ReoSizeL"] = result[:, 4]
        structured_array["ReoCtsL"] = result[:, 5]
        structured_array["ReoSizeW"] = result[:, 6]
        structured_array["ReoCtsW"] = result[:, 7]

        return structured_array

    def calculate_all_properties(self) -> None:
        self.sizes = self.calculator.calculate(self.sizes)

    def remove_fails(self) -> None:
        valid_mask = (
            (self.sizes["Mrat"] <= 1)
            & (self.sizes["VLrat"] <= 1)
            & (self.sizes["VPrat"] <= 1)
            & (self.sizes["Astact"] >= self.sizes["Astreq"])
            & (self.sizes["Astact"] >= self.sizes["AstreqW"])
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
        max_name_length = max(
            len(calc.long_name) for calc in self.calculator.calculations.values()
        )
        column_width = max_name_length + 15  # Add some extra space for the short name

        print(f"Details for foundation size at row {row_index}:")

        # Print the foundation properties first
        for field, short_name, long_name, format_func in self.properties:
            value = row[field]
            formatted_value = format_func(value)
            print(
                "{:<{width}}  {}".format(
                    f"{long_name} ({short_name}):", formatted_value, width=column_width
                )
            )

        # Print the calculations
        for name, calc in self.calculator.calculations.items():
            value = row[name]
            formatted_value = calc.format_value(value)
            print(
                "{:<{width}}  {}".format(
                    f"{calc.long_name} ({calc.short_name}):",
                    formatted_value,
                    width=column_width,
                )
            )

    def print_array(self, num_rows: int) -> None:
        if num_rows <= 0:
            print("Number of rows must be greater than 0.")
            return

        if num_rows > len(self.sizes):
            print(f"Only {len(self.sizes)} rows available. Printing all rows.")
            num_rows = len(self.sizes)

        # Determine which fields to display
        display_fields = [prop[0] for prop in self.properties] + [
            "SWt"
        ]  # Add other fields as needed

        # Create headers
        headers = []
        for field in display_fields:
            prop = next((p for p in self.properties if p[0] == field), None)
            if prop:
                headers.append(prop[1])  # Use short name
            else:
                calc = self.calculator.calculations.get(field)
                if calc:
                    headers.append(calc.short_name)
                else:
                    headers.append(field)  # Fallback to field name

        # Print column headers
        header = " | ".join(f"{header:^15}" for header in headers)
        print(header)
        print("-" * len(header))

        # Print rows
        for i in range(num_rows):
            row = self.sizes[i]
            formatted_row = []
            for field in display_fields:
                value = row[field]
                prop = next((p for p in self.properties if p[0] == field), None)
                if prop:
                    formatted_value = prop[3](value)  # Use the format function
                else:
                    calc = self.calculator.calculations.get(field)
                    if calc:
                        formatted_value = calc.format_value(value)
                    else:
                        formatted_value = f"{value}"
                formatted_row.append(formatted_value)
            print(" | ".join(f"{value:^15}" for value in formatted_row))


# Usage
params = FoundationParameters()
calculator = FoundationCalculator(params)

# Add calculations
calculator.add_calculation(
    "SWt",
    Calculation(
        func=lambda FtgDepth, FtgLength, FtgWidth: (
            6 * FtgDepth * FtgLength * FtgWidth
        ),
        short_name="SWt",
        long_name="Self Weight of Footing",
        format_func=lambda x: f"{x:.1f} kN",
    ),
)

calculator.add_calculation(
    "Pult",
    Calculation(
        func=lambda SWt: 1.2 * (params.LOAD_DEAD + SWt) + 1.5 * params.LOAD_LIVE,
        short_name="P*",
        long_name="Ultimate Load",
        format_func=lambda x: f"{x:.1f} kN",
    ),
)

calculator.add_calculation(
    "BPmax",
    Calculation(
        func=lambda FtgLength, FtgWidth, SWt: (
            params.LOAD_DEAD + SWt + params.LOAD_LIVE
        )
        / (FtgLength * FtgWidth),
        short_name="BPm",
        long_name="Maximum Bearing Pressure",
        format_func=lambda x: f"{x:.0f} kPa",
    ),
)

calculator.add_calculation(
    "BPult",
    Calculation(
        func=lambda Pult, FtgLength, FtgWidth: Pult / (FtgLength * FtgWidth),
        short_name="BP*",
        long_name="Ultimate Bearing Pressure",
        format_func=lambda x: f"{x:.0f} kPa",
    ),
)

calculator.add_calculation(
    "AstL",
    Calculation(
        func=lambda ReoSizeL, ReoCtsL: 250000 / ReoCtsL * ReoSizeL**2 * np.pi,
        short_name="AstL",
        long_name="Actual Steel Area (L)",
        format_func=lambda x: f"{x:.0f} mm²/m",
    ),
)

calculator.add_calculation(
    "AstW",
    Calculation(
        func=lambda ReoSizeW, ReoCtsW: 250000 / ReoCtsW * ReoSizeW**2 * np.pi,
        short_name="AstL",
        long_name="Actual Steel Area (W)",
        format_func=lambda x: f"{x:.0f} mm²/m",
    ),
)

calculator.add_calculation(
    "dsL",
    Calculation(
        func=lambda FtgDepth, ReoSizeL: FtgDepth - params.FTG_COVER - ReoSizeL / 2,
        short_name="ds",
        long_name="Effective Depth (L)",
        format_func=lambda x: f"{x*1000:.0f} mm",
    ),
)

calculator.add_calculation(
    "dsW",
    Calculation(
        func=lambda FtgDepth, ReoSizeL, ReoSizeW: FtgDepth
        - params.FTG_COVER
        - ReoSizeL
        - ReoSizeW / 2,
        short_name="dsW",
        long_name="Effective Depth (W)",
        format_func=lambda x: f"{x*1000:.0f} mm",
    ),
)

calculator.add_calculation(
    "AstminL",
    Calculation(
        func=lambda FtgDepth, fc, dsL: (228 * FtgDepth**2 * np.sqrt(fc)) / dsL,
        short_name="AstminL",
        long_name="Minimum Steel Area (L)",
        format_func=lambda x: f"{x:.0f} mm²/m",
    ),
)

calculator.add_calculation(
    "AstminW",
    Calculation(
        func=lambda FtgDepth, fc, dsW: (228 * FtgDepth**2 * np.sqrt(fc)) / dsW,
        short_name="AstminW",
        long_name="Minimum Steel Area (W)",
        format_func=lambda x: f"{x:.0f} mm²/m",
    ),
)

calculator.add_calculation(
    "alpha",
    Calculation(
        func=lambda fc: 0.85 - 0.0015 * fc,
        short_name="α",
        long_name="Alpha Factor",
        format_func=lambda x: f"{x:.3f}",
    ),
)

calculator.add_calculation(
    "gamma",
    Calculation(
        func=lambda fc: 0.97 - 0.0025 * fc,
        short_name="γ",
        long_name="Gamma Factor",
        format_func=lambda x: f"{x:.3f}",
    ),
)

calculator.add_calculation(
    "MultL",
    Calculation(
        func=lambda FtgLength, SWt, BPult, FtgWidth: (
            (7 * params.COLUMN_LENGTH - 10 * FtgLength) ** 2
            * (-9 * SWt + 10 * BPult * FtgLength * FtgWidth)
        )
        / (8000 * FtgLength * FtgWidth),
        short_name="M*L",
        long_name="Ultimate Moment (L)",
        format_func=lambda x: f"{x:.1f} kNm",
    ),
)

calculator.add_calculation(
    "MultW",
    Calculation(
        func=lambda FtgWidth, SWt, BPult, FtgLength: (
            (7 * params.COLUMN_WIDTH - 10 * FtgWidth) ** 2
            * (-9 * SWt + 10 * BPult * FtgLength * FtgWidth)
        )
        / (8000 * FtgLength * FtgWidth),
        short_name="M*W",
        long_name="Ultimate Moment (W)",
        format_func=lambda x: f"{x:.1f} kNm",
    ),
)

calculator.add_calculation(
    "AstshrL",
    Calculation(
        func=lambda MultL, dsL: 5 * MultL / (2 * dsL),
        short_name="AstshrL",
        long_name="Shrinkage Steel Area (L)",
        format_func=lambda x: f"{x:.0f} mm²/m",
    ),
)

calculator.add_calculation(
    "AstshrW",
    Calculation(
        func=lambda MultW, dsW: 5 * MultW / (2 * dsW),
        short_name="AstshrW",
        long_name="Shrinkage Steel Area (W)",
        format_func=lambda x: f"{x:.0f} mm²/m",
    ),
)

calculator.add_calculation(
    "AstreqL",
    Calculation(
        func=lambda AstminL, AstshrL: np.maximum(AstminL, AstshrL),
        short_name="AstreqL",
        long_name="Required Steel Area (L)",
        format_func=lambda x: f"{x:.0f} mm²/m",
    ),
)

calculator.add_calculation(
    "AstreqW",
    Calculation(
        func=lambda AstminW, AstshrW: np.maximum(AstminW, AstshrW),
        short_name="AstreqW",
        long_name="Required Steel Area (W)",
        format_func=lambda x: f"{x:.0f} mm²/m",
    ),
)

calculator.add_calculation(
    "kuL",
    Calculation(
        func=lambda AstreqL, alpha, dsL, fc, gamma, FtgLength: AstreqL
        / (2000 * alpha * dsL * fc * gamma * FtgLength),
        short_name="kuL",
        long_name="Neutral Axis Parameter (L)",
        format_func=lambda x: f"{x:.3f}",
    ),
)

calculator.add_calculation(
    "kuW",
    Calculation(
        func=lambda AstreqW, alpha, dsW, fc, gamma, FtgWidth: AstreqW
        / (2000 * alpha * dsW * fc * gamma * FtgWidth),
        short_name="kuW",
        long_name="Neutral Axis Parameter (W)",
        format_func=lambda x: f"{x:.3f}",
    ),
)

calculator.add_calculation(
    "phiL",
    Calculation(
        func=lambda kuL: np.minimum(0.85, np.maximum(0.65, 1.24 - 13 * kuL / 12)),
        short_name="φL",
        long_name="Capacity Reduction Factor (L)",
        format_func=lambda x: f"{x:.3f}",
    ),
)

calculator.add_calculation(
    "phiW",
    Calculation(
        func=lambda kuW: np.minimum(0.85, np.maximum(0.65, 1.24 - 13 * kuW / 12)),
        short_name="φW",
        long_name="Capacity Reduction Factor (W)",
        format_func=lambda x: f"{x:.3f}",
    ),
)

calculator.add_calculation(
    "fMuoL",
    Calculation(
        func=lambda AstL, dsL, phiL, alpha, fc: (AstL * dsL * phiL) / 2
        - (AstL**2 * phiL) / (8000 * alpha * fc),
        short_name="φMuoL",
        long_name="Design Moment Capacity (L)",
        format_func=lambda x: f"{x:.1f} kNm",
    ),
)

calculator.add_calculation(
    "fMuoW",
    Calculation(
        func=lambda AstW, dsW, phiW, alpha, fc: (AstW * dsW * phiW) / 2
        - (AstW**2 * phiW) / (8000 * alpha * fc),
        short_name="φMuoW",
        long_name="Design Moment Capacity (W)",
        format_func=lambda x: f"{x:.1f} kNm",
    ),
)

calculator.add_calculation(
    "CLR",
    Calculation(
        func=lambda BPult, dsL: BPult
        * (params.COLUMN_WIDTH + dsL)
        * (params.COLUMN_LENGTH + dsL),
        short_name="CLR",
        long_name="Critical Load Reaction",
        format_func=lambda x: f"{x:.1f} kN",
    ),
)

calculator.add_calculation(
    "VPult",
    Calculation(
        func=lambda Pult, CLR: Pult - CLR,
        short_name="VP*",
        long_name="Ultimate Punching Shear",
        format_func=lambda x: f"{x:.1f} kN",
    ),
)

calculator.add_calculation(
    "fcv",
    Calculation(
        func=lambda fc: 0.17
        * (
            1
            + 2
            / max(
                params.COLUMN_LENGTH / params.COLUMN_WIDTH,
                params.COLUMN_WIDTH / params.COLUMN_LENGTH,
                2,
            )
        )
        * np.sqrt(fc),
        short_name="fcv",
        long_name="Concrete Shear Strength",
        format_func=lambda x: f"{x:.2f} MPa",
    ),
)

calculator.add_calculation(
    "fVP",
    Calculation(
        func=lambda dsW, fcv: 1400
        * dsW
        * (params.COLUMN_LENGTH + params.COLUMN_WIDTH + 2 * dsW)
        * fcv,
        short_name="φVP",
        long_name="Design Punching Shear Capacity",
        format_func=lambda x: f"{x:.1f} kN",
    ),
)

calculator.add_calculation(
    "dvL",
    Calculation(
        func=lambda dsL, FtgDepth: np.maximum(0.9 * dsL, 0.72 * FtgDepth),
        short_name="dvL",
        long_name="Effective Shear Depth (L)",
        format_func=lambda x: f"{x*1000:.0f} mm",
    ),
)

calculator.add_calculation(
    "dvW",
    Calculation(
        func=lambda dsW, FtgDepth: np.maximum(0.9 * dsW, 0.72 * FtgDepth),
        short_name="dvW",
        long_name="Effective Shear Depth (W)",
        format_func=lambda x: f"{x*1000:.0f} mm",
    ),
)

calculator.add_calculation(
    "VOultL",
    Calculation(
        func=lambda FtgLength, dvW, BPult, SWt, FtgWidth: 0.5
        * (-params.COLUMN_LENGTH - 2 * dvW + FtgLength)
        * (BPult - (9 * SWt) / (10 * FtgLength * FtgWidth)),
        short_name="VO*L",
        long_name="Ultimate Shear Force (L)",
        format_func=lambda x: f"{x:.1f} kN",
    ),
)

calculator.add_calculation(
    "VOultW",
    Calculation(
        func=lambda FtgWidth, dvW, BPult, SWt, FtgLength: 0.5
        * (BPult - (9 * SWt) / (10 * FtgLength * FtgWidth))
        * (-params.COLUMN_WIDTH - 2 * dvW + FtgWidth),
        short_name="VO*W",
        long_name="Ultimate Shear Force (W)",
        format_func=lambda x: f"{x:.1f} kN",
    ),
)

calculator.add_calculation(
    "MOultL",
    Calculation(
        func=lambda FtgLength, dvW, SWt, BPult, FtgWidth: (
            (params.COLUMN_LENGTH + 2 * dvW - FtgLength) ** 2
            * (-9 * SWt + 10 * BPult * FtgLength * FtgWidth)
        )
        / (80 * FtgLength * FtgWidth),
        short_name="MO*L",
        long_name="Ultimate Moment (L)",
        format_func=lambda x: f"{x:.1f} kNm",
    ),
)

calculator.add_calculation(
    "MOultW",
    Calculation(
        func=lambda FtgWidth, dvW, SWt, BPult, FtgLength: (
            (params.COLUMN_WIDTH + 2 * dvW - FtgWidth) ** 2
            * (-9 * SWt + 10 * BPult * FtgLength * FtgWidth)
        )
        / (80 * FtgLength * FtgWidth),
        short_name="MO*W",
        long_name="Ultimate Moment (W)",
        format_func=lambda x: f"{x:.1f} kNm",
    ),
)

calculator.add_calculation(
    "ex1L",
    Calculation(
        func=lambda VOultL, dvL, MOultL, AstL: np.minimum(
            (np.maximum(VOultL * dvL, MOultL) / dvL + VOultL)
            / (2 * 200000 * AstL)
            * 1000,
            3 / 1000,
        ),
        short_name="ex1L",
        long_name="Eccentricity (L)",
        format_func=lambda x: f"{x*1000:.3f} mm",
    ),
)

calculator.add_calculation(
    "ex1W",
    Calculation(
        func=lambda VOultW, dvW, MOultW, AstW: np.minimum(
            (np.maximum(VOultW * dvW, MOultW) / dvW + VOultW)
            / (2 * 200000 * AstW)
            * 1000,
            3 / 1000,
        ),
        short_name="ex1W",
        long_name="Eccentricity (W)",
        format_func=lambda x: f"{x*1000:.3f} mm",
    ),
)

calculator.add_calculation(
    "kvL",
    Calculation(
        func=lambda dvL, ex1L: 13 / (25 * (1 + dvL) * (1 + 1500 * ex1L)),
        short_name="kvL",
        long_name="Shear Factor (L)",
        format_func=lambda x: f"{x:.3f}",
    ),
)

calculator.add_calculation(
    "kvW",
    Calculation(
        func=lambda dvW, ex1W: 13 / (25 * (1 + dvW) * (1 + 1500 * ex1W)),
        short_name="kvW",
        long_name="Shear Factor (W)",
        format_func=lambda x: f"{x:.3f}",
    ),
)

calculator.add_calculation(
    "AngleL",
    Calculation(
        func=lambda ex1L: 29 + 7000 * ex1L,
        short_name="AngleL",
        long_name="Angle (L)",
        format_func=lambda x: f"{x:.2f}°",
    ),
)

calculator.add_calculation(
    "AngleW",
    Calculation(
        func=lambda ex1W: 29 + 7000 * ex1W,
        short_name="AngleW",
        long_name="Angle (W)",
        format_func=lambda x: f"{x:.2f}°",
    ),
)

calculator.add_calculation(
    "ks",
    Calculation(
        func=lambda FtgDepth: np.maximum(1 / 2, (10 / 7) * (1 - FtgDepth)),
        short_name="ks",
        long_name="Size Factor",
        format_func=lambda x: f"{x:.3f}",
    ),
)

calculator.add_calculation(
    "fVucL",
    Calculation(
        func=lambda dvL, fc, ks, kvL: 700 * dvL * np.sqrt(fc) * ks * kvL,
        short_name="φVucL",
        long_name="Design Shear Strength (L)",
        format_func=lambda x: f"{x:.1f} kN",
    ),
)

calculator.add_calculation(
    "fVucW",
    Calculation(
        func=lambda dvW, fc, ks, kvW: 700 * dvW * np.sqrt(fc) * ks * kvW,
        short_name="φVucW",
        long_name="Design Shear Strength (W)",
        format_func=lambda x: f"{x:.1f} kN",
    ),
)

calculator.add_calculation(
    "Bpr",
    Calculation(
        func=lambda BPmax: BPmax / params.BEARINGPRESSURE,
        short_name="BPr",
        long_name="Bearing Pressure Ratio",
        format_func=lambda x: f"{x:.3f}",
    ),
)

calculator.add_calculation(
    "Mur",
    Calculation(
        func=lambda MultL, fMuoL, MultW, fMuoW: np.maximum(
            MultL / fMuoL, MultW / fMuoW
        ),
        short_name="Mur",
        long_name="Moment Utilization Ratio",
        format_func=lambda x: f"{x:.3f}",
    ),
)

calculator.add_calculation(
    "VPr",
    Calculation(
        func=lambda VPult, fVP: VPult / fVP,
        short_name="VPr",
        long_name="Punching Shear Ratio",
        format_func=lambda x: f"{x:.3f}",
    ),
)

calculator.add_calculation(
    "VOr",
    Calculation(
        func=lambda VOultL, fVucL, VOultW, fVucW: np.maximum(
            VOultL / fVucL, VOultW / fVucW
        ),
        short_name="VOr",
        long_name="One-way Shear Ratio",
        format_func=lambda x: f"{x:.3f}",
    ),
)

calculator.add_calculation(
    "Cost",
    Calculation(
        func=lambda AstW, AstL, FtgLength, FtgWidth, FtgDepth, fc: (
            AstW / 1000000 * FtgLength + AstL / 1000000 * FtgWidth
        )
        * 7850
        * 3.400
        + FtgLength * FtgWidth * FtgDepth * (130.866 * np.exp(fc * 0.0111) + 45 + 130)
        + 2 * FtgDepth * FtgLength * FtgWidth * 180,
        short_name="Cost",
        long_name="Total Cost",
        format_func=lambda x: f"${x:.2f}",
    ),
)

# Initialize FoundationSizer and perform calculations
foundation = FoundationSizer(params, calculator)
foundation.calculate_all_properties()

# Remove fails, sort by cost, and print results
# print(foundation.sizes)
# foundation.remove_fails()
# foundation.sort_by_cost()
foundation.print_foundation_details(0)
# foundation.print_array(10)
