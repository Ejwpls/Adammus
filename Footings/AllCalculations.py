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
            + FtgLength
            * FtgWidth
            * FtgDepth
            * (130.866 * np.exp(fc * 0.0111) + 45 + 130)
            + 2 * FtgDepth * FtgLength * FtgWidth * 180,
            short_name="Cost",
            long_name="Total Cost",
            format_func=lambda x: f"${x:.2f}",
        ),
    )