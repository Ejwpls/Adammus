import math

def gen_ranges():
    FTG_CONC_STRENGTHS = [25, 32, 40, 50, 65]

    FTG_REO_SIZES = [0.012, 0.016, 0.02, 0.024, 0.028, 0.032, 0.036, 0.04]

    FTG_REO_CTS = [0.100, 0.120, 0.125, 0.140, 0.150, 0.160, 0.175, 0.180, 0.200, 0.220, 0.225, 0.240, 0.250, 0.260, 0.275, 0.280, 0.300]


    return (
        FTG_CONC_STRENGTHS,
        FTG_REO_SIZES,
        FTG_REO_CTS
    )

def gen_starting_values(
    COLUMN_LENGTH,
    COLUMN_WIDTH,
    BEARINGPRESSURE,
    LOAD_DEAD,
    LOAD_LIVE,
):
    load_sum = LOAD_DEAD + LOAD_LIVE
    column_diff = COLUMN_LENGTH - COLUMN_WIDTH
    common_term = math.sqrt(BEARINGPRESSURE * (BEARINGPRESSURE * column_diff**2 + 4 * load_sum))

    START_WIDTH = (COLUMN_WIDTH - COLUMN_LENGTH) / 2 + common_term / (2 * BEARINGPRESSURE)
    START_LENGTH = (COLUMN_LENGTH - COLUMN_WIDTH) / 2 + common_term / (2 * BEARINGPRESSURE)

    FTG_LENGTH_MIN = math.ceil(START_LENGTH / 0.05) * 0.05
    FTG_WIDTH_MIN = math.ceil(START_WIDTH / 0.05) * 0.05

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

    return (
        FTG_LENGTH_MIN,
        FTG_WIDTH_MIN,
        FTG_DPTH_MIN
    )