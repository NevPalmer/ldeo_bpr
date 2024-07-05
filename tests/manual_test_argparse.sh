#!/usr/bin/env bash
python -m tests.manual_test_argparse \
    --infile        "TAN1907_GNS18-01_136309.apg" \
    --snapg         136309    \
    --loggerversion CSAC2013 \
    --tempsmth      1 \
    --decimate      0 \
    --clkstart      2018-10-12_07:11:30 \
    --beginwndw     2019-04-29_03:00:00.0 \
    --endwndw       2019-04-29_04:00:00.0 \
    --bininterval   24H \
    --gpssynctime   2019-307_19:49:30 \
    --synctickcount 0x03CBB30B84 \
    --plot          rd \
    --fmttime       d \
    --network       XX \
    --station       G20PE \
    --outfile       output_GNS18-01.csv \
    --mseedpath     ./mseed/

# This is the first and only rollover of the data
    # --beginwndw     2019-04-29_03:00:00.0 \
    # --endwndw       2019-04-29_04:00:00.0 \

# This is the first rollover of the data
# python ~/scripts/ldeo_bpr/apg_read.py \
#     --infile        TAN1907_GNS18-07_140338.apg \
#     --snapg         140338 \
#     --loggerversion Seascan2018 \
#     --tempsmth      1 \
#     --decimate      0 \
#     --clkstart      2019-02-25_14:00:30 \
#     --beginwndw     2019-04-16_07:00:30.0 \
#     --endwndw       2019-04-16_08:00:30.0 \
#     --plot          r \
#     --fmttime       d \
#     --outfile      output_GNS18-07.csv

# This is the first rollover of the data
#     --beginwndw     2019-04-16_07:00:30.0 \
#     --endwndw       2019-04-16_08:00:30.0 \

# This is the last rollover of the data
#     --beginwndw     2019-09-12_10:00:00.0 \
#     --endwndw       2019-09-12_11:00:00.0 \

# This is the end of the data
#     --beginwndw     2019-10-30_02:52:20.0 \
#     --endwndw       2019-10-30_02:52:40.0 \
