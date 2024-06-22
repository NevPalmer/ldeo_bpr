@call conda activate geod

@python C:\scripts\seafloor-geod\apg_read_ignore_timestamps.py ^
    --infile        D:/TAN2314_data/2023_TAN2314/raw/TAN2314_GNS22-PL_136305.apg ^
    --snapg         136305 ^
    --version       Seascan2018 ^
    --tempsmth      5001 ^
    --decimate      1 ^
    --clkstart      2022-10-13_02:57:30 ^
    --beginwndw     2023-03-01_00:00:00.0 ^
    --endwndw       2023-03-02_00:00:00.0 ^
    --bininterval   1d ^
    --gpssynctime   2023:276:21:31:30 ^
    --fmttime       d ^
    --plot          rd

    @REM --outfile       .\out\GNS22-PL.csv ^
    @REM --plot          n ^
    @REM | tee ./out/GNS22-PL.log


@REM @python C:\scripts\seafloor-geod\apg_read.py ^
@REM     --infile        D:/TAN2314_data/2023_TAN2314/raw/TAN2314_GNS22-PF2_140343.apg ^
@REM     --snapg         140343 ^
@REM     --version       Seascan2018 ^
@REM     --tempsmth      5001 ^
@REM     --decimate      1 ^
@REM     --clkstart      2023-01-06_13:08:30 ^
@REM     --beginwndw     2023-03-01_00:00:00.0 ^
@REM     --endwndw       2023-03-02_00:00:00.0 ^
@REM     --bininterval   1d ^
@REM     --gpssynctime   2023:277:23:56:00 ^
@REM     --fmttime       d ^
@REM     --plot          rd

    @REM --outfile       .\out\GNS22-PF2.csv ^
    @REM --plot          n ^
    @REM | tee ./out/GNS22-PF2.log

@REM A little noisy. ~ 500 Pa
@REM @python C:\scripts\seafloor-geod\apg_read.py ^
@REM     --infile        D:/TAN2314_data/2023_TAN2314/raw/TAN2314_GNS22-PH_140337.apg ^
@REM     --snapg         140337 ^
@REM     --version       Seascan2018 ^
@REM     --tempsmth      5001 ^
@REM     --decimate      1 ^
@REM     --clkstart      2023-01-05_21:02:00 ^
@REM     --beginwndw     2023-03-01_00:00:00.0 ^
@REM     --endwndw       2023-03-02_00:00:00.0 ^
@REM     --bininterval   1d ^
@REM     --gpssynctime   2023:278:15:31:30 ^
@REM     --fmttime       d ^
@REM     --plot          rd

    @REM --outfile       .\out\GNS22-PH.csv ^
    @REM --plot          n
    @REM | tee ./out/GNS22-PH.log

@call conda deactivate