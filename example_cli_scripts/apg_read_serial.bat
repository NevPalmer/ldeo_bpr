@call conda activate "C:\Users\nevillep\AppData\Local\miniforge3\envs\geod"

python "C:\scripts\ldeo_bpr\apg_read_serial.py" ^
    -s "136309" ^
    -v "CSAC2013" ^
    --serport COM3 ^
    --serbaud 9600

@call conda deactivate
