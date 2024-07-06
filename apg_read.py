"""A script for extracting raw data from LDEO type APG data loggers."""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import obspy
import scipy.signal as sig
import scipy.stats as stat

import ldeo_bpr as bpr
from ldeo_bpr import dt64_utils


def main():
    """The first function run when this script is run directly."""
    # Dictionary of flags for turning on/off trouble shooting outputs.
    # Assign False or '' to disable.
    trbl_sht: dict[str, bool] = {
        # Save raw data as a text file of binary 1s & 0s. (Very slow)
        "binary_out": False,
        # Save raw data as a text file of hexadecimal values. (Very slow)
        "hex_out": False,
        # Save raw records to file as integers without tick rollover removed.
        "raw_rlovr": False,
        # Save raw records to file as integers with tick rollover removed.
        "raw_no_rlovr": False,
        # Save time sync records to file as integers with tick rollover remved.
        "raw_sync": False,
        # Write a summary file showing syncronised tic counts and exit.
        "tic_sync": False,
    }

    # Retrieve CLI parameters.
    args = bpr.parse_arguments()

    print("=" * 80)
    print(f"STATS OF RAW FILE: {args.infile}")
    print(f"Time of first sample = {args.clkstart}")

    # Read Paros transducer coefficients from .ini file.
    paros = bpr.Paros.from_file(filename=args.apgini, paros_sn=args.snapg)

    # Read APG logger configuration parameters from .ini file.
    logger = bpr.Logger.from_file(
        filename=args.loggerini,
        logger_version=args.loggerversion,
    )

    # Calculate duration and end time of raw file.
    try:
        fsize = args.infile.stat().st_size  # file size in bytes
    except FileNotFoundError:
        sys.exit(f'Raw APG file "{args.infile}" does not exist.')
    print(f"Filesize = {fsize} bytes")
    nrecs = int((fsize - logger.head_len) / logger.rec_len) - 1
    print(f"Number of records = {nrecs:d}")
    file_duration_ms = nrecs * logger.record_epoch
    file_duration_secs = file_duration_ms / 1000
    file_duration_days = file_duration_secs / 3600 / 24
    print(f"File duration = {file_duration_days} days")
    clk_end_dt = args.clkstart + np.timedelta64(file_duration_ms, "ms")
    print(f"Time of last sample = {clk_end_dt}")
    print("=" * 80)
    print("STATS OF DATA WINDOW TO BE EXTRACTED:")

    if args.beginwndw is None:
        args.beginwndw = args.clkstart
    if args.endwndw is None:
        args.endwndw = clk_end_dt
    print(f"Window beginning: {args.beginwndw}")
    wndw_len = args.endwndw - args.beginwndw
    wndw_len_ms = dt64_utils.delta64_to_ms(wndw_len)
    wndw_len_days = wndw_len_ms / (24 * 3600000)
    print(f"Window length (days): {wndw_len_days}")

    # Clock drift
    if args.gpssynctime is not None:
        drift = clockdrift(
            args.infile,
            logger,
            args.clkstart,
            args.gpssynctime,
            args.synctickcount,
            trbl_sht,
        )
        print(
            f"Clock drift at end of recording (millisecs): {drift}\n"
            f"   (logged time - actual GPS time)"
        )
    else:
        drift = 0
        args.gpssynctime = clk_end_dt

    print(
        f"NOTE: All times given above are Nominal. \n"
        f'"Nominal" times are calculated by counting the number '
        f"of sample epochs multiplied by the sample period "
        f"({logger.sample_epoch} secs).\n"
        f'"Actual" times are the actual recorded tick count values '
        f"before any adjustment for clock drift.\n"
        f"For some APG loggers, Nominal and Actual times correspond "
        f"precisely, some do not. Generally CSAC loggers do and Seascan "
        f"loggers do not.\n"
        f"If a difference is noted below, the nominal epochs are not "
        f"precise and the tick count values have been used for precise "
        f"timing.\n"
        f"All subsequent output times are base on Actual times, plus "
        f"correction for clock drift, where this is provided.\n"
    )

    print("=" * 80)
    print("STATS FOR EACH TIME BIN:")

    # Loop to extract data in time bins
    if args.outfile:
        # Create empty file, overwrite if exists.
        open(args.outfile, "w", encoding="utf8").close()
    wndw_beg_unix_ms = dt64_utils.dt64_to_ms(args.beginwndw)
    wndw_end_unix_ms = dt64_utils.dt64_to_ms(args.endwndw)
    bin_beg_ms = wndw_beg_unix_ms
    if args.bininterval == 0:
        bin_end_ms = wndw_end_unix_ms
    else:
        bin_end_ms = bin_beg_ms - (bin_beg_ms % args.bininterval) + args.bininterval
    while bin_beg_ms < wndw_end_unix_ms:
        if bin_end_ms > wndw_end_unix_ms:
            bin_end_ms = wndw_end_unix_ms
        bin_begin_dt = dt64_utils.ms_to_dt64(bin_beg_ms)
        bin_end_dt = dt64_utils.ms_to_dt64(bin_end_ms)
        print("-" * 80)
        print(f"Processing time bin: start: {bin_begin_dt} | end: {bin_end_dt}")

        generate_results(
            logger,
            paros,
            args.network,
            args.station,
            args.clkstart,
            bin_begin_dt,
            bin_end_dt,
            file_duration_ms,
            args.gpssynctime,
            drift,
            args.infile,
            args.decimate,
            args.tempsmth,
            args.noisefilt,
            args.fmttime,
            args.plot,
            args.outfile,
            args.mseedpath,
            trbl_sht,
        )

        bin_beg_ms = bin_end_ms
        bin_end_ms = bin_beg_ms + args.bininterval


###############################################################################
def generate_results(
    logger: bpr.Logger,
    paros: bpr.Paros,
    nwk_name,
    stn_name,
    clk_start_dt,
    bin_begin_dt,
    bin_end_dt,
    file_duration_ms,
    gpssync_dt,
    drift,
    apg_filename: Path,
    decmt_intvl,
    tmptr_smth_fctr,
    noisefilt,
    time_format,
    plot_flags,
    out_filename: Path,
    mseed_path: Path,
    trbl_sht: dict[str, bool],
):
    """This is the primary function used to extract and output results."""
    # Calculate window for extracting data.
    if noisefilt:
        bin_padding = 1_000_000  # milliseconds
    else:
        bin_padding = (tmptr_smth_fctr - 1) * logger.smpls_per_rec * 10
        # bin_padding = 0
        # bin_padding = 600000  # milliseconds

    bin_begin_ms = dt64_utils.delta64_to_ms(bin_begin_dt - clk_start_dt)
    bin_len_ms = dt64_utils.delta64_to_ms(bin_end_dt - bin_begin_dt)
    bin_end_ms = bin_begin_ms + bin_len_ms

    # Make sure the requested times don't fall outside available records.
    if (bin_begin_ms - bin_padding) < 0:
        padded_bin_begin_ms = 0
    else:
        padded_bin_begin_ms = bin_begin_ms - bin_padding

    padded_bin_len_ms = bin_end_ms - padded_bin_begin_ms + bin_padding
    avail_bin_len_ms = file_duration_ms - bin_begin_ms
    if (avail_bin_len_ms - padded_bin_len_ms) <= 0:
        padded_bin_len_ms = padded_bin_len_ms - bin_padding

    rec_begin = int(padded_bin_begin_ms / logger.record_epoch)
    nrecs_want = int(padded_bin_len_ms / logger.record_epoch) + 1

    # Extract records from APG file for the time window specified.
    print("Extracting raw records:\n", end="")

    records = extractrecords(
        apg_filename, logger, nrecs_want, rec_begin, trbl_sht, clk_start_dt
    )

    # Save raw records to file as integers with tick rollover removed.
    if trbl_sht["raw_no_rlovr"]:
        np.savetxt(
            "raw_records_no-rollover.txt", records, fmt="%d", header="", comments=""
        )

    # Assign names to column numbers of raw data array.
    # Note that raw data array columns are reverse order to raw binary.
    last_field = len(logger.rec_fmt) - 1
    pcore_col = last_field - logger.fmt_field["pcore"]
    pn_col = last_field - logger.fmt_field["pn"]
    tptr_col = last_field - logger.fmt_field["tptr"]

    # Create an array for each raw observable (pressure, temperature, ticks_ms)
    if pn_col > pcore_col:
        press_raw = records[:, pcore_col : pn_col + 1]
    else:
        press_raw = records[:, pcore_col : pn_col - 1 : -1]
    press_raw = np.cumsum(press_raw, axis=1)
    press_raw = press_raw.reshape(nrecs_want * logger.smpls_per_rec)
    temp_raw = records[:, tptr_col]
    ticks_ms = records[:, last_field - logger.fmt_field["tic"]]

    actual_end_tick = ticks_ms[-1]
    actual_begin_tick = ticks_ms[0]
    nominal_begin_tick = rec_begin * logger.record_epoch
    nominal_end_tick = (rec_begin + nrecs_want - 1) * logger.record_epoch
    # print(f"\nActual beginning tick:  {actual_begin_tick}")
    # print(f"Nominal beginning tick: {nominal_begin_tick}")
    # print(f"Actual end tick:        {actual_end_tick}")
    # print(f"Nominal end tick:       {nominal_end_tick}\n")

    nom_ticks_t = np.linspace(
        nominal_begin_tick, nominal_end_tick, num=nrecs_want, endpoint=True
    ).astype("int64")
    nom_ticks_p = np.linspace(
        nominal_begin_tick,
        nominal_end_tick + logger.record_epoch,
        num=nrecs_want * logger.smpls_per_rec,
        endpoint=False,
    ).astype("int64")

    # Write a summary file showing syncronised tic counts and exit.
    if trbl_sht["tic_sync"]:
        time_diff = ticks_ms - nom_ticks_t
        ticks_ms = np.column_stack((ticks_ms, nom_ticks_t, time_diff))
        summary_ticks_ms = [ticks_ms[0]]
        for i in range(1, len(time_diff)):
            if time_diff[i] != time_diff[i - 1]:
                summary_ticks_ms.append(ticks_ms[i])
        np.savetxt(
            "summary_ticks_ms.txt",
            summary_ticks_ms,
            fmt="%d",
            header="actual,nominal,difference",
            comments="",
        )
        np.savetxt(
            "ticks_ms.txt",
            ticks_ms,
            fmt="%d",
            header="actual,nominal,difference",
            comments="",
        )
        sys.exit()

    # If the nominal tick count and actual recorded tick count are not
    # precisely aligned then use linear interpolation to generate a precisely
    # periodic record of raw temperature and pressure values.
    if actual_begin_tick == nominal_begin_tick and actual_end_tick == nominal_end_tick:
        millisecs_t = nom_ticks_t
        millisecs_p = nom_ticks_p
    else:
        beg_diff = nominal_begin_tick - actual_begin_tick
        if beg_diff < 0:
            dirn = "behind"
        else:
            dirn = "ahead"
        print(
            f"Nominal time at start of window is {abs(beg_diff)/1000} "
            f"seconds {dirn} Actual recorded time ticks.",
            flush=True,
        )
        end_diff = nominal_end_tick - actual_end_tick
        if end_diff < 0:
            dirn = "behind"
        else:
            dirn = "ahead"
        print(
            f"Nominal time at end of window is {abs(end_diff)/1000} "
            f"seconds {dirn} Actual recorded time ticks.",
            flush=True,
        )

        # Determine first tick count to achieve fixed period epochs.
        final_begin_tick = actual_begin_tick + (
            logger.record_epoch - actual_begin_tick % logger.record_epoch
        )
        # Determine final tick count to achieve fixed period epochs.
        final_end_tick = actual_end_tick - (actual_end_tick % logger.record_epoch)
        epoch_count = int((final_end_tick - final_begin_tick) / logger.record_epoch) + 1
        millisecs_t = np.linspace(
            final_begin_tick, final_end_tick, num=epoch_count, endpoint=True
        ).astype("int64")
        millisecs_p = np.linspace(
            final_begin_tick,
            final_end_tick + logger.record_epoch,
            num=epoch_count * logger.smpls_per_rec,
            endpoint=False,
        ).astype("int64")
        ticks_ext = np.append(ticks_ms, ticks_ms[-1] + logger.record_epoch)
        nom_ticks_t = np.append(nom_ticks_t, nom_ticks_t[-1] + logger.record_epoch)

        ticks_p = np.interp(nom_ticks_p, nom_ticks_t, ticks_ext)

        # Interpolate to generate fixed period observation epochs.
        temp_raw = np.interp(millisecs_t, ticks_ms, temp_raw)
        press_raw = np.interp(millisecs_p, ticks_p, press_raw)

    # Apply clock drift to time values
    # Clock drift is fixed at the mid-point of the period of extracted data
    # and any change is assumed to be insignificant over that period.
    millisecs_logged = dt64_utils.delta64_to_ms(gpssync_dt - clk_start_dt)
    drift_beg = drift * (millisecs_t[0] / millisecs_logged)
    drift_end = drift * (millisecs_t[-1] / millisecs_logged)
    drift_applied = (drift_beg + drift_end) / 2
    # Round the drift to be applied to the  nearest whole sample epoch.
    drift_applied = int(
        logger.sample_epoch * round(drift_applied / logger.sample_epoch, 0)
    )
    print(
        f"Clock drift applied to the extracted block:  {drift_applied} ms.",
        flush=True,
    )
    # Apply clock drift to time records.
    millisecs_p = millisecs_p - drift_applied
    millisecs_t = millisecs_t - drift_applied

    # Temperature period (usec)
    tmptr_period_usec = (temp_raw / (logger.tp_fctr) + logger.tp_cnst) / (
        logger.clock_freq
    )
    # Uncomment one of the lines below to ignore logged temperature values and
    # assume a fixed value instead.
    # tmptr_period_usec.fill(5.8224) #Fixed temp for Paros 140344
    # tmptr_period_usec.fill(5.7900) #Fixed temp for Paros 140346
    # tmptr_period_usec.fill(5.7875) #Fixed temp of +3.06°C for Paros 140339
    # tmptr_period_usec.fill(5.7830) #Fixed temp of +20.65°C for Paros 140339
    # tmptr_period_usec.fill(5.7530) #Fixed temp of +2.27°C for Paros 140338
    # tmptr_period_usec.fill(5.8475) #Fixed temp of +2.55°C for Paros 136309
    # tmptr_period_usec.fill(5.8430) #Fixed temp of +20.08°C for Paros 136309

    # Pressure period (usec)
    presr_period_usec = (press_raw / (logger.pp_fctr) + logger.pp_cnst) / (
        logger.clock_freq
    )

    tmptr_period_usec_raw = tmptr_period_usec
    coef_uv_raw = tmptr_period_usec_raw - paros.u[0]
    temperature_raw = np.polyval(paros.y, coef_uv_raw)

    if noisefilt:
        # Temperature spike/noise removal
        print("Temperature spike/noise removal.")
        # Very course first pass removal based on impossible extremes.
        mask = np.where((temperature_raw < 0) | (temperature_raw > 30))
        temperature_del = np.delete(temperature_raw, mask)
        tmptr_period_usec_del = np.delete(tmptr_period_usec, mask)
        millisecs_t_del = np.delete(millisecs_t, mask)

        # More refined noise removal based on binned median values
        tmptr_period_usec_filt, millisecs_filt = remove_noise_meddiff(
            raw_data=tmptr_period_usec_del,
            mask_data=temperature_del,
            millisecs=millisecs_t_del,
            millisecs_all=millisecs_t,
            bin_size=1_000_000,
            tolerance=0.02,
        )
        coef_uv_filt = tmptr_period_usec_filt - paros.u[0]
        temperature_filt = np.polyval(paros.y, coef_uv_filt)
        tmptr_period_usec_filt, millisecs_filt = remove_noise_meddiff(
            raw_data=tmptr_period_usec_filt,
            mask_data=temperature_filt,
            millisecs=millisecs_filt,
            millisecs_all=millisecs_t,
            bin_size=1_000,
            tolerance=0.005,
        )
        tmptr_period_usec = np.interp(
            millisecs_t, millisecs_filt, tmptr_period_usec_filt
        )

    # Apply smoothing filter to temperature before calculating pressure.
    # This eliminates significant noise from the pressure values.
    if tmptr_smth_fctr >= 5:
        print("Applying temperature smoothing filter.", flush=True)
        tmptr_period_usec = sig.savgol_filter(
            tmptr_period_usec, tmptr_smth_fctr, 3, axis=0, mode="mirror"
        )

    # Calculate temperature array
    print("Calculating temperatures and pressures.", flush=True)
    coef_uv = tmptr_period_usec - paros.u[0]
    # Upsample temperatures to match frequency of pressure samples by
    #  linear interpolation.
    coef_uv_expnd = np.interp(millisecs_p, millisecs_t, coef_uv)
    temperature_upsmpld = np.polyval(paros.y, coef_uv_expnd)

    # Calculate pressure array
    coef_cv = np.polyval(paros.c, coef_uv_expnd)
    coef_dv = np.polyval(paros.d, coef_uv_expnd)
    coef_t0 = np.polyval(paros.t, coef_uv_expnd)

    facts = 1 - (coef_t0**2) / (presr_period_usec**2)
    pressure = coef_cv * facts * (1 - coef_dv * facts)  # pressure in PSIA
    pressure = pressure * bpr.PRESS_CONV_FCTR  # Convert pressure units
    pressure_raw = pressure

    if noisefilt:
        # Pressure spike/noise removal
        print("Pressure spike/noise removal.")

        # Very course first pass removal based on likely spread from overall median
        difference = np.abs(pressure - np.median(pressure))
        mask = np.where(difference > 20_000)
        # Very course first pass removal based on impossible extremes.
        # mask = np.where((pressure > 50_000_000) | (pressure < 0))
        pressure_del = np.delete(pressure, mask)
        millisecs_del = np.delete(millisecs_p, mask)

        # Refined noise removal based on binned median values
        pressure_filt, millisecs_filt = remove_noise_meddiff(
            raw_data=pressure_del,
            mask_data=pressure_del,
            millisecs=millisecs_del,
            millisecs_all=millisecs_p,
            bin_size=10_000,
            tolerance=50,
        )

        # Uncomment line below to plot first iteration pressure instead of raw pressure.
        # pressure_raw = pressure

        # Second iteration noise removal with tighter tollerances
        pressure_filt, millisecs_filt = remove_noise_meddiff(
            raw_data=pressure_filt,
            mask_data=pressure_filt,
            millisecs=millisecs_filt,
            millisecs_all=millisecs_p,
            bin_size=500,
            tolerance=15,
        )
        pressure = np.interp(millisecs_p, millisecs_filt, pressure_filt)

    # Decimate results
    # To produce sensible decimation results when ftype='iir',
    #    ensure n=5 and iterate using downsampling factor (q) <= 13.
    pressure_dcmtd = []
    if decmt_intvl > 0:
        print(
            f"Decimating results to {decmt_intvl} second epochs by "
            f"iteration, using factors;",
            end="",
            flush=True,
        )
        pressure_dcmtd = pressure
        temperature_dcmtd = temperature_upsmpld
        millisecs_dcmtd = millisecs_p
        # Ensure first record in data arrays starts at a whole multiple of the
        # decimation factor.
        actual_first_tick = millisecs_dcmtd[1]
        actual_last_tick = millisecs_dcmtd[-1]
        intvl_ms = decmt_intvl * 1000
        dcmtd_first_tick = actual_first_tick + intvl_ms - actual_first_tick % intvl_ms
        dcmtd_last_tick = actual_last_tick - (actual_last_tick % intvl_ms)
        mask = np.logical_and(
            millisecs_dcmtd >= dcmtd_first_tick, millisecs_dcmtd <= dcmtd_last_tick
        )
        millisecs_dcmtd = millisecs_dcmtd[mask]
        pressure_dcmtd = pressure_dcmtd[mask]
        temperature_dcmtd = temperature_dcmtd[mask]

        # First decimate to whole seconds
        sample_freq = int(1000 / logger.sample_epoch)
        decmt_intvl_pre = sample_freq
        first = True
        while decmt_intvl_pre > 1:
            if decmt_intvl_pre % 5 == 0:
                decmt_fctr = 5
            else:
                decmt_fctr = decmt_intvl_pre
            if not first:
                print(" :", end="")
            first = False
            print(f" {decmt_fctr}", end="", flush=True)
            pressure_dcmtd = sig.decimate(pressure_dcmtd, decmt_fctr, n=5, ftype="iir")
            temperature_dcmtd = sig.decimate(
                temperature_dcmtd, decmt_fctr, n=5, ftype="iir"
            )
            millisecs_dcmtd = millisecs_dcmtd[::decmt_fctr]
            decmt_intvl_pre = decmt_intvl_pre // decmt_fctr
        # Now decimate to number of whole seconds requested
        while decmt_intvl > 1:
            if decmt_intvl % 10 == 0:
                decmt_fctr = 10
            elif decmt_intvl % 6 == 0:
                decmt_fctr = 6
            elif decmt_intvl % 5 == 0:
                decmt_fctr = 5
            elif decmt_intvl < 10 and decmt_intvl % 1 == 0:
                decmt_fctr = decmt_intvl
            else:
                sys.exit(
                    "\nDecimation failed! The interval specified must be "
                    "a single digit number of minutes or seconds, or be "
                    "divisible by 5 or 10."
                )
            print(f": {decmt_fctr}", end="", flush=True)
            pressure_dcmtd = sig.decimate(pressure_dcmtd, decmt_fctr, n=5, ftype="iir")
            temperature_dcmtd = sig.decimate(
                temperature_dcmtd, decmt_fctr, n=5, ftype="iir"
            )
            millisecs_dcmtd = millisecs_dcmtd[::decmt_fctr]
            decmt_intvl = decmt_intvl // decmt_fctr
        print()

    if len(pressure_dcmtd) != 0:
        time = millisecs_dcmtd
        pressure_out = pressure_dcmtd
        temperature_out = temperature_dcmtd
    else:
        time = millisecs_p
        pressure_out = pressure
        temperature_out = temperature_upsmpld

    # Trim results to the originally specified time bin.
    bin_end_ms = bin_begin_ms + bin_len_ms
    mask = np.logical_and(time >= bin_begin_ms, time < bin_end_ms)
    time = time[mask]
    pressure_out = pressure_out[mask]
    temperature_out = temperature_out[mask]

    # Convert timestamp from seconds to datetime if specified.
    if time_format == "d":
        print("Calculating a timestamp array from the milliseconds array.", flush=True)
        time = (clk_start_dt + dt64_utils.ms_to_delta64(time)).astype("datetime64[ms]")
        xlabel = "Date-Time"
    else:
        time = time / 1000
        xlabel = "Seconds"

    # Save results to CSV file
    if out_filename:
        print(f'Appending results to file "{out_filename}".')
        if time_format == "d":
            time = time.astype("datetime64[ms]")
        content = zip(time, pressure_out, temperature_out)
        with open(out_filename, "a", newline="", encoding="utf8") as csvfile:
            for row in content:
                csvfile.write(f"{row[0]},{row[1]:14.3f},{row[2]:10.6f}\n")

    # Save results to MiniSEED file
    if mseed_path:
        # Trim results to the originally specified time bin.
        bin_end_ms = bin_begin_ms + bin_len_ms
        mask = np.logical_and(millisecs_p >= bin_begin_ms, millisecs_p < bin_end_ms)
        pressure_mseed = pressure[mask]
        mask = np.logical_and(millisecs_t >= bin_begin_ms, millisecs_t < bin_end_ms)
        temperature_mseed = temperature_raw[mask]
        bin_begin = dt64_utils.dt64_to_pydt(bin_begin_dt)
        stats = {
            "network": nwk_name,
            "station": stn_name,
            "location": "",
            "starttime": bin_begin,
            "mseed": {"dataquality": "R"},
        }

        Path(mseed_path).mkdir(parents=True, exist_ok=True)

        dt_text = bin_begin.strftime("%Y-%m-%dT%H-%M-%Sz")

        stats["channel"] = bpr.P_CHNL_CODE
        stats["sampling_rate"] = 1000 / logger.sample_epoch
        trace_p = obspy.Trace(data=pressure_mseed, header=stats)
        stream_p = obspy.Stream(traces=[trace_p])
        mseed_filename = f"{nwk_name}_{stn_name}_{stats['channel']}_{dt_text}.mseed"
        mseed_filename = mseed_path / mseed_filename
        print(f'Writing pressure data to MiniSEED file "{mseed_filename}".')
        stream_p.write(mseed_filename, format="MSEED")

        stats["channel"] = bpr.T_CHNL_CODE
        stats["sampling_rate"] = 1000 / logger.record_epoch
        trace_t = obspy.Trace(data=temperature_mseed, header=stats)
        stream_t = obspy.Stream(traces=[trace_t])
        mseed_filename = f"{nwk_name}_{stn_name}_{stats['channel']}_{dt_text}.mseed"
        mseed_filename = mseed_path / mseed_filename
        print(f'Writing temperature data to MiniSEED file "{mseed_filename}".')
        stream_t.write(mseed_filename, format="MSEED")

        # Uncomment line(s) below to see MiniSEED plot output.
        # stream_t.plot()
        # stream_p.plot()

    # Generate and output a plot
    if plot_flags["format"] != "n":
        # Set min-max values for plot Y-axes
        print("Generating plot.", flush=True)
        p_min = np.min(pressure_out)
        p_max = np.max(pressure_out)
        p_range = p_max - p_min
        if p_range == 0:
            intvl = 10
        else:
            int_log = int(np.log10(p_range))
            if 10**int_log / p_range >= 0.5:
                intvl = 10**int_log / 10
            elif 10**int_log / p_range >= 0.2:
                intvl = 10**int_log / 5
            elif 10**int_log / p_range >= 0.1:
                intvl = 10**int_log / 2
            else:
                intvl = 10**int_log
        p_min = p_min - p_min % intvl - intvl
        p_max = p_max - p_max % intvl + 2 * intvl

        t_min = np.min(temperature_out)
        t_max = np.max(temperature_out)
        t_range = t_max - t_min
        if t_range == 0:
            intvl = 0.1
        else:
            int_log = int(np.log10(t_range))
            if 10**int_log / t_range >= 0.5:
                intvl = 10**int_log / 10
            elif 10**int_log / t_range >= 0.2:
                intvl = 10**int_log / 5
            elif 10**int_log / t_range >= 0.1:
                intvl = 10**int_log / 2
            else:
                intvl = 10**int_log
        t_min = t_min - t_min % intvl - intvl
        t_max = t_max - t_max % intvl + 2 * intvl

        # Plot Results
        fig, ax2 = plt.subplots(figsize=(15, 9))
        plt.ylim(t_min, t_max)
        ax1 = ax2.twinx()
        plt.ylim(p_min, p_max)

        # Plot raw pressure values if requested
        if plot_flags["format"] == "r":
            if time_format == "d":
                time_p = clk_start_dt + dt64_utils.ms_to_delta64(millisecs_p)
            else:
                time_p = millisecs_p / 1000
            ax1.plot(
                time_p,
                pressure_raw,
                color="pink",
                marker=".",
                markersize=1.0,
                linestyle="",
            )

        # Plot raw temperature values if requested
        if plot_flags["format"] == "r":
            if time_format == "d":
                time_t = clk_start_dt + dt64_utils.ms_to_delta64(millisecs_t)
            else:
                time_t = millisecs_t / 1000
            ax2.plot(
                time_t,
                temperature_raw,
                color="lightblue",
                marker=".",
                markersize=1.0,
                linestyle="",
            )

        # Plot final pressure values
        color = "red"
        ax1.set_ylabel("Pressure (Pa)", color=color)
        ax1.tick_params(axis="y", labelcolor=color)
        ax1.grid(axis="x")
        ax1.plot(
            time,
            pressure_out,
            color=color,
            marker=".",
            markersize=1.0,
            linestyle="solid",
            linewidth=0.5,
        )
        ax1.set_xlabel(xlabel)

        # Plot final temperature values
        color = "blue"
        ax2.set_ylabel("Temperature (°C)", color=color)
        ax2.tick_params(axis="y", labelcolor=color)
        ax2.plot(
            time,
            temperature_out,
            color=color,
            marker=".",
            markersize=1.0,
            linestyle="solid",
            linewidth=0.5,
        )

        fig.suptitle(f"{apg_filename}\nBegin: {bin_begin_dt}  -  " f"End: {bin_end_dt}")
        plt.tight_layout()
        if time_format == "d":
            # Rotates and aligns the X-axis labels.
            plt.gcf().autofmt_xdate(bottom=0.2, rotation=30, ha="right")
        if plot_flags["output"] in ["s", "b"]:
            basename = apg_filename.stem
            fig.savefig(f"./{basename}.png", dpi=200)
        if plot_flags["output"] in ["d", "b"]:
            plt.show()
        plt.close(fig)

    return


###########################################################################
def extractrecords(
    apg_filename: Path,
    logger: bpr.Logger,
    nrecs_want: int,
    rec_begin: int,
    trbl_sht: dict[str, bool],
    clk_start_dt: np.datetime64,
):
    """Extracts binary records from a raw APG data logger file."""
    if trbl_sht["binary_out"]:
        binary_filename = "raw_binary.txt"
        # Create empty file, overwrite if exists.
        open(binary_filename, "w", encoding="utf8").close()
    if trbl_sht["hex_out"]:
        hex_filename = "raw_hex.txt"
        # Create empty file, overwrite if exists.
        open(hex_filename, "w", encoding="utf8").close()

    with open(apg_filename, "rb") as apgfile:
        begin_byte = logger.head_len + rec_begin * logger.rec_len
        apgfile.seek(begin_byte, 0)
        records = []
        for i in range(0, nrecs_want):
            binary_record = apgfile.read(logger.rec_len)

            # Print record as a string of Binary values to file.
            if trbl_sht["binary_out"]:
                bin_str = ""
                for ch in binary_record:
                    bin_str += f"{ch:08b}"
                cum_rec_fmt = np.cumsum(list(map(abs, logger.rec_fmt)))
                bin_str_dlmtd = ""
                for count, char in enumerate(bin_str):
                    if count in cum_rec_fmt:
                        bin_str_dlmtd = bin_str_dlmtd + " "
                    bin_str_dlmtd = bin_str_dlmtd + char
                with open(binary_filename, "a", encoding="utf8") as binfile:
                    binfile.write(f"{bin_str_dlmtd}\n")

            # Write record as a string of Hexadecimal values to file.
            if trbl_sht["hex_out"]:
                hex_str = ""
                for ch in binary_record:
                    hex_str += hex(ch) + " "
                with open(hex_filename, "a", encoding="utf8") as hexfile:
                    hexfile.write(f"{hex_str}\n")

            # Split record into array of ints defined as groups of bits
            # by logger['rec_fmt'] .
            record_int = int.from_bytes(binary_record, byteorder="big", signed=False)
            record = []
            for signed_bit_len in reversed(logger.rec_fmt):
                bit_len = int(abs(signed_bit_len))
                # Read right most bit_len bits
                field = record_int & (2**bit_len - 1)
                # Check for sign bit and convert as a 2s-compliment negative.
                if (signed_bit_len < 0) and (field & (2 ** (bit_len - 1))):
                    full_bit_len = bit_len + (bit_len % 8)
                    field = field | ((2**full_bit_len - 1) ^ (2**bit_len - 1))
                    field = field.to_bytes(
                        full_bit_len // 8, byteorder="big", signed=False
                    )
                    field = int.from_bytes(field, byteorder="big", signed=True)
                record.append(field)
                # Shift to right bit_len bits
                record_int = record_int >> (bit_len)

            records.append(record)
            # Print a "." every 10000 records to indicate script is running.
            if i % 10000 == 0:
                print(".", end="", flush=True)
        print()

        records = np.array(records)

        # Save raw records to file as integers without tick rollover removed.
        if trbl_sht["raw_rlovr"]:
            np.savetxt(
                "raw_records_rollover.txt", records, fmt="%d", header="", comments=""
            )

        # Shift the tick count if necessary, so that it relates to the first
        # sample in each record (instead of the last).
        if logger.timing == "first":
            first_tic = 0
        elif logger.timing == "last":
            first_tic = int((logger.smpls_per_rec - 1) * logger.sample_epoch)
        last_field = len(logger.rec_fmt) - 1
        ticks_ms = records[:, last_field - logger.fmt_field["tic"]] - first_tic

        nominal_begin_tick = rec_begin * logger.record_epoch
        # print(f'Tick field length (in bits): {tic_field_len}')

        # If time tick values are not pressent then populate tick values with
        # assumed nominal tick count.
        if ticks_ms[-1] <= 0 and ticks_ms[1] <= 0:
            print(
                "ATTENTION!!! It appears that time-tick values were not recorded "
                "in the raw data file. All time values in the output are only "
                "as accurate as the PCB oscillator. Values from the  precision "
                "clock are not available!"
            )
            record_epoch = logger.sample_epoch * logger.smpls_per_rec
            stop = nominal_begin_tick + (ticks_ms.size) * record_epoch
            ticks_ms = np.arange(nominal_begin_tick, stop, record_epoch)
            records[:, last_field - logger.fmt_field["tic"]] = ticks_ms
            return records

        # Remove tick count rollovers and make actual ticks continuously
        # increasing.
        rollover_period = 2**logger.tic_bit_len  # in millisec
        # print(f'Rollover length (in millisec/ticks): {rollover_period}')
        # The number of rollovers prior to the beginning of the specified data
        # window.
        nom_rollovers_begin = int(nominal_begin_tick / rollover_period)
        nom_rollover_balance = nominal_begin_tick % rollover_period
        # Does the nominal rollover count of the first record align with the
        # actual count. (Within 1e7 millisec or approx 166 minutes.)
        if abs(ticks_ms[0] - nom_rollover_balance) < 1e7:
            actl_rollovers_begin = nom_rollovers_begin
        elif ticks_ms[0] - nom_rollover_balance > 1e7:
            # This indicates that a nominal rollover should have occurred but
            # the actual rollover hasn't yet.
            actl_rollovers_begin = nom_rollovers_begin - 1
        else:
            # This indicates that a nominal rollover should not have occurred
            # yet but the actual rollover has already.
            actl_rollovers_begin = nom_rollovers_begin + 1

        # {rollovers} contains index of the first record after each rollover.
        rollovers = np.where(ticks_ms[:-1] > ticks_ms[1:])[0] + 1
        cumtv_rollovers = actl_rollovers_begin
        if cumtv_rollovers != 0:
            ticks_ms = ticks_ms + rollover_period * cumtv_rollovers

        for idx, rollover in np.ndenumerate(rollovers):
            if rollover == rollovers[-1]:
                nxt_rollover = nrecs_want
            else:
                nxt_rollover = rollovers[idx[0] + 1]
            # print(rollover, nxt_rollover)
            if (ticks_ms[rollover + 1] - ticks_ms[rollover]) != 0:
                # Two consecutive identical tick counts indicates recording
                # has stopped.
                if nxt_rollover - rollover == 1:
                    # If the tick count does not rollover cleanly two
                    # consecutive tick records indicate a rollover (ie it takes
                    # more than onerecord to reset to zero), then for first
                    # record after the rollover calc as the previous cumulative
                    # tick count plus a std record period.
                    ticks_ms[rollover] = ticks_ms[rollover - 1] + logger.record_epoch
                elif (
                    abs(
                        (ticks_ms[rollover + 1] - ticks_ms[rollover - 1])
                        - 2 * logger.record_epoch
                    )
                    < 2
                ):
                    # If the record immediately before and after the currently
                    # indicated rollover are within 2ms of the expected time
                    # diff of 2 epochs, then the current single time tick is
                    # corrupt and not an actual rollover.
                    ticks_ms[rollover] = ticks_ms[rollover - 1] + logger.record_epoch
                elif (
                    abs(
                        (ticks_ms[rollover] - ticks_ms[rollover - 2])
                        - 2 * logger.record_epoch
                    )
                    < 2
                ):
                    # If the currently indicated rollover record and two records
                    # previous are within 2ms of the expected time diff of
                    # 2 epochs, then the previous single time tick is
                    # corrupt and not an actual rollover.
                    ticks_ms[rollover - 1] = (
                        ticks_ms[rollover - 2] + logger.record_epoch
                    )
                else:
                    cumtv_rollovers = cumtv_rollovers + 1
                    ticks_ms[rollover:nrecs_want] = (
                        ticks_ms[rollover:nrecs_want] + rollover_period
                    )
                    rollover_dt = clk_start_dt + np.timedelta64(
                        ticks_ms[rollover], "ms"
                    )
                    print(f"A time tick rollover occurred at " f"{rollover_dt}.")
    records[:, last_field - logger.fmt_field["tic"]] = ticks_ms

    return records


###############################################################################
def clockdrift(
    apg_filename, logger, clk_start_dt, gpssync_dt, sync_tick_count, trbl_sht
):
    """Calculates the clock drift using one of two methods.

    The expected number of tick counts between clk_start_dt and gpssync_dt
    will be calculated. The size of the tick record in logger['tic_bit_len'] is
    used to determine when the tick count 'rolls over' to zero.
    If sync_tick_count is provided, the difference between this and the
    expected value gives the drift in ticks (milliseconds).
    If sync_tick_count is not present, it is assumed that a fixed frequency has
    been injected into the raw data precisely at gpssync_dt. The precise tick
    count when this frequency starts is detected in the data and this value
    is used in place of sync_tick_count.
    """
    millisecs_logged = dt64_utils.delta64_to_ms(gpssync_dt - clk_start_dt)

    if sync_tick_count is None:
        # Assign names to column numbers of raw data array.
        # Note that raw data array columns are reverse order to raw binary.
        last_field = len(logger.rec_fmt) - 1
        tick_col = last_field - logger.fmt_field["tic"]
        pcore_col = last_field - logger.fmt_field["pcore"]
        pn_col = last_field - logger.fmt_field["pn"]

        # Window for identifying sync_tick_count is +/-5 minutes long.
        sync_wndw_ms = 5 * 60000
        # GPS sync time (gpssync_dt) is mid point of  window for sync.
        wndw_begin_ms = dt64_utils.delta64_to_ms(gpssync_dt - clk_start_dt)
        wndw_begin_ms = wndw_begin_ms - sync_wndw_ms
        rec_begin = int(wndw_begin_ms / logger.record_epoch)
        nrecs_want = int(sync_wndw_ms * 2 / logger.record_epoch) + 1

        sync_records = extractrecords(
            apg_filename, logger, nrecs_want, rec_begin, trbl_sht, clk_start_dt
        )

        # Save timesync records to file as integers with tick rollover removed.
        if trbl_sht["raw_sync"]:
            np.savetxt(
                "raw_sync_records.txt", sync_records, fmt="%d", header="", comments=""
            )

        # Identify the start of the record block where the pressure values
        # start changing again (ie This is where frequency injection for time
        # sync occurs).
        # Identify all consecutive row pairs where p_core changes.
        pcore_diff = np.diff(sync_records[:, pcore_col])
        pcore_diff_row = (np.where(pcore_diff != 0)[0]) + 1
        # Select the final instance where p_core starts changing
        # This exculdes any single noise values occuring before actual
        # frequency injection.
        diff_row_increments = np.diff(pcore_diff_row)
        x = np.where(diff_row_increments > 1)[0] + 1
        x = np.insert(x, 0, 0)
        poss_sync_row = pcore_diff_row[x] - 1

        # For each poss_sync_row check:
        #   - immed prev row has all Pn values as zero,
        #      (Indicates two consecutive p_core values to be identical
        #       althoughsurrounding values continue changing)
        #   - immed next row does not have all Pn values as zero.
        #       (Indicates noise value occuring before actual frequency
        #        injection.)
        # It is possible for two consecutive p_core values
        # to be identical although surrounding values continue changing.
        for n in range(np.size(poss_sync_row) - 1, -1, -1):
            sync_row = poss_sync_row[n]
            prev_sync_row = sync_row - 1
            prev_block = sync_records[prev_sync_row, :]
            next_sync_row = sync_row + 1
            next_block = sync_records[next_sync_row, :]

            if pn_col > pcore_col:
                prev_pn = prev_block[pcore_col + 1 : pn_col + 1]
                next_pn = next_block[pcore_col + 1 : pn_col + 1]
            else:
                prev_pn = prev_block[pcore_col - 1 : pn_col - 1 : -1]
                next_pn = next_block[pcore_col - 1 : pn_col - 1 : -1]

            prev_nonzero = np.where(prev_pn != 0)[0]
            next_nonzero = np.where(next_pn != 0)[0]
            if not prev_nonzero.any() and next_nonzero.any():
                sync_row = poss_sync_row[n]
                break

        try:
            sync_block = sync_records[sync_row, :]
        except UnboundLocalError:
            sys.exit(
                f"Unable to determine clock drift.\n"
                f"The raw data in during the period {gpssync_dt} "
                f"+/-{sync_wndw_ms/1000} seconds does not contain "
                f"a frequency injection for syncing to."
            )

        if pn_col > pcore_col:
            pn = sync_block[pcore_col + 1 : pn_col + 1]
        else:
            pn = sync_block[pcore_col - 1 : pn_col - 1 : -1]

        nonzero = np.where(pn != 0)[0]
        if nonzero.any():
            i = logger.smpls_per_rec - nonzero.size
        else:
            i = logger.smpls_per_rec

        sync_tick_count = sync_block[tick_col] + (i * logger.sample_epoch)
        sync_tick_count = int(sync_tick_count)

    else:
        # Number of ticks_ms until rollover and restart at zero.
        tick_rollover = 2**logger.tic_bit_len  # 1 tick = 1 millisecond
        millisecs_logged = millisecs_logged % tick_rollover

    # APG logger clock offset relative to GPS time (positive = APG > GPS)
    clk_drift_at_end = int(sync_tick_count - millisecs_logged)

    return clk_drift_at_end


def remove_noise_meddiff(
    raw_data,
    mask_data,
    millisecs,
    millisecs_all,
    bin_size,
    tolerance,
):
    """Generate refined spike removal.

    Bin data and take median of each bin (bin_size in milliseconds), then
    interpolate back to size of full dataset.
    Take difference of binned median and raw values. Where difference (tollerance)
    is greater than specified amount delete value and corresponding time stamp.
    Replace deleted data points with interpolated values.
    """
    bins = int((millisecs_all.max() - millisecs_all.min()) / bin_size)
    binned_tptr, bin_edges, _ = stat.binned_statistic(
        millisecs, mask_data, "median", bins
    )
    bin_centers = (bin_edges[1:] + bin_edges[:-1]) / 2
    tptr_medsmth = np.interp(millisecs, bin_centers, binned_tptr)
    difference = np.abs(mask_data - tptr_medsmth)
    mask = np.where(difference > tolerance)
    raw_data = np.delete(raw_data, mask)
    millisecs = np.delete(millisecs, mask)
    return raw_data, millisecs
    # return np.interp(millisecs_all, millisecs, raw_data)


###############################################################################
if __name__ == "__main__":
    main()
