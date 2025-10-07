"""A script for extracting raw data from LDEO type APG data loggers."""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import obspy
import scipy.signal as sig
import scipy.stats as stat

from ldeo_bpr import arg_parser, dt64_utils, raw_data
from ldeo_bpr import constants as const
from ldeo_bpr.logger import Logger
from ldeo_bpr.paros import Paros


def main():
    """The first function run when this script is run directly."""
    # Retrieve CLI parameters.
    args = arg_parser.parse_args_apg_read()

    # Read Paros transducer coefficients from .ini file.
    paros = Paros.from_file(filename=args.apgini, paros_sn=args.snapg)

    # Read APG logger configuration parameters from .ini file.
    logger = Logger.from_file(
        filename=args.loggerini, logger_version=args.loggerversion
    )

    # Create a BPR raw data file object.
    raw_file = raw_data.RawFile(
        filename=args.infile,
        logger=logger,
        start_clk=args.clkstart,
        gpssync_dt=args.gpssynctime,
        sync_tick_count=args.synctickcount,
    )

    print("=" * 80)
    print(f"STATS OF RAW FILE: {raw_file.filename}")
    print(f"Time of first sample = {raw_file.start_clk}")
    print(f"Time of last sample = {raw_file.end_clk}")
    file_duration_days = raw_file.actl_file_tics_ms / (24 * 3600000)
    print(f"File duration = {file_duration_days} days")
    print(f"Filesize = {raw_file.filesize_b} bytes")
    print(f"Number of records = {raw_file.num_rcrds:d}")
    print(
        f"NOTE: Times and dates given above are based on "
        f"the Precise tick count values adjusted for the "
        f"drift of the precision time base module.\n"
        f"Precise tick counts are generated from the precision time base "
        f"module/clock (eg CSAC, Seascan, etc).\n"
        f"'Nominal' tick counts are calculated by counting the number "
        f"of sample epochs multiplied by the sample period "
        f"({logger.sample_epoch} millisecs).\n"
        f"For some APG loggers (eg CSAC) the Nominal and Precise tick counts are "
        f"identical because they are both derived from the precision time base "
        f"module.\n"
        f"For other loggers (eg Seascan) the Nominal timing is derived from an "
        f"imprecise PCB oscillator which drifts relative to the Precise "
        f"tick count. This 'Nominal' tick difference is corrected contuously in the "
        f"results by interpolating every record.\n"
        f"Nominal tick count difference at end of recording (millisecs):         "
        f"{raw_file.nom_tick_diff_ms}\n"
        f"   (Nominal - Precise tick counts)"
    )
    if not raw_file.clockdrift_ms:
        print(
            "No clock drift for precision time base has been calculated or applied.\n"
            "   Insufficient parameters supplied."
        )
    else:
        print(
            f"Clock drift for precision time base at end of recording (millisecs):   "
            f"{raw_file.clockdrift_ms}\n"
            f"   (Precise tick counts - actual GPS time)"
        )

    print("=" * 80)
    print("STATS OF DATA WINDOW TO BE EXTRACTED:")
    if args.beginwndw is None:
        args.beginwndw = raw_file.start_clk
    if args.endwndw is None:
        args.endwndw = raw_file.end_clk
    print(f"Window beginning: {args.beginwndw}")
    wndw_len = args.endwndw - args.beginwndw
    wndw_len_ms = dt64_utils.delta64_to_ms(wndw_len)
    wndw_len_days = wndw_len_ms / (24 * 3600000)
    print(f"Window length (days): {wndw_len_days}")

    print("=" * 80)
    print("STATS FOR EACH TIME BIN:")

    # Loop to extract data in time bins
    if args.outfile:
        # Create empty file, overwrite if exists.
        open(args.outfile, "w", encoding="utf8").close()
    wndw_beg_unix_ms = dt64_utils.dt64_to_ms(args.beginwndw)
    wndw_end_unix_ms = dt64_utils.dt64_to_ms(args.endwndw)
    bin_beg_ms = wndw_beg_unix_ms
    bininterval_ms = dt64_utils.delta64_to_ms(args.bininterval)
    if bininterval_ms == 0:
        bin_end_ms = wndw_end_unix_ms
    else:
        bin_end_ms = bin_beg_ms - (bin_beg_ms % bininterval_ms) + bininterval_ms
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
            raw_file,
            args.network,
            args.station,
            bin_begin_dt,
            bin_end_dt,
            args.decimate,
            args.tempsmth,
            args.noisefilt,
            args.bitshift,
            args.fmttime,
            args.plot,
            args.outfile,
            args.mseedpath,
        )

        bin_beg_ms = bin_end_ms
        bin_end_ms = bin_beg_ms + bininterval_ms


###############################################################################
def generate_results(
    logger: Logger,
    paros: Paros,
    raw_file: raw_data.RawFile,
    nwk_name,
    stn_name,
    bin_begin_dt,
    bin_end_dt,
    decmt_intvl,
    tmptr_smth_fctr,
    noisefilt,
    bitshift,
    time_format,
    plot_flags,
    out_filename: Path,
    mseed_path: Path,
):
    """This is the primary function used to extract and output results."""
    # Calculate window for extracting data.
    if noisefilt:
        bin_padding = 2_500_000  # milliseconds
    else:
        bin_padding = (tmptr_smth_fctr - 1) * logger.smpls_per_rec * 10

    bin_begin_ms = dt64_utils.delta64_to_ms(bin_begin_dt - raw_file.start_clk)
    bin_len_ms = dt64_utils.delta64_to_ms(bin_end_dt - bin_begin_dt)
    bin_end_ms = bin_begin_ms + bin_len_ms
    if raw_file.clockdrift_ms:
        clockdrift_ms = raw_file.clockdrift_ms
    else:
        clockdrift_ms = 0
    nom_tick_correction = np.floor(
        raw_file.nom_tick_diff_ms * (bin_begin_ms / raw_file.actl_file_tics_ms)
    )
    drift_correction = np.floor(
        clockdrift_ms * (bin_begin_ms / raw_file.actl_file_tics_ms)
    )
    bin_begin_ms = bin_begin_ms + nom_tick_correction + drift_correction
    nom_tick_correction = np.ceil(
        raw_file.nom_tick_diff_ms * (bin_end_ms / raw_file.actl_file_tics_ms)
    )
    drift_correction = np.ceil(
        clockdrift_ms * (bin_begin_ms / raw_file.actl_file_tics_ms)
    )
    bin_end_ms = bin_end_ms + nom_tick_correction + drift_correction

    # Make sure the requested times don't fall outside available records.
    if (bin_begin_ms - bin_padding) < 0:
        padded_bin_begin_ms = 0
    else:
        padded_bin_begin_ms = bin_begin_ms - bin_padding

    padded_bin_len_ms = bin_end_ms - padded_bin_begin_ms + bin_padding
    avail_bin_len_ms = raw_file.nom_file_duration_ms - bin_begin_ms
    if (avail_bin_len_ms - padded_bin_len_ms) <= 0:
        padded_bin_len_ms = padded_bin_len_ms - bin_padding

    start_rcrd = int(padded_bin_begin_ms / logger.record_epoch)
    num_rcrds_wanted = int(padded_bin_len_ms / logger.record_epoch) + 1

    # Extract records from APG file for the time window specified.
    print("Extracting raw records:\n", end="")

    records = raw_data.extract_records(
        raw_file.filename,
        raw_file.start_clk,
        logger,
        start_rcrd,
        num_rcrds_wanted,
        bitshift,
    )

    # Save raw records to file as integers with tick rollover removed.
    if const.TROUBLE_SHOOT["raw_no_rlovr"]:
        np.savetxt(
            "raw_records_no-rollover.txt", records, fmt="%d", header="", comments=""
        )

    # Assign names to column numbers of raw data array.
    # Note that raw data array columns are reverse order to raw binary.
    last_col = len(logger.rec_fmt) - 1
    pcore_col = last_col - logger.fmt_field["pcore"]
    pn_col = last_col - logger.fmt_field["pn"]
    tptr_col = last_col - logger.fmt_field["tptr"]
    tics_col = last_col - logger.fmt_field["tic"]

    # Create an array for each raw observable (pressure, temperature, ticks_ms)
    if pn_col > pcore_col:
        press_raw = records[:, pcore_col : pn_col + 1]
    else:
        press_raw = records[:, pcore_col : pn_col - 1 : -1]
    press_raw = np.cumsum(press_raw, axis=1)
    press_raw = press_raw.reshape(num_rcrds_wanted * logger.smpls_per_rec)
    temp_raw = records[:, tptr_col]
    ticks_ms = records[:, tics_col]

    actual_end_tick = ticks_ms[-1]
    actual_begin_tick = ticks_ms[0]
    nominal_begin_tick = start_rcrd * logger.record_epoch
    nominal_end_tick = (start_rcrd + num_rcrds_wanted - 1) * logger.record_epoch
    # print(f"\nActual beginning tick:  {actual_begin_tick}")
    # print(f"Nominal beginning tick: {nominal_begin_tick}")
    # print(f"Actual end tick:        {actual_end_tick}")
    # print(f"Nominal end tick:       {nominal_end_tick}\n")

    nom_ticks_t = np.linspace(
        nominal_begin_tick, nominal_end_tick, num=num_rcrds_wanted, endpoint=True
    ).astype("int64")
    nom_ticks_p = np.linspace(
        nominal_begin_tick,
        nominal_end_tick + logger.record_epoch,
        num=num_rcrds_wanted * logger.smpls_per_rec,
        endpoint=False,
    ).astype("int64")

    # Write a summary file showing syncronised tic counts and exit.
    if const.TROUBLE_SHOOT["tic_sync"]:
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
        end_diff = nominal_end_tick - actual_end_tick
        print(
            f"Nominal tick count difference at start of time bin (millisecs):         "
            f"{beg_diff}\n"
            f"Nominal tick count difference at end of time bin (millisecs):           "
            f"{end_diff}\n"
            f"   (Nominal - Precise tick counts)\n"
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
    millisecs_logged = dt64_utils.delta64_to_ms(
        raw_file.gpssync_dt - raw_file.start_clk
    )
    if not raw_file.clockdrift_ms:
        raw_file.clockdrift_ms = 0
    drift_beg = raw_file.clockdrift_ms * (millisecs_t[0] / millisecs_logged)
    drift_end = raw_file.clockdrift_ms * (millisecs_t[-1] / millisecs_logged)
    drift_applied = (drift_beg + drift_end) / 2
    # Round the drift to be applied to the  nearest whole sample epoch.
    drift_applied = int(
        logger.sample_epoch * round(drift_applied / logger.sample_epoch, 0)
    )
    print(
        f"Clock drift of precision time base applied to the time bin (millisecs): "
        f"{drift_applied}\n"
        f"   (Precise tick count - drift = final time).\n",
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
    # tmptr_period_usec.fill(5.78714)  # Fixed temp of +5.71°C for Paros 150181
    # tmptr_period_usec.fill(5.822)  # Fixed temp of +29.2°C for Paros 150337

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
            bin_size=2_500_000,
            tolerance=0.125,
        )
        coef_uv_filt = tmptr_period_usec_filt - paros.u[0]
        temperature_filt = np.polyval(paros.y, coef_uv_filt)
        tmptr_period_usec_filt, millisecs_filt = remove_noise_meddiff(
            raw_data=tmptr_period_usec_filt,
            mask_data=temperature_filt,
            millisecs=millisecs_filt,
            millisecs_all=millisecs_t,
            bin_size=500_000,
            tolerance=0.025,
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
    pressure = pressure * const.PSIA_TO_PASCAL  # Convert pressure units
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
            bin_size=2_500_000,
            tolerance=200,
        )

        # Uncomment line below to plot first iteration pressure instead of raw pressure.
        # pressure_raw = pressure

        # Second iteration noise removal with tighter tollerances
        pressure_filt, millisecs_filt = remove_noise_meddiff(
            raw_data=pressure_filt,
            mask_data=pressure_filt,
            millisecs=millisecs_filt,
            millisecs_all=millisecs_p,
            bin_size=10_000,
            tolerance=50,
        )

        # Third iteration noise removal with tighter tollerances
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
    mask = np.logical_and(time >= bin_begin_ms, time < bin_end_ms)
    time = time[mask]
    pressure_out = pressure_out[mask]
    temperature_out = temperature_out[mask]

    # Convert timestamp from seconds to datetime if specified.
    if time_format == "d":
        print("Calculating a timestamp array from the milliseconds array.", flush=True)
        time = (raw_file.start_clk + dt64_utils.ms_to_delta64(time)).astype(
            "datetime64[ms]"
        )
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

        stats["channel"] = const.P_CHNL_CODE
        stats["sampling_rate"] = 1000 / logger.sample_epoch
        trace_p = obspy.Trace(data=pressure_mseed, header=stats)
        stream_p = obspy.Stream(traces=[trace_p])
        mseed_filename = f"{nwk_name}_{stn_name}_{stats['channel']}_{dt_text}.mseed"
        mseed_filename = mseed_path / mseed_filename
        print(f'Writing pressure data to MiniSEED file "{mseed_filename}".')
        stream_p.write(mseed_filename, format="MSEED")

        stats["channel"] = const.T_CHNL_CODE
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
        print("Generating plot.", flush=True)
        p_kpa = pressure_out / 1000
        ## Set min-max values for plot Y-axes
        ## Uncomment below to plot full range of data.
        # p_min = np.min(p_kpa)
        # p_max = np.max(p_kpa)
        # p_range = p_max - p_min

        # if p_range == 0:
        #     intvl = 10
        # else:
        #     int_log = int(np.log10(p_range))
        #     if 10**int_log / p_range >= 0.5:
        #         intvl = 10**int_log / 10
        #     elif 10**int_log / p_range >= 0.2:
        #         intvl = 10**int_log / 5
        #     elif 10**int_log / p_range >= 0.1:
        #         intvl = 10**int_log / 2
        #     else:
        #         intvl = 10**int_log
        # p_min = p_min - p_min % intvl - intvl
        # p_max = p_max - p_max % intvl + 2 * intvl

        # t_min = np.min(temperature_out)
        # t_max = np.max(temperature_out)
        # t_range = t_max - t_min
        # if t_range == 0:
        #     intvl = 0.1
        # else:
        #     int_log = int(np.log10(t_range))
        #     if 10**int_log / t_range >= 0.5:
        #         intvl = 10**int_log / 10
        #     elif 10**int_log / t_range >= 0.2:
        #         intvl = 10**int_log / 5
        #     elif 10**int_log / t_range >= 0.1:
        #         intvl = 10**int_log / 2
        #     else:
        #         intvl = 10**int_log
        # t_min = t_min - t_min % intvl - intvl
        # t_max = t_max - t_max % intvl + 2 * intvl

        ## Uncomment below to plot specified range of data about the mean.
        p_median = np.median(p_kpa)
        # pres_median = 22270
        p_spread = 20
        p_min = p_median - (p_spread / 2)
        p_max = p_median + (p_spread / 2)

        t_median = np.median(temperature_out)
        # tptr_median = 5.7
        t_spread = 1
        t_min = t_median - (t_spread / 2)
        t_max = t_median + (t_spread / 2)

        # Plot Results
        fig, ax1 = plt.subplots(figsize=(15, 9))
        plt.ylim(p_max, p_min)
        ax2 = ax1.twinx()
        plt.ylim(t_min, t_max)

        # Plot raw pressure values if requested
        if plot_flags["format"] == "r":
            if time_format == "d":
                time_p = raw_file.start_clk + dt64_utils.ms_to_delta64(millisecs_p)
            else:
                time_p = millisecs_p / 1000
            ax1.plot(
                time_p,
                pressure_raw / 1000,
                color="pink",
                marker=".",
                markersize=1.0,
                linestyle="",
            )

        # Plot raw temperature values if requested
        if plot_flags["format"] == "r":
            if time_format == "d":
                time_t = raw_file.start_clk + dt64_utils.ms_to_delta64(millisecs_t)
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
        ax1.set_ylabel("Pressure (kPa)", color=color)
        ax1.tick_params(axis="y", labelcolor=color)
        ax1.grid(axis="x")
        ax1.plot(
            time,
            p_kpa,
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

        fig.suptitle(
            f"{raw_file.filename}\nBegin: {bin_begin_dt}  -  End: {bin_end_dt}"
        )
        plt.tight_layout()
        if time_format == "d":
            # Rotates and aligns the X-axis labels.
            plt.gcf().autofmt_xdate(bottom=0.2, rotation=30, ha="right")
        if plot_flags["output"] in ["s", "b"]:
            basename = raw_file.filename.stem
            fig.savefig(f"./{basename}.png", dpi=200)
        if plot_flags["output"] in ["d", "b"]:
            plt.show()
        plt.close(fig)

    return


###############################################################################
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
    try:
        binned_data, bin_edges, _ = stat.binned_statistic(
            millisecs, mask_data, "median", bins
        )
    except ValueError as err:
        print(
            "The data window or bin selected is too short for the filter "
            "parameters currently specified in the script, or there is "
            "insufficient usable remianing after filtering. Try selecting "
            "a longer data sample or edit the bin size of 'remove_noise_meddiff' "
            "filter in module 'apg_read.py'. Raw filtered data has been "
            "processed as zeros."
        )
        raw_data = millisecs_all.copy()
        raw_data.fill(0.0)
        return raw_data, millisecs_all
    bin_centers = (bin_edges[1:] + bin_edges[:-1]) / 2
    tptr_medsmth = np.interp(millisecs, bin_centers, binned_data)
    difference = np.abs(mask_data - tptr_medsmth)
    mask = np.where(difference > tolerance)
    raw_data = np.delete(raw_data, mask)
    millisecs = np.delete(millisecs, mask)
    return raw_data, millisecs


###############################################################################
if __name__ == "__main__":
    main()
