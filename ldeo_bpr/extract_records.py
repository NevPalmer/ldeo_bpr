"""Module for extracting data records from raw BPR/APG file."""

from pathlib import Path

import numpy as np

import ldeo_bpr as bpr
from ldeo_bpr.logger import Logger


def extract_records(
    raw_filename: Path,
    file_start_clk: np.datetime64,
    logger: Logger,
    start_rcrd: int,
    num_rcrds_wanted: int,
):
    """Extracts binary records from a raw APG data logger file."""
    if bpr.TROUBLE_SHOOT["binary_out"]:
        binary_filename = "raw_binary.txt"
        # Create empty file, overwrite if exists.
        open(binary_filename, "w", encoding="utf8").close()
    if bpr.TROUBLE_SHOOT["hex_out"]:
        hex_filename = "raw_hex.txt"
        # Create empty file, overwrite if exists.
        open(hex_filename, "w", encoding="utf8").close()

    with open(raw_filename, "rb") as apgfile:
        begin_byte = logger.head_len + start_rcrd * logger.rec_len
        apgfile.seek(begin_byte, 0)
        records = []
        for i in range(0, num_rcrds_wanted):
            binary_record = apgfile.read(logger.rec_len)

            # Print record as a string of Binary values to file.
            if bpr.TROUBLE_SHOOT["binary_out"]:
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
            if bpr.TROUBLE_SHOOT["hex_out"]:
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
        if bpr.TROUBLE_SHOOT["raw_rlovr"]:
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

        nominal_begin_tick = start_rcrd * logger.record_epoch
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
                nxt_rollover = num_rcrds_wanted
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
                    ticks_ms[rollover:num_rcrds_wanted] = (
                        ticks_ms[rollover:num_rcrds_wanted] + rollover_period
                    )
                    rollover_dt = file_start_clk + np.timedelta64(
                        ticks_ms[rollover], "ms"
                    )
                    print(f"A time tick rollover occurred at " f"{rollover_dt}.")
    records[:, last_field - logger.fmt_field["tic"]] = ticks_ms

    return records
