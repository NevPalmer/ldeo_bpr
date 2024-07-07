"""Module for handling BPR/APG raw data."""

import sys
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from . import constants as const
from . import dt64_utils
from .logger import Logger


@dataclass(frozen=False)
class RawFile:
    """Class for managing BPR/APG raw data file."""

    filename: Path
    logger: Logger
    start_clk: np.datetime64
    gpssync_dt: np.datetime64 = None
    sync_tick_count: int = None
    end_clk: np.datetime64 = field(init=False)
    filesize_b: int = field(init=False)
    num_rcrds: int = field(init=False)
    file_duration_ms: int = field(init=False)
    clockdrift_ms: int = field(init=False)

    def __post_init__(self):
        """Generate calculated attributes."""
        # Calculate duration and end time of raw file.
        try:
            self.filesize_b = self.filename.stat().st_size
        except FileNotFoundError:
            sys.exit(f'Raw APG file "{self.filename}" does not exist.')

        self.num_rcrds = (
            int((self.filesize_b - self.logger.head_len) / self.logger.rec_len) - 1
        )
        self.file_duration_ms = self.num_rcrds * self.logger.record_epoch
        self.end_clk = self.start_clk + np.timedelta64(self.file_duration_ms, "ms")
        self._clockdrift()

    def _clockdrift(self):
        """Calculates the clock drift using one of two methods.

        The expected number of tick counts between self.start_clk and gpssync_dt
        will be calculated. The size of the tick record in logger.tic_bit_len is
        used to determine when the tick count 'rolls over' to zero.
        If sync_tick_count is provided, the difference between this and the
        expected value gives the drift in ticks (milliseconds).
        If sync_tick_count is not present, it is assumed that a fixed frequency has
        been injected into the raw data precisely at gpssync_dt. The precise tick
        count when this frequency starts is detected in the data and this value
        is used in place of sync_tick_count.
        """
        if self.gpssync_dt is None:
            # No clock drift calculated, as no sync parameters provided.
            self.clockdrift_ms = None
            self.gpssync_dt = self.end_clk
            return

        millisecs_logged = dt64_utils.delta64_to_ms(self.gpssync_dt - self.start_clk)

        if self.sync_tick_count is not None:
            # Calculate clock drift from sync_tick_count provided.
            # Number of ticks_ms until rollover and restart at zero.
            tick_rollover = 2**self.logger.tic_bit_len  # 1 tick = 1 millisecond
            millisecs_logged = millisecs_logged % tick_rollover
            # APG logger clock offset relative to GPS time (positive = APG > GPS)
            self.clockdrift_ms = int(self.sync_tick_count - millisecs_logged)
            return

        # Calculate clock drift by detecting frequency injection.
        # Assign names to column numbers of raw data array.
        # Note that raw data array columns are reverse order to raw binary.
        last_field = len(self.logger.rec_fmt) - 1
        tick_col = last_field - self.logger.fmt_field["tic"]
        pcore_col = last_field - self.logger.fmt_field["pcore"]
        pn_col = last_field - self.logger.fmt_field["pn"]

        # Window for identifying sync_tick_count is 30 seconds long.
        sync_wndw_ms = 30_000
        # GPS sync time (gpssync_dt) is mid point of  window for sync.
        wndw_begin_ms = dt64_utils.delta64_to_ms(self.gpssync_dt - self.start_clk)
        wndw_begin_ms = wndw_begin_ms - sync_wndw_ms
        start_rcrd = int(wndw_begin_ms / self.logger.record_epoch)
        num_rcrds_wanted = int(sync_wndw_ms * 2 / self.logger.record_epoch) + 1

        sync_records = extract_records(
            self.filename,
            self.logger,
            num_rcrds_wanted,
            start_rcrd,
            const.TROUBLE_SHOOT,
            self.start_clk,
        )

        # Save timesync records to file as integers with tick rollover removed.
        if const.TROUBLE_SHOOT["raw_sync"]:
            np.savetxt(
                "raw_sync_records.txt",
                sync_records,
                fmt="%d",
                header="",
                comments="",
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
                f"The raw data in during the period {self.gpssync_dt} "
                f"+/-{sync_wndw_ms/1000} seconds does not contain "
                f"a frequency injection for syncing to."
            )

        if pn_col > pcore_col:
            pn = sync_block[pcore_col + 1 : pn_col + 1]
        else:
            pn = sync_block[pcore_col - 1 : pn_col - 1 : -1]

        nonzero = np.where(pn != 0)[0]
        if nonzero.any():
            i = self.logger.smpls_per_rec - nonzero.size
        else:
            i = self.logger.smpls_per_rec

        self.sync_tick_count = int(
            sync_block[tick_col] + (i * self.logger.sample_epoch)
        )

        # APG logger clock offset relative to GPS time (positive = APG > GPS)
        self.clockdrift_ms = int(self.sync_tick_count - millisecs_logged)
        return


def extract_records(
    raw_filename: Path,
    file_start_clk: np.datetime64,
    logger: Logger,
    start_rcrd: int,
    num_rcrds_wanted: int,
):
    """Extracts binary records from a raw APG data logger file."""
    if const.TROUBLE_SHOOT["binary_out"]:
        binary_filename = "raw_binary.txt"
        # Create empty file, overwrite if exists.
        open(binary_filename, "w", encoding="utf8").close()
    if const.TROUBLE_SHOOT["hex_out"]:
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
            if const.TROUBLE_SHOOT["binary_out"]:
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
            if const.TROUBLE_SHOOT["hex_out"]:
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
        if const.TROUBLE_SHOOT["raw_rlovr"]:
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
