"""Module for handling BPR/APG raw data."""

import sys
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

import ldeo_bpr as bpr
from ldeo_bpr import dt64_utils


@dataclass(frozen=False)
class RawFile:
    """Class for managing BPR/APG raw data file."""

    filename: Path
    logger: bpr.Logger
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

        sync_records = bpr.extract_records(
            self.filename,
            self.logger,
            num_rcrds_wanted,
            start_rcrd,
            bpr.TROUBLE_SHOOT,
            self.start_clk,
        )

        # Save timesync records to file as integers with tick rollover removed.
        if bpr.TROUBLE_SHOOT["raw_sync"]:
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
