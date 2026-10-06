"""Osmosis plots and readouts (counts per side, chemical potentials, wells)."""

from __future__ import annotations

import math

import numpy as np
from PySide6 import QtCore, QtWidgets
import pyqtgraph as pg

from ..measurements.osmosis import OsmosisFrame
from .instruments import BACKGROUND, MUTED, TEXT, _plot, fmt

LEFT_COLOR = "#76d6c4"       # host side; bound discs use the same colour
RIGHT_COLOR = "#e8b44d"      # disc-only side
BOUND_COLOR = "#ab9edb"
REFERENCE_STYLE = QtCore.Qt.PenStyle.DashLine


def _binomial(n: int, p: float, k: np.ndarray) -> np.ndarray:
    logs = (math.lgamma(n + 1) - np.array([math.lgamma(x + 1) + math.lgamma(n - x + 1)
                                            for x in k])
            + k*math.log(max(p, 1e-300)) + (n - k)*math.log(max(1 - p, 1e-300)))
    return np.exp(logs)


class OsmosisPanel(QtWidgets.QWidget):
    """Four plots and a summary line; the ideal-point reference is dashed."""

    def __init__(self):
        super().__init__()
        self.setStyleSheet(f"background:{BACKGROUND};color:{TEXT}")
        layout = QtWidgets.QVBoxLayout(self)
        layout.setContentsMargins(4, 4, 4, 4)
        self.summary = QtWidgets.QLabel("—")
        self.summary.setWordWrap(True)
        layout.addWidget(self.summary)
        grid = QtWidgets.QGridLayout()
        layout.addLayout(grid, stretch=1)
        dash = dict(width=1, style=REFERENCE_STYLE)

        self.counts = _plot("Discs on each side", "Physical time", "count")
        self.left_line = self.counts.plot(pen=pg.mkPen(LEFT_COLOR, width=2))
        self.right_line = self.counts.plot(pen=pg.mkPen(RIGHT_COLOR, width=2))
        self.left_ref = self.counts.plot(pen=pg.mkPen(LEFT_COLOR, **dash))
        self.right_ref = self.counts.plot(pen=pg.mkPen(RIGHT_COLOR, **dash))

        self.potentials = _plot("Chemical potential μ/T", "Physical time", "ln c")
        self.mu_right = self.potentials.plot(pen=pg.mkPen(RIGHT_COLOR, width=2))
        self.mu_free = self.potentials.plot(pen=pg.mkPen(LEFT_COLOR, width=2))
        self.mu_bound = self.potentials.plot(pen=pg.mkPen(BOUND_COLOR, width=1.5))

        self.bound = _plot("Discs bound in wells", "Physical time", "count")
        self.bound_line = self.bound.plot(pen=pg.mkPen(BOUND_COLOR, width=2))
        self.bound_ref = self.bound.plot(pen=pg.mkPen(BOUND_COLOR, **dash))

        self.histogram = _plot("Discs per well", "n", "fraction of samples")
        self.bars = pg.BarGraphItem(x=[], height=[], width=.7, brush=BOUND_COLOR)
        self.histogram.addItem(self.bars)
        self.histogram_ref = self.histogram.plot(
            pen=pg.mkPen(TEXT, **dash), symbol="o", symbolSize=5, symbolBrush=TEXT)

        # Wall-temperature changes, marked on the time plots.
        self._steps: list[tuple[float, list]] = []

        for i, plot in enumerate((self.counts, self.potentials, self.bound, self.histogram)):
            plot.setMinimumHeight(150)
            grid.addWidget(plot, i // 2, i % 2)
        legend = QtWidgets.QLabel(
            f"<span style='color:{LEFT_COLOR}'>━</span> host side (left)  "
            f"<span style='color:{RIGHT_COLOR}'>━</span> discs only (right)  "
            f"<span style='color:{BOUND_COLOR}'>━</span> bound (μ: ideal-well estimate)  "
            "┄ ideal-point reference  ┆ wall temperature change")
        legend.setStyleSheet(f"color:{MUTED}")
        layout.addWidget(legend)

    def clear(self):
        for line in (self.left_line, self.right_line, self.left_ref, self.right_ref,
                     self.mu_right, self.mu_free, self.mu_bound, self.bound_line,
                     self.bound_ref, self.histogram_ref):
            line.setData([], [])
        self.bars.setOpts(x=[], height=[])
        self._set_steps([])
        self.summary.setText("—")

    def _set_steps(self, steps: list[tuple[float, float]]):
        """Show a labelled vertical line at each (time, new temperature)."""
        if [t for t, _ in self._steps] == [t for t, _ in steps]:
            return
        for _, lines in self._steps:
            for plot, line in lines:
                plot.removeItem(line)
        self._steps = []
        for t, temperature in steps:
            lines = []
            for plot in (self.counts, self.potentials, self.bound):
                line = pg.InfiniteLine(t, angle=90, pen=pg.mkPen(MUTED, width=1,
                                                                 style=QtCore.Qt.PenStyle.DotLine),
                                       label=f"T={temperature:g}" if plot is self.counts else None,
                                       labelOpts={"position": .92, "color": MUTED})
                plot.addItem(line)
                lines.append((plot, line))
            self._steps.append((t, lines))

    def update_frame(self, frame: OsmosisFrame | None, discs: int, plots: bool = True):
        if frame is None:
            self.clear()
            return
        current, reference, averages = frame.current, frame.reference, frame.averages
        pressure = averages.get("osmotic_pressure")
        mean = (f"mean since t = {averages['since']:.0f}: "
                f"L:R {averages['left']:.1f} : {averages['right']:.1f}, "
                f"bound {averages['bound']:.1f}" if averages else "averaging…")
        self.summary.setText(
            f"t = {current.time:.1f}   left : right = {current.left} : {current.right}   "
            f"bound {current.bound}   μ/T right {fmt(current.mu_right, 3)}, "
            f"free left {fmt(current.mu_free_left, 3)}<br>"
            f"<span style='color:{MUTED}'>{mean}; osmotic pressure "
            f"{fmt(pressure, 3)} (host contacts {fmt(averages.get('host_contact_part'), 3)}, "
            f"disc contacts {fmt(averages.get('disc_contact_part'), 2)}; ideal hosts "
            f"{fmt(reference['osmotic_pressure_ideal'], 3)}); "
            f"reference L:R {reference['left']:.1f} : {reference['right']:.1f}, "
            f"bound {reference['bound']:.1f}</span>")
        if not plots:
            return
        history = frame.history
        t = np.array([s.time for s in history])
        self.left_line.setData(t, [s.left for s in history])
        self.right_line.setData(t, [s.right for s in history])
        span = np.array([t[0], t[-1]]) if len(t) else np.array([])
        self.left_ref.setData(span, np.full(len(span), reference["left"]))
        self.right_ref.setData(span, np.full(len(span), reference["right"]))
        self.mu_right.setData(t, [s.mu_right for s in history], connect="finite")
        self.mu_free.setData(t, [s.mu_free_left for s in history], connect="finite")
        self.mu_bound.setData(t, [s.mu_bound for s in history], connect="finite")
        self.bound_line.setData(t, [s.bound for s in history])
        self.bound_ref.setData(span, np.full(len(span), reference["bound"]))
        self._set_steps([(b.time, b.temperature) for a, b in zip(history[:-1], history[1:])
                         if b.temperature != a.temperature])
        counts = frame.occupancy_histogram
        if counts.sum():
            n = np.arange(len(counts))
            self.bars.setOpts(x=n, height=counts/counts.sum())
            k = np.arange(len(counts) + 3)
            self.histogram_ref.setData(k, _binomial(discs, reference["per_host"]/discs, k))
