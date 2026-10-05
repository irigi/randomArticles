"""Carnot cycle diagrams (p–V, T–S, scoreboard) and the readout column."""

from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np
from PySide6 import QtCore, QtGui, QtWidgets
import pyqtgraph as pg


BRANCHES = ("hot", "adiabatic_expansion", "cold", "adiabatic_compression")
BRANCH_COLORS = {
    "hot":"#f2a65a", "adiabatic_expansion":"#c9bddf",
    "cold":"#76bce3", "adiabatic_compression":"#9fd3a6",
}
BRANCH_LABELS = {
    "hot":"hot isotherm", "adiabatic_expansion":"adiabatic expansion",
    "cold":"cold isotherm", "adiabatic_compression":"adiabatic compression",
}
TEXT = "#d5dce4"
MUTED = "#8e99a8"
BACKGROUND = "#15191f"
GOOD = "#9fd3a6"
BAD = "#e88a8a"


def fmt(value, digits=4):
    if value is None or not math.isfinite(value):
        return "—"
    return f"{value:.{digits}g}"


@dataclass(frozen=True)
class InstrumentSeries:
    """Plot history as arrays; the live window and replay share this shape."""

    time: np.ndarray
    area: np.ndarray
    pressure: np.ndarray
    entropy: np.ndarray
    gas_energy: np.ndarray
    temperature_trans: np.ndarray
    temperature_rot: np.ndarray
    branch: np.ndarray

    @classmethod
    def from_history(cls, history) -> "InstrumentSeries":
        def column(name):
            return np.array([math.nan if getattr(s, name) is None else
                             getattr(s, name) for s in history], dtype=np.float64)
        return cls(column("time"), column("area"), column("pressure"),
                   column("entropy"), column("gas_energy"),
                   column("temperature_trans"), column("temperature_rot"),
                   np.array([s.branch or "" for s in history]))

    @classmethod
    def from_columns(cls, columns, stop: int, max_points: int = 6000
                     ) -> "InstrumentSeries":
        """Frames [0, stop), thinned to about ``max_points`` for drawing.

        The most recent ``max_points/2`` frames are always kept at full
        resolution; older frames are strided. Long replays stay responsive.
        """
        recent = max_points//2
        if stop <= max_points:
            index = slice(0, stop)
        else:
            older = stop-recent
            stride = -(-older//(max_points-recent))
            index = np.concatenate((np.arange(0, older, stride),
                                    np.arange(older, stop)))
        return cls(*(columns[name][index] for name in (
            "time", "area", "pressure", "entropy", "gas_energy",
            "temperature_trans", "temperature_rot", "branch")))


def _plot(title, bottom, left):
    plot = pg.PlotWidget(background=BACKGROUND)
    plot.setTitle(title, color=TEXT, size="10pt")
    plot.setLabel("bottom", bottom, color=MUTED)
    plot.setLabel("left", left, color=MUTED)
    plot.showGrid(x=True, y=True, alpha=.15)
    plot.setMenuEnabled(False)
    for axis in ("bottom", "left"):
        plot.getAxis(axis).enableAutoSIPrefix(False)
    return plot


def _branch_lines(plot, width=2):
    lines = {"older": plot.plot(pen=pg.mkPen(OLDER_CYCLES, width=1), connect="finite")}
    lines.update({name: plot.plot(pen=pg.mkPen(color, width=width), connect="finite")
                  for name, color in BRANCH_COLORS.items()})
    return lines


OLDER_CYCLES = "#3a4350"


def _set_branch_lines(lines, x, y, branch, recent):
    """Colour recent samples by branch; draw older cycles as one dim line."""
    def joined(mask):
        # Include the next sample so consecutive segments meet without gaps.
        return mask | np.concatenate(([False], mask[:-1]))

    for name, line in lines.items():
        mask = ~recent if name == "older" else (branch == name) & recent
        line.setData(x, np.where(joined(mask), y, np.nan), connect="finite")


class CycleDiagrams(QtWidgets.QWidget):
    """p–V and T–S diagrams over a per-cycle scoreboard."""

    def __init__(self):
        super().__init__()
        self.setStyleSheet(f"background:{BACKGROUND};color:{TEXT}")
        layout = QtWidgets.QVBoxLayout(self)
        layout.setContentsMargins(4, 0, 4, 4)
        layout.setSpacing(4)
        plots = QtWidgets.QHBoxLayout()
        layout.addLayout(plots, stretch=1)

        self.pv = _plot("p–V  ·  pressure on the piston", "Area A (2D volume)",
                        "Pressure P")
        self.pv_reference = {name: self.pv.plot(pen=pg.mkPen(
            color, width=1.5, style=QtCore.Qt.PenStyle.DashLine))
            for name, color in BRANCH_COLORS.items()}
        self.pv_lines = _branch_lines(self.pv)
        self.pv_current = pg.ScatterPlotItem(size=11, pen=pg.mkPen("#ffffff"),
                                             brush=pg.mkBrush("#f3e7c3"))
        self.pv.addItem(self.pv_current)
        plots.addWidget(self.pv)

        self.ts = _plot("T–S  ·  ideal-gas entropy estimate",
                        "Entropy S/N  (k_B = 1, from A and gas energy)",
                        "Gas temperature 2E/(fN)")
        self.ts_reference = {name: self.ts.plot(pen=pg.mkPen(
            color, width=1.5, style=QtCore.Qt.PenStyle.DashLine))
            for name, color in BRANCH_COLORS.items()}
        self.ts_lines = _branch_lines(self.ts)
        self.ts_current = pg.ScatterPlotItem(size=11, pen=pg.mkPen("#ffffff"),
                                             brush=pg.mkBrush("#f3e7c3"))
        self.ts.addItem(self.ts_current)
        plots.addWidget(self.ts)

        self.legend = QtWidgets.QLabel()
        self.legend.setStyleSheet(f"color:{MUTED}")
        self.legend.setText("   ".join(
            f"<span style='color:{BRANCH_COLORS[b]}'>━</span> {BRANCH_LABELS[b]}"
            for b in BRANCHES) + "   ┄ ideal Carnot cycle")
        layout.addWidget(self.legend)

        board = QtWidgets.QHBoxLayout()
        layout.addLayout(board)
        self.current_cycle = QtWidgets.QLabel()
        self.current_cycle.setTextFormat(QtCore.Qt.TextFormat.RichText)
        self.current_cycle.setMinimumWidth(260)
        self.current_cycle.setAlignment(QtCore.Qt.AlignmentFlag.AlignTop)
        board.addWidget(self.current_cycle)
        self.table = QtWidgets.QTableWidget(0, 6)
        self.table.setHorizontalHeaderLabels(
            ["Cycle", "Q_H in", "Q_C in", "W by gas", "η = W/Q_H", "Gas ΔE"])
        self.table.verticalHeader().setVisible(False)
        self.table.horizontalHeader().setSectionResizeMode(
            QtWidgets.QHeaderView.ResizeMode.Stretch)
        self.table.setEditTriggers(QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.setSelectionMode(QtWidgets.QAbstractItemView.SelectionMode.NoSelection)
        self.table.setFixedHeight(130)
        self.table.setStyleSheet(
            f"QTableWidget{{background:{BACKGROUND};color:{TEXT};gridline-color:#2a323c}}"
            f"QHeaderView::section{{background:#222a33;color:{MUTED};border:0;padding:2px}}")
        board.addWidget(self.table, stretch=1)

        self._reference = None
        self._entropy_origin = None
        self._shown_cycles = None

    def set_reference(self, reference):
        """Draw the analytical point-gas cycle; None clears it."""
        if reference is self._reference:
            return
        self._reference = reference
        if reference is None:
            for line in (*self.pv_reference.values(), *self.ts_reference.values()):
                line.setData([], [])
            self._entropy_origin = None
            return
        branches = {branch.name: branch for branch in reference.branches}
        f = reference.degrees_of_freedom
        hot, cold = reference.hot_temperature, reference.cold_temperature
        a1 = branches["hot"].area[0]
        a2 = branches["hot"].area[-1]
        self._entropy_origin = math.log(a1)+f/2*math.log(hot)
        width = math.log(a2/a1)
        rectangle = {"hot": ([0., width], [hot, hot]),
                     "adiabatic_expansion": ([width, width], [hot, cold]),
                     "cold": ([width, 0.], [cold, cold]),
                     "adiabatic_compression": ([0., 0.], [cold, hot])}
        for name, line in self.pv_reference.items():
            branch = branches[name]
            line.setData(branch.area, branch.pressure)
        for name, line in self.ts_reference.items():
            line.setData(*rectangle[name])
        self.ts.setTitle(f"T–S  ·  ideal-gas entropy estimate  ·  η_Carnot = "
                         f"{reference.efficiency:.3f}", color=TEXT, size="10pt")

    def update_series(self, series: InstrumentSeries, current, particles, dof,
                      recent_from=-math.inf):
        """Draw the history; samples before ``recent_from`` are dimmed."""
        recent = series.time >= recent_from
        # The trailing window is zero-padded until it first fills; those
        # partial pressures are not drawn.
        pressure = np.where(series.time >= current.pressure_window,
                            series.pressure, np.nan)
        _set_branch_lines(self.pv_lines, series.area, pressure,
                          series.branch, recent)
        temperature = 2*series.gas_energy/(dof*max(particles, 1))
        origin = self._entropy_origin or 0.
        _set_branch_lines(self.ts_lines, series.entropy-origin, temperature,
                          series.branch, recent)
        if current.area is not None and current.pressure is not None:
            self.pv_current.setData([current.area], [current.pressure])
        else:
            self.pv_current.setData([], [])
        if current.entropy is not None:
            self.ts_current.setData([current.entropy-origin],
                                    [2*current.gas_energy/(dof*max(particles, 1))])
        else:
            self.ts_current.setData([], [])

    def update_scoreboard(self, current, cycles, cycle_start, transient_cycles,
                          ideal_efficiency):
        if cycle_start is not None:
            heat_hot = current.heat_hot-cycle_start.heat_hot
            heat_cold = current.heat_cold-cycle_start.heat_cold
            # Gas first law: W by gas = Q - ΔE, exact up to the logged residual.
            work = heat_hot+heat_cold-(current.gas_energy-cycle_start.gas_energy)
            colour = GOOD if work >= 0 else BAD
            branch = BRANCH_LABELS.get(current.branch, current.branch or "—")
            self.current_cycle.setText(
                f"<div style='color:{MUTED}'>CYCLE {current.completed_cycles+1} · {branch}</div>"
                f"<div style='font-size:15pt'>"
                f"<span style='color:{BRANCH_COLORS['hot']}'>Q_H {heat_hot:+.4g}</span><br>"
                f"<span style='color:{BRANCH_COLORS['cold']}'>Q_C {heat_cold:+.4g}</span><br>"
                f"<span style='color:{colour}'>W {work:+.4g}</span></div>"
                f"<div style='color:{MUTED}'>η_Carnot {fmt(ideal_efficiency, 3)} · "
                "W &gt; 0: work delivered</div>")
        else:
            self.current_cycle.setText(f"<span style='color:{MUTED}'>No cycle data</span>")
        key = (len(cycles), transient_cycles)
        if key == self._shown_cycles:
            return
        self._shown_cycles = key
        self.table.setRowCount(len(cycles))
        for row, cycle in enumerate(cycles):
            work = -cycle.piston_work_on_gas
            efficiency = work/cycle.heat_hot if cycle.heat_hot > 0 else math.nan
            transient = row < transient_cycles
            values = (f"{row+1}{' (warm-up)' if transient else ''}",
                      fmt(cycle.heat_hot), fmt(cycle.heat_cold), fmt(work),
                      fmt(efficiency, 3), fmt(cycle.delta_gas_energy))
            for column, text in enumerate(values):
                item = QtWidgets.QTableWidgetItem(text)
                item.setTextAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
                colour = MUTED if transient else TEXT
                if column in (3, 4) and math.isfinite(work):
                    colour = GOOD if work >= 0 else BAD
                item.setForeground(QtGui.QColor(colour))
                self.table.setItem(row, column, item)
        if cycles:
            self.table.scrollToBottom()


class _Section(QtWidgets.QWidget):
    """A heading with label/value rows, optionally collapsible."""

    def __init__(self, title, rows, collapsible=False, tooltips=None):
        super().__init__()
        layout = QtWidgets.QVBoxLayout(self)
        layout.setContentsMargins(0, 4, 0, 0)
        layout.setSpacing(1)
        if collapsible:
            self.toggle = QtWidgets.QToolButton()
            self.toggle.setText(title)
            self.toggle.setCheckable(True)
            self.toggle.setArrowType(QtCore.Qt.ArrowType.RightArrow)
            self.toggle.setToolButtonStyle(
                QtCore.Qt.ToolButtonStyle.ToolButtonTextBesideIcon)
            self.toggle.setStyleSheet(f"color:{MUTED};border:0")
            layout.addWidget(self.toggle)
        else:
            self.toggle = None
            heading = QtWidgets.QLabel(title)
            heading.setStyleSheet(f"color:{MUTED};font-weight:bold;"
                                  "border-bottom:1px solid #2a323c")
            layout.addWidget(heading)
        self.body = QtWidgets.QWidget()
        grid = QtWidgets.QGridLayout(self.body)
        grid.setContentsMargins(6 if collapsible else 0, 0, 0, 0)
        grid.setVerticalSpacing(1)
        self.values = {}
        for row, name in enumerate(rows):
            label = QtWidgets.QLabel(name)
            value = QtWidgets.QLabel("—")
            value.setAlignment(QtCore.Qt.AlignmentFlag.AlignRight)
            value.setTextInteractionFlags(
                QtCore.Qt.TextInteractionFlag.TextSelectableByMouse)
            if tooltips and name in tooltips:
                label.setToolTip(tooltips[name])
                value.setToolTip(tooltips[name])
            grid.addWidget(label, row, 0)
            grid.addWidget(value, row, 1)
            self.values[name] = value
        layout.addWidget(self.body)
        if collapsible:
            self.body.setVisible(False)
            self.toggle.toggled.connect(self._expand)

    def _expand(self, expanded):
        self.body.setVisible(expanded)
        self.toggle.setArrowType(QtCore.Qt.ArrowType.DownArrow if expanded
                                 else QtCore.Qt.ArrowType.RightArrow)


class ReadoutPanel(QtWidgets.QWidget):
    """Grouped numeric readouts and the temperature strip, one screen tall."""

    def __init__(self):
        super().__init__()
        self.setMinimumWidth(320)
        self.setStyleSheet(f"background:{BACKGROUND};color:{TEXT}")
        outer = QtWidgets.QVBoxLayout(self)
        outer.setContentsMargins(4, 4, 4, 4)
        scroll = QtWidgets.QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QtWidgets.QFrame.Shape.NoFrame)
        outer.addWidget(scroll)
        content = QtWidgets.QWidget()
        layout = QtWidgets.QVBoxLayout(content)
        layout.setContentsMargins(4, 0, 4, 4)
        layout.setSpacing(2)
        scroll.setWidget(content)

        self.inspector = QtWidgets.QLabel()
        self.inspector.setWordWrap(True)
        self.inspector.setStyleSheet("background:#222a33;padding:6px")
        self.inspector.setVisible(False)
        layout.addWidget(self.inspector)

        tips = {
            "Work by gas: ledger / ∮P dA":
                "Event-ledger piston work, and the integral of the trailing "
                "pressure window over the cycle with its smoothing bound "
                "(not a statistical uncertainty).",
            "η_net ± SE / η_Carnot":
                "Net external output over hot heat for matched cycles after the "
                "warm-up cycles; SE is a four-block jackknife standard error.",
        }
        self.sections = [
            _Section("STATE", ("Branch / phase", "Time / cycles", "Area",
                               "Pressure (window)", "T translation / rotation",
                               "Reservoirs T_H / T_C", "Gas energy")),
            _Section("CUMULATIVE", ("Heat in from hot Q_H", "Heat in from cold Q_C",
                                    "Motor work / load output",
                                    "First-law residual")),
            _Section("LAST CYCLE", ("Q_H / Q_C", "Work by gas: ledger / ∮P dA",
                                    "External output", "Total residual"),
                     tooltips=tips),
            _Section("EFFICIENCY", ("η_net ± SE / η_Carnot", "Gate"),
                     tooltips=tips),
        ]
        for section in self.sections:
            layout.addWidget(section)
        self.values = {name: value for section in self.sections
                       for name, value in section.values.items()}

        self.temperatures = _plot("Temperature", "Physical time", "T")
        self.temperatures.setMinimumHeight(170)
        self.temp_trans = self.temperatures.plot(pen=pg.mkPen("#f0c46c", width=2))
        self.temp_rot = self.temperatures.plot(pen=pg.mkPen("#ab9edb", width=2))
        dash = QtCore.Qt.PenStyle.DashLine
        self.temp_hot = self.temperatures.plot(
            pen=pg.mkPen(BRANCH_COLORS["hot"], width=1, style=dash))
        self.temp_cold = self.temperatures.plot(
            pen=pg.mkPen(BRANCH_COLORS["cold"], width=1, style=dash))
        layout.addWidget(self.temperatures, stretch=1)
        legend = QtWidgets.QLabel(
            "<span style='color:#f0c46c'>━</span> translation  "
            "<span style='color:#ab9edb'>━</span> rotation  "
            f"<span style='color:{BRANCH_COLORS['hot']}'>┄</span> T_H  "
            f"<span style='color:{BRANCH_COLORS['cold']}'>┄</span> T_C")
        legend.setStyleSheet(f"color:{MUTED}")
        layout.addWidget(legend)

        self.mechanism = _Section("Mechanism", (
            "Piston energy", "Shaft inertia energy", "Spring energy"),
            collapsible=True)
        layout.addWidget(self.mechanism)
        self.diagnostics = _Section("Numerical diagnostics", (
            "Max penetration", "CCD refinements / failures", "Contact clusters",
            "Max cluster residual", "Last cycle gas residual", "P dA − ledger work",
            "Efficiency cycles / storage ΔE", "Event throughput",
            "Achieved playback"), collapsible=True)
        layout.addWidget(self.diagnostics)
        self.diagnostics_toggle = self.diagnostics.toggle
        self.diagnostic_values = self.diagnostics.values
        tooltip = self.diagnostic_values
        tooltip["Max penetration"].setToolTip(
            "Largest overlap measured during this run, in simulation length units.")
        tooltip["Event throughput"].setToolTip(
            "Processed physics events per elapsed wall-clock second during playback.")
        tooltip["Achieved playback"].setToolTip(
            "Physical simulation seconds advanced per elapsed wall-clock second.")

    def show_inspection(self, text):
        self.inspector.setText(text)
        self.inspector.setVisible(bool(text))

    def set_rates(self, rates, replay=False):
        if replay:
            text = ("Recorded", "Replay")
        elif rates:
            text = (f"{rates.events_per_wall:.1f} events/s",
                    f"{rates.physical_per_wall:.3f}×")
        else:
            text = ("Paused", "Paused")
        self.diagnostic_values["Event throughput"].setText(text[0])
        self.diagnostic_values["Achieved playback"].setText(text[1])

    def set_temperatures(self, series: InstrumentSeries, hot, cold):
        self.temp_trans.setData(series.time, series.temperature_trans)
        self.temp_rot.setData(series.time, series.temperature_rot)
        flat = np.ones_like(series.time)
        self.temp_hot.setData(series.time, flat*hot if hot is not None else flat*np.nan)
        self.temp_cold.setData(series.time, flat*cold if cold is not None else flat*np.nan)

    def set_sample(self, current, phase_fraction=None, ideal_efficiency=None):
        cycle = current.latest_cycle
        pressure_area = current.latest_pressure_area
        efficiency = current.efficiency_report
        estimate = efficiency.estimate
        branch = BRANCH_LABELS.get(current.branch, current.branch or "—")
        if phase_fraction is not None:
            branch += f" · {100*phase_fraction:.0f}%"
        values = {
            "Branch / phase": branch,
            "Time / cycles": f"{current.time:.2f} s / {current.completed_cycles}",
            "Area": fmt(current.area),
            "Pressure (window)": f"{fmt(current.pressure)} ({current.pressure_window:.2f} s)",
            "T translation / rotation":
                f"{fmt(current.temperature_trans)} / {fmt(current.temperature_rot)}",
            "Reservoirs T_H / T_C":
                f"{fmt(current.reservoir_hot)} / {fmt(current.reservoir_cold)}",
            "Gas energy": fmt(current.gas_energy),
            "Heat in from hot Q_H": fmt(current.heat_hot),
            "Heat in from cold Q_C": fmt(current.heat_cold),
            "Motor work / load output":
                f"{fmt(current.motor_work)} / {fmt(current.load_output)}",
            "First-law residual": f"{current.first_law_residual:.2e}",
            "Q_H / Q_C": (f"{fmt(cycle.heat_hot)} / {fmt(cycle.heat_cold)}"
                          if cycle else "—"),
            "Work by gas: ledger / ∮P dA": (
                f"{fmt(-cycle.piston_work_on_gas)} / "
                + (f"{fmt(pressure_area.integrated_work_by_gas)} "
                   f"± {fmt(pressure_area.smoothing_error_bound, 2)}"
                   if pressure_area else "—") if cycle else "—"),
            "External output": fmt(cycle.external_output) if cycle else "—",
            "Total residual": f"{cycle.total_first_law_residual:.2e}" if cycle else "—",
            "η_net ± SE / η_Carnot": (
                (f"{estimate.value:.3g} ± {estimate.block_standard_error:.2g}"
                 if estimate else "—") + f" / {fmt(ideal_efficiency, 3)}"),
            "Gate": (f"ready · {efficiency.transient_cycles} warm-up skipped"
                     if estimate else
                     f"{efficiency.status.replace('_', ' ')} · "
                     f"{efficiency.matched_cycles}/{efficiency.required_cycles}"),
        }
        for name, value in values.items():
            self.values[name].setText(value)
        self.mechanism.values["Piston energy"].setText(fmt(current.piston_energy))
        self.mechanism.values["Shaft inertia energy"].setText(fmt(current.flywheel_energy))
        self.mechanism.values["Spring energy"].setText(fmt(current.spring_energy))
        diagnostics = {
            "Max penetration": f"{current.max_penetration:.2e}",
            "CCD refinements / failures":
                f"{current.ccd_refinements} / {current.ccd_failures}",
            "Contact clusters": str(current.cluster_count),
            "Max cluster residual": f"{current.max_cluster_residual:.2e}",
            "Last cycle gas residual":
                f"{cycle.gas_first_law_residual:.2e}" if cycle else "—",
            "P dA − ledger work": fmt(pressure_area.difference) if pressure_area else "—",
            "Efficiency cycles / storage ΔE": (
                f"{estimate.sample_count} / {estimate.storage_change:.3g}"
                if estimate else "—"),
        }
        for name, value in diagnostics.items():
            self.diagnostic_values[name].setText(value)
