"""The laboratory layout shared by the live window and the replay window.

Left: the apparatus, then the p–V and T–S diagrams and the cycle scoreboard
(Carnot scenes, with a cycle dial) or the osmosis plots (osmosis scenes).
Right: the readout column.
"""

from __future__ import annotations

import math
import time

from PySide6 import QtCore, QtGui, QtWidgets

from .instruments import (BRANCH_COLORS, BRANCH_LABELS, CycleDiagrams,
                          InstrumentSeries, ReadoutPanel)
from .main_window import ApparatusView
from .osmosis_panel import OsmosisPanel


class CycleDial(QtWidgets.QWidget):
    """Shaft phase on a ring of the four cam sectors, clockwise from 12 o'clock."""

    def __init__(self, sim):
        super().__init__()
        mechanism = sim.world.mechanism
        cam = mechanism.cam
        reverse = bool(getattr(mechanism, "reversed_cycle", False))
        bounds = [float(b) for b in cam.boundaries]
        self.sectors = tuple((start, end, cam.branch((start+end)/2, reverse))
                             for start, end in zip(bounds[:-1], bounds[1:]))
        self.direction = -1. if reverse else 1.
        shaft = mechanism.shaft
        self.mode = (f"motor {abs(shaft.prescribed_omega):g} rad/s"
                     if shaft.prescribed_omega is not None else "free shaft")
        self.phase = None
        self.branch = None
        self.setFixedWidth(150)
        self.setMinimumHeight(150)

    def set_snapshot(self, snapshot):
        self.phase = snapshot.shaft_phase
        self.branch = snapshot.branch
        self.update()

    def phase_fraction(self):
        """Fraction of the current cam sector already travelled."""
        if self.phase is None:
            return None
        local = (self.direction*self.phase) % (2*math.pi)
        for start, end, _ in self.sectors:
            if start <= local < end:
                return (local-start)/(end-start)
        return 1.

    def paintEvent(self, event):
        painter = QtGui.QPainter(self)
        painter.fillRect(self.rect(), QtGui.QColor("#15191f"))
        painter.setRenderHint(QtGui.QPainter.RenderHint.Antialiasing)
        radius = max(20., min(self.width()/2-12, (self.height()-48)/2))
        center = QtCore.QPointF(self.width()/2, 8+radius)
        box = QtCore.QRectF(center.x()-radius, center.y()-radius, 2*radius, 2*radius)
        painter.setPen(QtCore.Qt.PenStyle.NoPen)
        # Qt pie angles run counter-clockwise from 3 o'clock in 1/16 degree.
        for start, end, branch in self.sectors:
            color = QtGui.QColor(BRANCH_COLORS.get(branch, "#8997a7"))
            color.setAlpha(255 if branch == self.branch else 90)
            painter.setBrush(color)
            painter.drawPie(box, int((90-math.degrees(start))*16),
                            int(-math.degrees(end-start)*16))
        painter.setBrush(QtGui.QColor("#15191f"))
        painter.drawEllipse(center, radius*.58, radius*.58)
        if self.phase is not None:
            angle = (self.direction*self.phase) % (2*math.pi)
            tip = center+QtCore.QPointF(radius*math.sin(angle), -radius*math.cos(angle))
            painter.setPen(QtGui.QPen(QtGui.QColor("#ffffff"), 2.5))
            painter.drawLine(center, tip)
            fraction = self.phase_fraction()
            painter.setPen(QtGui.QColor("#d5dce4"))
            painter.drawText(QtCore.QRectF(center.x()-radius*.55, center.y()-10,
                                           radius*1.1, 20),
                             QtCore.Qt.AlignmentFlag.AlignCenter, f"{100*fraction:.0f}%")
        painter.setPen(QtGui.QColor("#d5dce4"))
        label = BRANCH_LABELS.get(self.branch, self.branch or "")
        painter.drawText(QtCore.QRectF(0, center.y()+radius+4, self.width(), 18),
                         QtCore.Qt.AlignmentFlag.AlignCenter, label)
        painter.setPen(QtGui.QColor("#8e99a8"))
        painter.drawText(QtCore.QRectF(0, center.y()+radius+22, self.width(), 18),
                         QtCore.Qt.AlignmentFlag.AlignCenter, self.mode)
        painter.end()


class LabView(QtWidgets.QSplitter):
    # Plots redraw at most this often during playback; the particles, dial
    # and numbers update every frame. Redrawing several thousand plot points
    # costs tens of milliseconds and would otherwise stall the animation.
    PLOT_INTERVAL = .15

    def __init__(self, sim):
        super().__init__(QtCore.Qt.Orientation.Horizontal)
        self.left = QtWidgets.QSplitter(QtCore.Qt.Orientation.Vertical)
        self.top = QtWidgets.QWidget()
        self.top_layout = QtWidgets.QHBoxLayout(self.top)
        self.top_layout.setContentsMargins(0, 0, 0, 0)
        self.top_layout.setSpacing(0)
        self.diagrams = CycleDiagrams()
        self.osmosis = OsmosisPanel()
        self.readout = ReadoutPanel()
        self.left.addWidget(self.top)
        self.left.addWidget(self.diagrams)
        self.left.addWidget(self.osmosis)
        self.left.setStretchFactor(0, 2)
        self.left.setStretchFactor(1, 3)
        self.addWidget(self.left)
        self.addWidget(self.readout)
        self.setStretchFactor(0, 3)
        self.setStretchFactor(1, 1)
        self.view = None
        self.dial = None
        self._plotted_at = -math.inf
        self.set_simulation(sim)

    def set_simulation(self, sim):
        """Install the apparatus for a new scene; diagrams appear for Carnot."""
        for widget in (self.view, self.dial):
            if widget is not None:
                self.top_layout.removeWidget(widget)
                widget.deleteLater()
        self.view = ApparatusView(sim)
        self.view.inspected.connect(self.readout.show_inspection)
        self.top_layout.addWidget(self.view, stretch=1)
        carnot = self.view.cylinder_only
        self.dial = CycleDial(sim) if carnot else None
        if self.dial is not None:
            self.top_layout.addWidget(self.dial)
        self.diagrams.setVisible(carnot)
        self.diagrams.set_reference(None)
        osmosis = bool(sim.world.rings) and "membrane_x" in sim.world.metadata
        self.osmosis.setVisible(osmosis)
        self.osmosis.clear()
        self._discs = int(sum(1 for i in range(sim.world.bodies.n) if i not in sim.world.rings))
        # Cycle, efficiency and mechanism readouts belong to Carnot scenes.
        for section in self.readout.sections[2:] + [self.readout.mechanism]:
            section.setVisible(not osmosis)
        state, cumulative = self.readout.sections[:2]
        for name in ("Branch / phase", "Area", "Pressure (window)", "Reservoirs T_H / T_C"):
            state.set_row(name, not osmosis)
        for name in ("Heat in from cold Q_C", "Motor work / load output"):
            cumulative.set_row(name, not osmosis)
        cumulative.set_row("Heat in from hot Q_H", label="Heat in from the walls"
                           if osmosis else None)
        state.set_row("Gas energy", label="Kinetic energy" if osmosis else None)
        self.readout.show_inspection("")
        if carnot:
            self.left.setSizes([260, 520])
        elif osmosis:
            self.left.setSizes([520, 0, 380])

    def plots_due(self, force=False):
        """Whether the next frame should redraw the plots."""
        return force or time.perf_counter()-self._plotted_at >= self.PLOT_INTERVAL

    def show_frame(self, snapshot, frame, series: InstrumentSeries | None,
                   rates=None, replay=False, osmosis=None):
        """Update everything; ``series`` None skips the plots this frame."""
        if self.osmosis.isVisible() or osmosis is not None:
            self.osmosis.update_frame(osmosis, self._discs, plots=series is not None)
        current = frame.current
        self.view.set_snapshot(snapshot)
        fraction = None
        if self.dial is not None:
            self.dial.set_snapshot(snapshot)
            fraction = self.dial.phase_fraction()
        hot, cold = current.reservoir_hot, current.reservoir_cold
        ideal = 1-cold/hot if hot and cold else None
        self.readout.set_sample(current, fraction, ideal)
        self.readout.set_rates(rates, replay=replay)
        if self.dial is not None:
            self.diagrams.update_scoreboard(
                current, frame.cycles, frame.cycle_start,
                current.efficiency_report.transient_cycles, ideal)
        if series is None:
            return
        self._plotted_at = time.perf_counter()
        self.readout.set_temperatures(series, hot, cold)
        if self.dial is not None:
            self.diagrams.set_reference(frame.ideal_reference)
            # Colour the previous and the current cycle; dim older ones.
            recent_from = (frame.cycles[-1].start_time if frame.cycles else -math.inf)
            self.diagrams.update_series(series, current, frame.particles,
                                        frame.degrees_of_freedom, recent_from)
