"""Live Carnot plots and event-ledger readouts."""

from __future__ import annotations

import math
from PySide6 import QtCore, QtWidgets
import pyqtgraph as pg


BRANCH_COLORS = {
    "hot":"#f2a65a", "adiabatic_expansion":"#c9bddf",
    "cold":"#76bce3", "adiabatic_compression":"#c9bddf",
}


class InstrumentPanel(QtWidgets.QWidget):
    def __init__(self):
        super().__init__()
        self.setMinimumWidth(320)
        self.setStyleSheet("background:#15191f;color:#d5dce4")
        outer=QtWidgets.QVBoxLayout(self)
        outer.setContentsMargins(4,4,4,4)
        scroll=QtWidgets.QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QtWidgets.QFrame.Shape.NoFrame)
        outer.addWidget(scroll)
        content=QtWidgets.QWidget()
        layout=QtWidgets.QVBoxLayout(content)
        layout.setContentsMargins(4,4,4,4)
        scroll.setWidget(content)

        heading=QtWidgets.QLabel("LIVE INSTRUMENTS")
        heading.setStyleSheet("font-weight:bold;color:#f0c46c")
        layout.addWidget(heading)
        self.summary=QtWidgets.QLabel("Waiting for worker snapshot")
        self.summary.setWordWrap(True)
        self.summary.setStyleSheet("color:#d5dce4")
        layout.addWidget(self.summary)
        self.pa=pg.PlotWidget(background="#15191f")
        self.pa.setFixedHeight(160)
        self.pa.setLabel("bottom","Area A")
        self.pa.setLabel("left","Pressure P")
        self.pa.showGrid(x=True,y=True,alpha=.2)
        self.pa_reference={name:self.pa.plot(
            pen=pg.mkPen(color,width=2,style=QtCore.Qt.PenStyle.DashLine))
            for name,color in BRANCH_COLORS.items()}
        self._shown_reference=None
        self.pa_line=self.pa.plot(pen=pg.mkPen("#697481",width=1))
        self.pa_points=pg.ScatterPlotItem(size=5,pen=None)
        self.pa.addItem(self.pa_points)
        self.pa_current=pg.ScatterPlotItem(size=10,
                                          pen=pg.mkPen("#f3e7c3",width=1),
                                          brush=pg.mkBrush("#f3e7c3"))
        self.pa.addItem(self.pa_current)
        layout.addWidget(self.pa)
        self.pa_legend=QtWidgets.QLabel("● measured impulse window   ┄ analytical point gas")
        self.pa_legend.setWordWrap(True)
        self.pa_legend.setStyleSheet("color:#aeb9c8")
        layout.addWidget(self.pa_legend)

        self.temperatures=pg.PlotWidget(background="#15191f")
        self.temperatures.setFixedHeight(160)
        self.temperatures.setLabel("bottom","Physical time")
        self.temperatures.setLabel("left","Temperature")
        self.temperatures.showGrid(x=True,y=True,alpha=.2)
        self.temp_trans=self.temperatures.plot(pen=pg.mkPen("#f0c46c",width=2))
        self.temp_rot=self.temperatures.plot(pen=pg.mkPen("#ab9edb",width=2))
        self.temp_hot=self.temperatures.plot(pen=pg.mkPen("#f2a65a",width=1,
                                           style=QtCore.Qt.PenStyle.DashLine))
        self.temp_cold=self.temperatures.plot(pen=pg.mkPen("#76bce3",width=1,
                                            style=QtCore.Qt.PenStyle.DashLine))
        layout.addWidget(self.temperatures)
        legend=QtWidgets.QLabel("● translation   ● rotation   ┄ hot   ┄ cold")
        legend.setStyleSheet("color:#aeb9c8")
        layout.addWidget(legend)

        grid=QtWidgets.QGridLayout()
        layout.addLayout(grid)
        rows=("Branch","Area","Pressure / window","T translation","T rotation",
              "Reservoir hot / cold","Gas energy","Piston energy",
              "Shaft inertia energy","Spring energy","Hot heat QH","Cold heat QC",
              "Motor work","Load output","First-law residual","Events",
              "Completed cycles","Last cycle QH / QC","Last cycle gas work",
              "Cycle ∮P dA / smoothing bound","P dA − event work",
              "Last cycle external output","Last cycle gas residual",
              "Last cycle total residual","Efficiency gate",
              "Measured ηnet / block SE","Efficiency cycles / storage ΔE")
        self.values={}
        for row,name in enumerate(rows):
            label=QtWidgets.QLabel(name)
            value=QtWidgets.QLabel("—")
            if name == "Cycle ∮P dA / smoothing bound":
                label.setToolTip("Work by gas from the fixed trailing pressure window; "
                                 "bound covers window smoothing and cycle clipping, "
                                 "not statistical sampling uncertainty.")
            if name == "Measured ηnet / block SE":
                label.setToolTip("Signed net external output divided by hot heat "
                                 "over matched completed cycles. SE is a four-block "
                                 "jackknife standard error, not a confidence interval.")
            value.setAlignment(QtCore.Qt.AlignmentFlag.AlignRight)
            value.setTextInteractionFlags(
                QtCore.Qt.TextInteractionFlag.TextSelectableByMouse)
            grid.addWidget(label,row,0)
            grid.addWidget(value,row,1)
            self.values[name]=value

        self.diagnostics_toggle=QtWidgets.QToolButton()
        self.diagnostics_toggle.setText("Numerical diagnostics")
        self.diagnostics_toggle.setCheckable(True)
        self.diagnostics_toggle.setChecked(False)
        self.diagnostics_toggle.setToolButtonStyle(
            QtCore.Qt.ToolButtonStyle.ToolButtonTextBesideIcon)
        self.diagnostics_toggle.setArrowType(QtCore.Qt.ArrowType.RightArrow)
        layout.addWidget(self.diagnostics_toggle)
        self.diagnostics=QtWidgets.QWidget()
        diagnostic_grid=QtWidgets.QGridLayout(self.diagnostics)
        diagnostic_grid.setContentsMargins(8,0,0,0)
        diagnostic_names=("Max penetration", "CCD refinements / failures",
                          "Contact clusters", "Max cluster residual",
                          "Event throughput", "Achieved playback")
        self.diagnostic_values={}
        for row,name in enumerate(diagnostic_names):
            label=QtWidgets.QLabel(name)
            value=QtWidgets.QLabel("—")
            value.setAlignment(QtCore.Qt.AlignmentFlag.AlignRight)
            value.setTextInteractionFlags(
                QtCore.Qt.TextInteractionFlag.TextSelectableByMouse)
            diagnostic_grid.addWidget(label,row,0)
            diagnostic_grid.addWidget(value,row,1)
            self.diagnostic_values[name]=value
        self.diagnostic_values["Max penetration"].setToolTip(
            "Largest overlap measured during this run, in simulation length units.")
        self.diagnostic_values["Max cluster residual"].setToolTip(
            "Largest energy residual of a simultaneous elastic contact solve.")
        self.diagnostic_values["Event throughput"].setToolTip(
            "Processed physics events per elapsed wall-clock second during playback.")
        self.diagnostic_values["Achieved playback"].setToolTip(
            "Physical simulation seconds advanced per elapsed wall-clock second.")
        self.diagnostics.setVisible(False)
        layout.addWidget(self.diagnostics)
        self.diagnostics_toggle.toggled.connect(self._toggle_diagnostics)
        layout.addStretch(1)

    def _toggle_diagnostics(self, expanded):
        self.diagnostics.setVisible(expanded)
        self.diagnostics_toggle.setArrowType(
            QtCore.Qt.ArrowType.DownArrow if expanded else
            QtCore.Qt.ArrowType.RightArrow)

    def set_rates(self, rates):
        self.diagnostic_values["Event throughput"].setText(
            f"{rates.events_per_wall:.1f} events/s" if rates else "Paused")
        self.diagnostic_values["Achieved playback"].setText(
            f"{rates.physical_per_wall:.3f}×" if rates else "Paused")

    def set_frame(self, frame):
        current=frame.current
        history=frame.history
        reference=frame.ideal_reference
        if reference is not self._shown_reference:
            self._shown_reference=reference
            branches={branch.name:branch for branch in reference.branches} if reference else {}
            for name,line in self.pa_reference.items():
                branch=branches.get(name)
                line.setData(branch.area,branch.pressure) if branch else line.setData([],[])
            self.pa_legend.setText(
                f"● measured impulse window   ┄ analytical point gas "
                f"(N={reference.particles}, f={reference.degrees_of_freedom}, "
                f"ηC={reference.efficiency:.3f})" if reference else
                "● measured impulse window")
        pa=[sample for sample in history if sample.area is not None and
            sample.pressure is not None and math.isfinite(sample.pressure)]
        self.pa_line.setData([sample.area for sample in pa],
                             [sample.pressure for sample in pa])
        self.pa_points.setData(x=[sample.area for sample in pa],
                               y=[sample.pressure for sample in pa],
                               brush=[pg.mkBrush(BRANCH_COLORS.get(sample.branch,
                                                                    "#d5dce4"))
                                      for sample in pa])
        if current.area is not None and current.pressure is not None:
            self.pa_current.setData(x=[current.area],y=[current.pressure])
        else:
            self.pa_current.setData([],[])
        self.pa.setTitle(f"Pressure–area · Δt={current.pressure_window:.3f} s"
                         if current.area is not None else "Pressure–area · Carnot")
        times=[sample.time for sample in history]
        self.temp_trans.setData(times,[sample.temperature_trans for sample in history])
        self.temp_rot.setData(times,[sample.temperature_rot for sample in history])
        self.temp_hot.setData(times,[sample.reservoir_hot if sample.reservoir_hot
                                     is not None else math.nan for sample in history])
        self.temp_cold.setData(times,[sample.reservoir_cold if sample.reservoir_cold
                                      is not None else math.nan for sample in history])
        self.temperatures.setTitle("Temperature · translation / rotation")

        def fmt(value):
            if value is None or not math.isfinite(value):
                return "—"
            return f"{value:.4g}"
        self.summary.setText(
            f"A {fmt(current.area)}   P {fmt(current.pressure)}   "
            f"T {fmt(current.temperature_trans)}\n"
            f"QH {fmt(current.heat_hot)}   QC {fmt(current.heat_cold)}   "
            f"Wload {fmt(current.load_output)}   Cycles {current.completed_cycles}")
        cycle=current.latest_cycle
        pressure_area=current.latest_pressure_area
        efficiency=current.efficiency_report
        estimate=efficiency.estimate
        values={
            "Branch":current.branch.replace("_"," ") if current.branch else "—",
            "Area":fmt(current.area),
            "Pressure / window":f"{fmt(current.pressure)} / {current.pressure_window:.3f} s",
            "T translation":fmt(current.temperature_trans),
            "T rotation":fmt(current.temperature_rot),
            "Reservoir hot / cold":f"{fmt(current.reservoir_hot)} / {fmt(current.reservoir_cold)}",
            "Gas energy":fmt(current.gas_energy),
            "Piston energy":fmt(current.piston_energy),
            "Shaft inertia energy":fmt(current.flywheel_energy),
            "Spring energy":fmt(current.spring_energy),
            "Hot heat QH":fmt(current.heat_hot),
            "Cold heat QC":fmt(current.heat_cold),
            "Motor work":fmt(current.motor_work),
            "Load output":fmt(current.load_output),
            "First-law residual":f"{current.first_law_residual:.2e}",
            "Events":str(current.event_count),
            "Completed cycles":str(current.completed_cycles),
            "Last cycle QH / QC":(
                f"{fmt(cycle.heat_hot)} / {fmt(cycle.heat_cold)}" if cycle else "—"),
            "Last cycle gas work":fmt(cycle.piston_work_on_gas) if cycle else "—",
            "Cycle ∮P dA / smoothing bound":(
                f"{fmt(pressure_area.integrated_work_by_gas)} / "
                f"±{fmt(pressure_area.smoothing_error_bound)}"
                if pressure_area else "—"),
            "P dA − event work":fmt(pressure_area.difference) if pressure_area else "—",
            "Last cycle external output":fmt(cycle.external_output) if cycle else "—",
            "Last cycle gas residual":(
                f"{cycle.gas_first_law_residual:.2e}" if cycle else "—"),
            "Last cycle total residual":(
                f"{cycle.total_first_law_residual:.2e}" if cycle else "—"),
            "Efficiency gate":(
                f"ready · skip {efficiency.transient_cycles}" if estimate else
                f"{efficiency.status.replace('_',' ')} · "
                f"{efficiency.matched_cycles}/{efficiency.required_cycles}"),
            "Measured ηnet / block SE":(
                f"{estimate.value:.4g} / {estimate.block_standard_error:.2g}"
                if estimate else "—"),
            "Efficiency cycles / storage ΔE":(
                f"{estimate.sample_count} / {estimate.storage_change:.3g}"
                if estimate else "—"),
        }
        for name,value in values.items():
            self.values[name].setText(value)
        diagnostics={
            "Max penetration":f"{current.max_penetration:.2e}",
            "CCD refinements / failures":(
                f"{current.ccd_refinements} / {current.ccd_failures}"),
            "Contact clusters":str(current.cluster_count),
            "Max cluster residual":f"{current.max_cluster_residual:.2e}",
        }
        for name,value in diagnostics.items():
            self.diagnostic_values[name].setText(value)
