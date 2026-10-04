#!/usr/bin/env bash
# TEMPORARY: delete this launcher and TEMP_gap_desktop_controls_batch_*
# after review, retaining only compact results and an accepted screenshot.
# Run from the graphical desktop session:
# bash runs/TEMP_gap_desktop_controls_acceptance.sh

repo=/home/blue-whale/randomArticles/carnot
cd "$repo" || exit 1
batch=$(mktemp -d "$repo/runs/TEMP_gap_desktop_controls_batch_$(date +%Y%m%d_%H%M%S)_XXXXXX") || exit 1
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
printf 'Results: %s\n' "$batch"

(
  QT_QPA_PLATFORM=offscreen /usr/bin/time -v .venv/bin/python -m unittest discover -s tests -q
  printf '%s\n' "$?" > "$batch/test_suite.exit"
) > "$batch/test_suite.stdout" 2> "$batch/test_suite.stderr" &
(
  /usr/bin/time -v .venv/bin/python -m microthermo validate --suite scientific
  printf '%s\n' "$?" > "$batch/scientific_validation.exit"
) > "$batch/scientific_validation.stdout" 2> "$batch/scientific_validation.stderr" &
wait

(
  GAP_BATCH="$batch" /usr/bin/time -v .venv/bin/python - <<'PY'
import json, math, os, time
from pathlib import Path
import numpy as np
from PySide6 import QtCore, QtTest, QtWidgets
from microthermo.api import load_preset
from microthermo.config import RunConfig
from microthermo.ui.main_window import MainWindow

app=QtWidgets.QApplication([])
platform=app.platformName().lower()
if platform in ('offscreen','minimal'):
    raise RuntimeError(f'real desktop Qt required; got {platform}')
window=MainWindow('carnot_discs',autoplay=False)
window.show()
def until(condition,timeout=30):
    deadline=time.perf_counter()+timeout
    while time.perf_counter()<deadline:
        app.processEvents()
        if condition(): return
        time.sleep(.005)
    raise RuntimeError('timed out waiting for desktop command')
until(lambda: window.instruments.values['Area'].text() != '—')
config=RunConfig(preset='carnot_discs',particles=48,seed=123,
                 max_horizon=.02,shaft_speed=.15)
reference=load_preset(config)

window.playback_rate.setCurrentIndex(window.playback_rate.findData(4.))
until(lambda: abs(window.worker._physical_step-.064)<1e-12)
window.duration_step.setValue(.05)
window.step_duration()
until(lambda: window.view.snapshot.time >= .05)
reference.advance_to(.05)
np.testing.assert_array_equal(window.view.snapshot.position,
                              reference.snapshot().position)
window.shaft_speed.setValue(.3)
window.change_shaft_speed()
until(lambda: window.view.snapshot.event_count > reference.event_count)
reference.apply_command({'name':'set_shaft_speed','speed':.3})
window.duration_step.setValue(.08)
window.step_duration()
until(lambda: window.view.snapshot.time >= .13-1e-12)
reference.advance_to(.13)
np.testing.assert_array_equal(window.view.snapshot.position,
                              reference.snapshot().position)

pings=[]
timer=QtCore.QTimer()
timer.setInterval(20)
timer.timeout.connect(lambda: pings.append(time.perf_counter()))
timer.start()
start=time.perf_counter()
window.step_branch()
until(lambda: window.view.snapshot.branch == 'adiabatic_expansion',timeout=60)
branch_wall_seconds=time.perf_counter()-start
timer.stop()
reference.advance_to(.05+(math.pi/2-.05*.15)/.3)
actual=window.view.snapshot
expected=reference.snapshot()
np.testing.assert_array_equal(actual.position,expected.position)
np.testing.assert_array_equal(actual.velocity,expected.velocity)
assert actual.event_count==expected.event_count
assert abs(actual.energy_residual-expected.energy_residual)<1e-10
assert len(pings)>0
assert window.status.currentMessage().startswith('PAUSED')

view=window.view
before_pan=(view.pan_x,view.pan_y)
QtTest.QTest.mousePress(view,QtCore.Qt.MouseButton.LeftButton,
                        pos=QtCore.QPoint(8,8))
QtTest.QTest.mouseMove(view,QtCore.QPoint(38,28))
QtTest.QTest.mouseRelease(view,QtCore.Qt.MouseButton.LeftButton,
                          pos=QtCore.QPoint(38,28))
assert (view.pan_x,view.pan_y)!=before_pan
view.reset_camera()
point=view.camera().map(*actual.position[0])
QtTest.QTest.mouseClick(view,QtCore.Qt.MouseButton.LeftButton,
                        pos=QtCore.QPoint(round(point[0]),round(point[1])))
assert view.selected_particle==0

window.preset_choice.setCurrentText('gas_box')
until(lambda: window._generation==1 and window.worker._generation==1 and
      window.view.snapshot.time==0)
assert not window.shaft_speed.isEnabled()
window.seed.setValue(124)
window.reset_to_seed()
gas=load_preset(RunConfig(preset='gas_box',particles=48,seed=124,
                           max_horizon=.02)).snapshot()
until(lambda: window._generation==2 and window.worker._generation==2 and
      window.view.snapshot.time==0 and
      np.array_equal(window.view.snapshot.position,gas.position))
np.testing.assert_array_equal(window.view.snapshot.position,gas.position)
window.resize(850,520); app.processEvents()
assert window.view.width()>=500 and window.instruments.width()>=320
shot=Path(os.environ['GAP_BATCH'])/'desktop_controls.png'
assert window.grab().save(str(shot))
result=dict(qt_platform=platform,branch_step_wall_seconds=branch_wall_seconds,
            gui_pings_during_branch_step=len(pings),
            matched_branch_time=actual.time,matched_event_count=actual.event_count,
            matched_energy_residual=actual.energy_residual,
            panned=True,selected_particle=0,preset_switched=True,
            reset_seed=124,minimum_window=[window.width(),window.height()],
            screenshot=str(shot))
print(json.dumps(result,indent=2))
window.close()
PY
  printf '%s\n' "$?" > "$batch/desktop.exit"
) > "$batch/desktop.stdout" 2> "$batch/desktop.stderr"

printf '\nExit codes (0 means success):\n'
for file in "$batch"/*.exit; do printf '%s: ' "${file##*/}"; cat "$file"; done
printf '\nResults: %s\n' "$batch"
