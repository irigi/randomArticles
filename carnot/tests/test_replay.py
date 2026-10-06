"""Chunked replay archive round trip and bounded-memory seek."""

import json
import math
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
from types import SimpleNamespace

import numpy as np

from microthermo.api import load_preset
from microthermo.cli import main
from microthermo.config import RunConfig
from microthermo.core.apparatus import ApparatusComponent
from microthermo.io.exports import _json
from microthermo.io.replay import ReplayReader, _interpolate_apparatus, precalculate_replay
from microthermo.measurements.live import LiveInstruments


class ReplayTests(unittest.TestCase):
    def test_selector_shoe_switches_at_branch_event(self):
        old = SimpleNamespace(apparatus=(
            ApparatusComponent("selector_shoe", "point", ((.22, 0.),), state="hot"),
            ApparatusComponent("piston", "point", ((1., 0.),)),
        ))
        new = SimpleNamespace(apparatus=(
            ApparatusComponent("selector_shoe", "point", ((.19, 0.),),
                               state="adiabatic_expansion"),
            ApparatusComponent("piston", "point", ((2., 0.),)),
        ))
        before = _interpolate_apparatus(old, new, .9, False)
        after = _interpolate_apparatus(old, new, .9, True)
        self.assertEqual(before[0].points, ((.22, 0.),))
        self.assertEqual(before[0].state, "hot")
        self.assertEqual(after[0].points, ((.19, 0.),))
        self.assertEqual(after[0].state, "adiabatic_expansion")
        self.assertAlmostEqual(after[1].points[0][0], 1.9)

    def test_draining_events_preserves_future_history_and_state(self):
        config=RunConfig(preset="carnot_discs",particles=16,seed=123,
                         duration=.5)
        streaming=load_preset(config)
        reference=load_preset(config)
        streaming.advance_to(.25)
        first=streaming.drain_events()
        self.assertEqual(len(first),streaming.event_count)
        self.assertEqual(streaming.events,[])
        later=streaming.advance_to(.5)
        expected=reference.advance_to(.5)
        self.assertEqual(first+streaming.events,reference.events)
        np.testing.assert_array_equal(later.position,expected.position)
        self.assertEqual(later.event_count,expected.event_count)

    def test_recorded_instruments_match_live_observer(self):
        # A fast shaft completes a cycle, so the drained-event observer must
        # reproduce the live pressure window and the cycle P dA integral.
        config=RunConfig(preset="carnot_discs",particles=8,seed=123,
                         shaft_speed=3.,duration=2.4)
        with tempfile.TemporaryDirectory() as directory:
            reader=ReplayReader(precalculate_replay(
                config,Path(directory)/"replay",fps=20.,chunk_frames=8))
            self.assertTrue(reader.has_instruments)
            direct=load_preset(config)
            live=LiveInstruments()
            for index in range(len(reader)):
                saved=reader.instrument_sample(index)
                snapshot=(direct.snapshot() if index==0 else
                          direct.advance_to(saved.time))
                expected=live.observe(direct,snapshot).current
                for name in ("time","area","pressure","entropy","gas_energy",
                             "heat_hot","heat_cold","motor_work",
                             "completed_cycles","branch"):
                    self.assertEqual(getattr(saved,name),getattr(expected,name),name)
            self.assertEqual(saved.completed_cycles,1)
            self.assertEqual(saved.latest_cycle,expected.latest_cycle)
            self.assertEqual(saved.latest_pressure_area,expected.latest_pressure_area)
            self.assertEqual(saved.efficiency_report,expected.efficiency_report)
            series=reader.instrument_series()
            self.assertEqual(len(series["time"]),len(reader))
            reader.close()

    def test_chunked_frames_and_events_match_source_run(self):
        config=RunConfig(preset="carnot_discs",particles=4,seed=123,
                         duration=.5)
        with tempfile.TemporaryDirectory() as directory:
            path=precalculate_replay(config,Path(directory)/"replay",fps=10.,
                                     chunk_frames=2)
            reader=ReplayReader(path)
            manifest=reader.manifest
            self.assertEqual(manifest["version"],3)
            # Carnot archives carry no osmosis columns or membership.
            self.assertNotIn("osmosis",manifest)
            self.assertFalse(reader.has_osmosis)
            self.assertIsNone(reader.frame(0).membership)
            self.assertEqual(manifest["config_hash"],config.digest())
            self.assertEqual(len(reader),6)
            self.assertEqual(len(manifest["chunks"]),3)
            self.assertEqual(manifest["geometry"]["shape"],[0]*4)
            direct=load_preset(config)
            for index in range(len(reader)):
                target=index/10
                source=direct.snapshot() if index==0 else direct.advance_to(target)
                saved=reader.frame(index)
                self.assertEqual(saved.time,source.time)
                np.testing.assert_array_equal(saved.position,source.position)
                np.testing.assert_array_equal(saved.velocity,source.velocity)
                np.testing.assert_array_equal(saved.angle,source.angle)
                np.testing.assert_array_equal(saved.omega,source.omega)
                self.assertEqual(saved.event_count,source.event_count)
                self.assertEqual(saved.energy,source.energy)
                self.assertEqual(saved.energy_residual,source.energy_residual)
                self.assertEqual(saved.apparatus,source.apparatus)
                self.assertEqual(reader.ledger(index)["heat_hot"],
                                 direct.ledger.heat_hot.value)
                self.assertEqual(reader.ledger(index)["work_on"],
                                 direct.ledger.work_on.value)
            events=[json.loads(line) for line in (path/"events.jsonl").read_text().splitlines()]
            self.assertEqual(events,[json.loads(json.dumps(event.as_dict(),
                                                       default=_json))
                                     for event in direct.events])
            self.assertEqual(manifest["event_count"],direct.event_count)
            self.assertEqual(manifest["final_ledger"]["heat_hot"]["value"],
                             direct.ledger.heat_hot.value)
            reader.prepare_playback()
            with patch("microthermo.io.replay.np.load",
                       side_effect=AssertionError("display seek read archive")):
                for time in (.05,.27,.45,.15):
                    reader.sample(time)
            self.assertEqual(reader.seek(.27).time,.2)
            interpolated, ledger_index=reader.sample(.27)
            self.assertEqual(ledger_index,2)
            self.assertEqual(interpolated.time,.27)
            self.assertEqual(interpolated.event_count,reader.frame(2).event_count)
            midpoint_source=load_preset(config).advance_to(.27)
            np.testing.assert_allclose(
                interpolated.position,midpoint_source.position,
                atol=1e-10,rtol=0)
            np.testing.assert_allclose(
                interpolated.velocity,midpoint_source.velocity,
                atol=1e-10,rtol=0)
            self.assertEqual(reader.seek(.5).time,.5)
            self.assertEqual(reader.seek(999).time,.5)
            self.assertEqual(reader.seek(-1).time,0.)
            with self.assertRaises(IndexError):
                reader.frame(6)
            with self.assertRaises(FileExistsError):
                precalculate_replay(config,path,fps=10.,chunk_frames=2)
            reader.close()

    def test_disc_impulses_reconstruct_collisions_between_frames(self):
        config=RunConfig(preset="carnot_discs",particles=16,seed=123,
                         duration=.5)
        with tempfile.TemporaryDirectory() as directory:
            path=precalculate_replay(config,Path(directory)/"replay",fps=20.,
                                     chunk_frames=2)
            reader=ReplayReader(path)
            direct=load_preset(config)
            try:
                for index in range(10):
                    midpoint=(index+.5)/20.
                    actual=direct.advance_to(midpoint)
                    visual,_=reader.sample(midpoint)
                    np.testing.assert_allclose(visual.position,actual.position,
                                               atol=1e-10,rtol=0)
                    np.testing.assert_allclose(visual.velocity,actual.velocity,
                                               atol=1e-10,rtol=0)
                self.assertGreater(reader.manifest["event_count"],0)
            finally:
                reader.close()

            # A previous version-1 archive lacks the new compact impulse and
            # mass columns, but retains the mandatory event log and config.
            manifest_path=path/"manifest.json"
            manifest=json.loads(manifest_path.read_text())
            manifest["geometry"].pop("mass")
            manifest["geometry"].pop("dynamic")
            manifest_path.write_text(json.dumps(manifest))
            index_path=path/"event_times.npz"
            with np.load(index_path,allow_pickle=False) as archive:
                old_columns={name:archive[name] for name in archive.files
                             if name!="impulse"}
            np.savez_compressed(index_path,**old_columns)
            legacy=ReplayReader(path)
            try:
                actual=load_preset(config).advance_to(.225)
                visual,_=legacy.sample(.225)
                np.testing.assert_allclose(visual.position,actual.position,
                                           atol=1e-10,rtol=0)
            finally:
                legacy.close()

    def test_prefetch_cache_stays_bounded_past_initial_chunks(self):
        config=RunConfig(preset="carnot_discs",particles=4,seed=123,
                         duration=.7)
        with tempfile.TemporaryDirectory() as directory:
            path=precalculate_replay(config,Path(directory)/"replay",fps=10.,
                                     chunk_frames=2)
            reader=ReplayReader(path)
            try:
                reader.prepare_playback()
                self.assertGreater(len(reader.chunks),reader.CACHE_CHUNKS)
                for index in range(len(reader)):
                    self.assertAlmostEqual(reader.frame(index).time,index/10)
                    self.assertLessEqual(len(reader._cache),reader.CACHE_CHUNKS)
            finally:
                reader.close()

    def test_triangle_event_points_reconstruct_angular_motion(self):
        config=RunConfig(preset="carnot_triangles",particles=16,seed=123,
                         duration=.5)
        with tempfile.TemporaryDirectory() as directory:
            path=precalculate_replay(config,Path(directory)/"replay",fps=20.,
                                     chunk_frames=2)
            reader=ReplayReader(path)
            direct=load_preset(config)
            try:
                self.assertTrue(reader._can_reconstruct_motion())
                for index in range(10):
                    midpoint=(index+.5)/20.
                    actual=direct.advance_to(midpoint)
                    visual,_=reader.sample(midpoint)
                    np.testing.assert_allclose(visual.position,actual.position,
                                               atol=1e-8,rtol=0)
                    np.testing.assert_allclose(visual.velocity,actual.velocity,
                                               atol=1e-8,rtol=0)
                    np.testing.assert_allclose(visual.omega,actual.omega,
                                               atol=1e-8,rtol=0)
                    angle_error=np.remainder(visual.angle-actual.angle+math.pi,
                                             2*math.pi)-math.pi
                    np.testing.assert_allclose(angle_error,0.,atol=1e-8,rtol=0)
                    direct.advance_to((index+1)/20.)
                self.assertGreater(reader.manifest["event_count"],0)
                with np.load(path/"event_times.npz",allow_pickle=False) as archive:
                    impacts=np.all(np.isfinite(archive["point"]),axis=1) & (
                        np.any(archive["impulse"] != 0,axis=1))
                    self.assertTrue(np.any(impacts))
            finally:
                reader.close()

    def test_precalculate_cli_writes_readable_archive(self):
        with tempfile.TemporaryDirectory() as directory:
            output=Path(directory)/"replay"
            self.assertEqual(main(["precalculate","--preset","carnot_discs",
                                   "--particles","2","--duration","0.12",
                                   "--fps","20","--chunk-frames","2",
                                   "--output",str(output)]),0)
            reader=ReplayReader(output)
            self.assertEqual(len(reader),4)
            self.assertAlmostEqual(reader.frame(3).time,.12)


if __name__=="__main__":
    unittest.main()
