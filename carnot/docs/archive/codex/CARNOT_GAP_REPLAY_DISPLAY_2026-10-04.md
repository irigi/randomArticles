# Real-display replay check — 2026-10-04

The existing ten-physical-second, 60-saved-frame/s disc archives were played
sequentially in an `xcb` Qt desktop window. The benchmark counted apparatus
paint events while playing each archive at 1×, and timed six seeks across its
three compressed frame chunks. It ran on the previously declared i7-11850H
machine. The small raw JSON summaries are saved as
`runs/gap_followup_2026-10-04/replay_display_{200,500}.json`.

| Discs | Paints | Paint rate | Median interval | 95th percentile | Longest gap | Gaps >25 ms | Window load | Median / max seek | Peak RSS |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 200 | 611 | 61.1/s | 16.0 ms | 16.7 ms | 149.4 ms | 3 | 0.216 s | 39.9 / 42.9 ms | 282 MiB |
| 500 | 607 | 60.7/s | 16.1 ms | 16.7 ms | 183.1 ms | 3 | 0.434 s | 50.5 / 54.9 ms | 297 MiB |

Both commands exited successfully and reached the archive's physical end
time. The median and 95th-percentile cadence is near the 60-frame/s target,
but the three visible-length pauses per archive mean the smooth-playback gate
is not yet closed. The benchmark did not record the timestamps of individual
pauses. Synchronous compressed-chunk changes and the reader's lazy event-index
load are plausible causes from code inspection, not established causes from
this timing data. Seeking also blocks the display thread for up to about
55 ms. The saved-frame fidelity tests passed in the prior unit suite; exact
between-frame trajectory error has not been bounded.

Next replay work: instrument pause timestamps against chunk boundaries and
first event-index load; move decompression and index preparation off the
display path or cache them ahead of playback; then repeat the same real-display
check and quantify interpolation error near collisions. Archive generation
still retains the full event history in memory until it finishes.

## After cache and prefetch change

The reader now loads the event index and the first three frame chunks before
showing the window, retains at most three decoded chunks, and prefetches later
chunks in a worker. The same saved archives and six seek positions were tested
again on the `xcb` desktop. Compact results are in
`runs/gap_followup_2026-10-04/replay_display_warm_{200,500}.json`.

| Discs | Paints | Paint rate | 95th-percentile interval | Longest gap | Gaps >25 ms | Window load | Median / max seek | Peak RSS |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 200 | 625 | 62.5/s | 16.70 ms | 17.54 ms | 0 | 0.282 s | 5.62 / 19.80 ms | 288 MiB |
| 500 | 625 | 62.5/s | 16.67 ms | 17.35 ms | 0 | 0.522 s | 6.97 / 25.14 ms | 304 MiB |

Both replay jobs exited successfully and reached physical time 10.0 seconds.
The full unittest suite passed 145 tests. For these ten-second archives at 1×,
the observed smoothness regression is resolved. Window load increased by
0.066 seconds at 200 and 0.089 seconds at 500; peak RSS increased by about
6.5 and 7.4 MiB. Random seeks outside the three-chunk cache in a longer
archive may still block, and the current interpolated particle coordinates
have no quantified error bound near collisions. Those remain part of the
replay acceptance work, along with removing full event-history retention
during precalculation.

## Midpoint trajectory check

New one-physical-second 60-frame/s archives were independently compared with
source simulations at all 60 half-frame times. Both jobs exited successfully;
all saved-frame particle positions and angles matched the independently
sampled simulations exactly. Every frame interval contained an event.

| Discs | Midpoint particle observations | 95th-percentile position error | Maximum position error | Maximum error / particle radius | Maximum apparatus point error |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 200 | 12,000 | 0.00645 | 0.02115 | 1.73 | 5.8e-7 |
| 500 | 30,000 | 0.00803 | 0.02623 | 3.39 | 5.8e-7 |

The current straight-line interpolation across particle collisions is too
inaccurate for faithful visual trajectories. The exact event archive already
contains impulse vectors and participant IDs, so disc positions can be
reconstructed piecewise between events. This is the next implementation
step; it will need a repeat of the midpoint and desktop-rate measurements.
Compact metrics are saved in
`runs/gap_followup_2026-10-04/replay_error_{200,500}.json`.

Event-aware disc reconstruction is now implemented, including a fallback for
the existing archives. Five focused replay tests pass. Independent high-count
midpoint and desktop-rate reassessment after this change are pending.

## Event-aware disc verification

The original ten-second archives were compared with independent source
simulations through their first physical second. Both source jobs exited
successfully, covered 60 midpoint times and event-containing intervals, and
matched all saved-frame positions exactly. Maximum midpoint position error
fell from 0.02115 to `9.49e-16` at 200 discs and from 0.02623 to `2.02e-15`
at 500; midpoint velocity error was zero at both counts. This is roundoff
level for these cases, though it is an observed maximum over one second and
one seed, not a universal analytic bound.

The repeated `xcb` desktop check still reached about 62.5 paints/s for both
ten-second archives, with 95th-percentile paint intervals of 16.60 and
16.77 ms. Each run had **one** interval above 25 ms, reaching 39.25 ms at
200 discs and 40.68 ms at 500. These occurred at physical times 7.95 and
3.52 seconds, away from the archive's 4.27 and 8.53-second chunk transitions.
That timing does not identify a cause; a scheduler pause or frame-specific
render cost remains possible. Window load rose to 0.444 and 1.073 seconds,
consistent with parsing old archives' JSONL impulse logs during preparation.
Maximum measured seek was 19.2 and 11.6 ms. The full 146-test suite passed.

The compact source and display metrics are saved in
`runs/gap_followup_2026-10-04/replay_exact_{source,display}_{200,500}.json`.
The ten-second disc replay is now accurate at measured midpoints and meets
the average paint-rate target, but the one-off long paint gaps need another
timed diagnosis before claiming consistently smooth playback. Longer archives
with seeks outside the three-chunk cache, other seeds, triangle angular
replay, and release evidence are still open.

## Batched replay diagnosis

Eight further jobs completed successfully: 32/96-triangle midpoint checks,
200/500-disc random seeks in 19-chunk repacked copies of the same ten-second
archives, and 200/500-disc desktop playback at 1× and 4×. The repacking was
only for cache stress and did not run the solver again.

| Job | Result |
| --- | --- |
| 200-disc uncached seek | 8.25 ms maximum; immediate repeat 1.09 ms maximum |
| 500-disc uncached seek | 9.61 ms maximum; immediate repeat 1.10 ms maximum |
| 200-disc display at 1× / 4× | 62.38 / 62.49 paints/s |
| 500-disc display at 1× / 4× | 62.39 / 62.47 paints/s |

At 1×, the longest paint gaps were 44.62 ms for 200 discs and 41.25 ms for
500. The corresponding render calls took 30.30 and 26.27 ms, coinciding with
measured garbage-collection pauses of 28.79 and 24.50 ms. This strongly
implicates Python cyclic collection as the cause of these long replay stalls.
The 4× runs had no garbage-collection pause over 5 ms; their longest paint
intervals were 21.41 and 26.26 ms. A 25.72 ms interval near the start of the
500-disc 1× run also lacked a slow render call, so occasional compositor or
scheduling jitter remains possible.

The triangle checks confirmed exact saved-frame positions but poor current
between-frame replay. Maximum midpoint position error was 0.523 particle
radii at 32 triangles and 1.024 radii at 96; maximum angular error was π
radians at both counts and midpoint velocity error reached 2.13 and 4.44.
The archives lack contact points needed to reconstruct angular impulses.
The next implementation should record those points for new archives and use
them for event-aware triangle motion. The raw metrics are in
`runs/gap_followup_2026-10-04/replay_next_{triangle,seek,display}_*.json`.

Playback now temporarily suspends cyclic garbage collection and restores its
prior state when paused or closed. New triangle archives record contact points
and use them with impulses to reconstruct angular motion. Focused tests pass;
the next batch must confirm the 32/96-triangle error reduction and repeat the
desktop timing before these changes count as acceptance evidence.

## Acceptance batch after garbage-collection and triangle changes

All ten jobs exited successfully. Eight compact JSON result files are in
`runs/gap_followup_2026-10-04/replay_acceptance_*.json`; the completed
temporary launcher and its raw logs were deleted. The two fresh triangle
archives used seed 123, 60 saved frames/s, and two physical seconds.

| Triangle count | Midpoints | Maximum position error | Maximum angle error | Maximum velocity error | Maximum angular-velocity error |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 32 | 120 | 4.58e-16 | 1.60e-14 rad | 0 | 3.13e-12 |
| 96 | 120 | 6.75e-16 | 6.08e-14 rad | 0 | 1.14e-11 |

Every saved-frame position matched the independently advanced source run.
Events occurred in 111 of 120 frame intervals at 32 triangles and all 120 at
96. Both readers reported event-aware reconstruction. These are observed
roundoff-level errors for one seed and duration, rather than a general bound.

| Archive | Speed | Paints/s | 95th-percentile interval | Longest gap | Gaps >25 ms | Window load | Peak RSS |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 200 discs, 10 s | 1× | 62.50 | 16.82 ms | 17.68 ms | 0 | 0.416 s | 288 MiB |
| 500 discs, 10 s | 1× | 62.50 | 16.82 ms | 23.14 ms | 0 | 1.198 s | 307 MiB |
| 200 discs, 10 s | 4× | 62.82 | 16.78 ms | 17.20 ms | 0 | 0.432 s | 288 MiB |
| 500 discs, 10 s | 4× | 62.52 | 16.88 ms | 22.93 ms | 0 | 1.196 s | 306 MiB |
| 32 triangles, 2 s | 1× | 62.50 | 16.56 ms | 17.21 ms | 0 | 0.178 s | 244 MiB |
| 96 triangles, 2 s | 1× | 62.51 | 16.58 ms | 17.46 ms | 0 | 0.172 s | 245 MiB |

No playback job recorded a cyclic-collection pause above 5 ms, and all
reported that cyclic collection was enabled after playback. The disc results
resolve the previously reproduced collection-related paint stalls in these
ten-second runs. The 147-test full unit suite and every scientific validation
check also passed. Remaining replay acceptance concerns include other seeds,
longer triangle archives, saved statistics after arbitrary seeks, and the
solver's full event-history retention during precalculation.

## Second-seed source and desktop verification

Seed-124 triangle archives were independently checked over five physical
seconds at 60 saved frames/s. Both 32- and 96-particle runs matched the
direct simulation at 300 midpoints. Maximum position errors were `4.97e-16`
and `8.90e-16`; maximum angle errors were `1.11e-14` and `1.00e-13` radians;
maximum angular-velocity errors were `2.05e-12` and `2.55e-11`. Midpoint
velocity error was zero. Saved-frame position, scalar statistics, event count,
branch, and cumulative ledger matched the source exactly in both runs.

The two `xcb` desktop runs reached 62.50 paints/s, with no interval above
25 ms and maxima of 17.28 and 17.22 ms. Six UI seeks per run displayed the
corresponding recorded event and heat values; the largest UI seek took 24.93
ms at 32 and 30.54 ms at 96 triangles. Garbage collection was restored after
playback. On the existing 200/500-disc archives, 60 shuffled frame/ledger
seeks per count matched a reference reader; their maximum warmed seek times
were 0.423 and 0.408 ms. These tests cover two triangle seeds and up to five
seconds, while archive-generation event memory and longer high-count playback
remain open. Compact JSON results are saved as
`runs/gap_followup_2026-10-04/replay_second_seed_*.json`.

## Streaming event export and 19-chunk acceptance

The new writer streams exact JSONL events after each saved frame and builds
the version-1 compact event index through temporary disk-backed columns.
Both ten-physical-second, seed-123 disc archives regenerated successfully at
200 and 500 particles, now with 19 frame chunks instead of three. The event
logs have the same SHA-256 hashes as the earlier archives. All 601 saved
frames, event counts, scalar statistics, apparatus states, and ledgers match
exactly; six sampled positions per count also match. Archive sizes were
14,054,318 and 32,720,838 bytes. Generation peak RSS was 153,788 and
158,248 KiB, down from approximately 214 and 243 MiB in the previous
ten-second generation measurements; wall times were 66.69 and 492.30 s.
The older archives used 256-frame chunks, so this RSS difference combines
event streaming with a smaller frame buffer; the data do not apportion the
reduction between those two changes.

| Archive | Speed | Paints/s | Longest interval | Gaps >25 ms | Max UI seek | Load time |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 200 discs | 1× | 62.50 | 16.85 ms | 0 | 30.49 ms | 0.233 s |
| 200 discs | 4× | 62.51 | 17.04 ms | 0 | 23.67 ms | 0.230 s |
| 500 discs | 1× | 62.51 | 17.64 ms | 0 | 17.20 ms | 0.527 s |
| 500 discs | 4× | 62.51 | 17.68 ms | 0 | 26.69 ms | 0.522 s |

All four desktop runs used `xcb`, reached physical time 10.0, and restored
garbage collection afterward. The full 148-test suite and every scientific
validation check passed. Compact JSON summaries are in
`runs/gap_followup_2026-10-04/replay_streaming_*.json`; the regenerated
archives and temporary logs were deleted after review. The measured replay
cases now satisfy the frame-rate and exact particle/statistics checks. At
that stage the apparatus interpolation bound was still pending; it is now
documented in `docs/REPLAY_FORMAT.md`.

## High-count triangle replay through a branch transition

New 200/500-triangle archives ran for 10.6 physical seconds, crossing the
first hot-to-expansion transition at about 10.472 seconds. Both have 637
saved frames in 20 chunks. The 200-triangle archive recorded 22,709 events
and occupied 19,084,285 bytes; the 500-triangle archive recorded 78,660
events and occupied 45,979,152 bytes. Their generation jobs exited with
zero status, peaking at 168,076 and 182,752 KiB RSS. The 500-triangle job
took 15:44 wall time. The first 60 midpoints matched independent source
motion to at most `2.03e-15` in position, `9.99e-13` radians in angle, and
`2.36e-10` in angular velocity; velocity error was zero. Saved frame
positions, scalar statistics, and ledgers matched exactly in that interval.

At 209 times spanning the archive and clustered around the branch change,
both counts matched the source branch and selector-shoe detent exactly.
Maximum apparatus point error was `7.92e-7`, below the conservative
`5.641e-6` controlled-cam bound in `docs/REPLAY_FORMAT.md`.

All four `xcb` desktop runs reached physical time 10.6, reported roughly
62.5 paints/s between their first and last counted paints, and had no
measured inter-paint gap above 25 ms. The longest measured gaps were
17.21/17.03 ms at 200 triangles (1×/4×) and 17.72/17.44 ms at 500
(1×/4×). Six UI seeks per run showed the expected branch and apparatus
state; maximum seek time was 18.89 ms. The full 149-test suite and
scientific validation passed. Compact results are preserved in
`runs/gap_followup_2026-10-04/replay_triangle_highcount_*.json`.

The 4× paint counts need a separate start-to-end timing check. In particular,
the 500-triangle 4× run counted 111 paints, covering about 1.76 seconds at
its measured cadence, while 10.6 physical seconds at 4× should take about
2.65 seconds. The current metric starts at the first paint and could miss a
long startup delay. The archives are temporarily retained for a display-only
diagnostic; no new simulation is needed. Until its first-paint and end-time
measurements are reviewed, the 4× smoothness claim and overall replay gate
remain open.

## Start-to-end timing diagnosis

Five display-only runs reused the same high-count triangle archives. Two
independent 4× repeats at each count started at physical time zero, entered
the Qt event loop within 0.25 ms, painted first at about 19 ms, and counted
166 ticks and 166 paints each. All four reached physical time 10.6 in
2.657–2.658 wall seconds against 2.65 seconds nominal. The longest tick
and paint gaps were below 17.5 ms. A 500-triangle 1× repeat counted 663
ticks/paints, reached the end in 10.610 wall seconds, and had no gap above
17.5 ms. Garbage collection was enabled again after all five runs.

The earlier 500-triangle 4× count of 111 paints did not recur. Its cause
cannot be identified from the original first-to-last-paint metric, which
did not record start or end timestamps; a one-off startup delay or measurement
artifact remains possible. The repeated start-to-end evidence supports the
display target for the tested 200/500-triangle archives. The measured replay
gate is accepted for these controlled-shaft presets, archive durations,
hardware, and Qt desktop; it is not a guarantee against every operating-
system scheduling pause. Compact timing JSON is in
`runs/gap_followup_2026-10-04/replay_triangle_timing_*.json`. Both retained
triangle archives, the raw logs, and all temporary launchers were deleted.
