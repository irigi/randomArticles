# Replay archive v1

`microthermo precalculate` writes a directory with `manifest.json`, one
compressed NumPy archive per frame chunk, `events.jsonl`, and a compressed
event-time/participant index. The manifest is
written last. The command requires an empty output directory and keeps at
most `chunk_frames` snapshots in memory while solving. It now drains committed
events to JSONL after every saved frame, then builds the version-1 event index
through temporary disk-backed arrays. The completed temporary arrays are
removed after the compressed index is written. Ten-second 200/500-disc
regeneration passed exact event-log, saved-frame, ledger, and sampled-motion
parity; generation peaked at 150/155 MiB, respectively. The measurements
are in `docs/CARNOT_GAP_REPLAY_DISPLAY_2026-10-04.md`.

`manifest.json` records the format/version, complete run configuration and
hash, seed, requested backend, duration, sampling rate, frame and event
counts, chunk offsets and time bounds, particle shape/radius/polygon data,
initial wall geometry, and the final compensated energy ledger. The JSONL
file stores every event in the existing exact export representation.

Each `frames_XXXXXX.npz` is a compressed columnar chunk. Its columns are
time, position, velocity, angle, angular velocity, energy, translational and
rotational temperature, branch, occupancy, event count, energy residual,
shaft phase, apparatus components, and cumulative hot heat, cold heat, motor
work, and load output. Physical numbers are stored as float64. Every saved
frame is an independent keyframe: no quantization or delta prediction is
used. Reconstruction error at saved timestamps is zero relative to the
source snapshot. This choice spends more storage to keep the first format
auditable. Interpolated display states will have a separate error model.

`ReplayReader.frame(index)` loads only the selected chunk;
`ReplayReader.seek(time)` returns the latest saved frame at or before that
time. `ReplayReader.ledger(index)` returns the cumulative recorded ledger
values at that frame. `ReplayReader.sample(time)` linearly interpolates
apparatus coordinates between saved frames. It reconstructs particle motion
from recorded impulses. New version-1 archives include mass, inertia, dynamic
flags, and a compact contact-point column so angular motion can also be
reconstructed. Older disc archives recover masses and impulses from the
saved configuration and exact JSONL event log. Older non-disc archives lack
contact points and retain approximate particle interpolation. Branch state
changes at the recorded transition time. Ledger and temperature readouts use
the preceding exact frame with its timestamp.

For a controlled Carnot shaft, the interpolated apparatus points have a
conservative real-arithmetic error bound. A twice-differentiable coordinate
with second derivative at most `M` differs from its linear interpolation by
at most `M Δt²/8`. The cam uses a quintic step whose maximum absolute second
derivative is `10√3/3`. With the default equal sectors, shaft speed 0.15,
and 60 saved frames/s, the moving piston/cylinder coordinates are bounded
by `3.085e-6` length units for discs and `5.641e-6` for triangles. The
rotating apparatus curves have a smaller `2.11e-7` point-distance bound,
using their maximum 0.27 radius. The selector shoe jumps exactly at the
logged branch event; it is held at the appropriate sector detent across a
transition frame. These bounds apply to the default controlled cam and
exclude float64 roundoff. Other cam fractions, shaft speeds, and frame rates
require the same formula with their own sector width, area change, and
`Δt`. They do not cover free-shaft acceleration or older archives read with
approximate triangle motion.

Short warmed seed-123 disc probes at 60 saved frames per physical second:

| Particles | Physical window | Frames | Archive bytes | Offline write | Cold seek |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 200 | 0.10 s | 7 | 175,466 | 0.688 s | 2.93 ms |
| 500 | 0.05 s | 4 | 208,433 | 2.681 s | 3.01 ms |

The desktop `microthermo replay PATH` command provides play, pause, seek,
and 0.25–4× physical-time speed controls. An offscreen 0.5-second archive
probe drew 31 frames at 131 frames/s for 200 discs and 82 frames/s for 500
discs, including forced widget grabs; archive window load took 0.065 and
0.280 seconds. These short offscreen probes do not establish sustained
real-display frame rate, long-run storage, or typical seek latency. The
older solver/export path retained its event history until the run ended;
the new streaming path passed ten-second 200/500-disc memory checks. Exact
event-time trajectory reconstruction for new triangle archives now passes
the measured source comparisons in the replay display report.

The later ten-second `xcb` desktop measurement averaged 60.7–61.1 painted
frames/s at 200 and 500 discs, but had three paint gaps over 25 ms in each
archive, reaching 149–183 ms. Cross-chunk seeks blocked for up to 55 ms.
Measurements and remaining work are in
`docs/CARNOT_GAP_REPLAY_DISPLAY_2026-10-04.md`.

The reader now warms the event index and up to three frame chunks before the
desktop window appears. Later chunks are prefetched in a background worker,
with a three-chunk cache. This moves the initial decompression cost into
window load and aims to avoid playback pauses; it does not guarantee instant
random seeks to chunks outside the cache. The repeated ten-second `xcb` check
had no paint gap above 18 ms at 200 or 500 discs; all 145 unit tests passed.
Longer random seeks and between-frame error remain open. Details are in
`docs/CARNOT_GAP_REPLAY_DISPLAY_2026-10-04.md`.

Event-aware disc reconstruction now gives roundoff-level midpoint position
error against independent 200/500-disc source runs over one physical second,
and exact midpoint velocity in those checks. A ten-second real-display replay
with the change averaged about 62.5 paints/s but showed one roughly 40 ms
gap per count. All 146 tests passed. Consistent smoothness, uncached seeks,
and triangle angular replay remain open.

A subsequent eight-job batch linked the two long 1× paint gaps to Python
cyclic garbage collections of 24–29 ms. The replay window now suspends cyclic
collection only while playback is active and restores its prior state on
pause or close. The batch also found triangle midpoint errors up to one
particle radius and π radians in archives without contact points. New
archives now record contact points and reconstruct translation and rotation
through impulses. Focused tests pass; independent high-count and desktop
checks after these changes are pending.
