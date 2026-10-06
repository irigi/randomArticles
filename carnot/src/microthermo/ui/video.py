"""Render a precalculated replay to MP4: the replay window drawn offscreen, piped to ffmpeg."""

from __future__ import annotations

import math
from pathlib import Path
import shutil
import subprocess

from PySide6 import QtGui, QtWidgets

from .replay_window import ReplayWindow


def export_video(path, output, *, fps: float = 30., speed: float = 1., start: float = 0.,
                 end: float | None = None, size: tuple[int, int] = (1480, 900),
                 crf: int = 20, progress=None) -> Path:
    """Write the replay window, plots included, from ``start`` to ``end``.

    ``speed`` is physical seconds per video second. Every video frame redraws
    the plots, so the video shows what a paused replay shows at that time.
    ``progress(frame, frames)`` is called after each written frame.
    """
    width, height = size
    if width <= 0 or height <= 0 or width % 2 or height % 2:
        raise ValueError("video width and height must be positive and even")
    if not (math.isfinite(fps) and fps > 0 and math.isfinite(speed) and speed > 0):
        raise ValueError("fps and speed must be positive")
    ffmpeg = shutil.which("ffmpeg")
    if ffmpeg is None:
        raise RuntimeError("ffmpeg is not on PATH")
    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    window = ReplayWindow(path)
    try:
        duration = window.duration
        end = duration if end is None else min(end, duration)
        if not 0 <= start < end:
            raise ValueError("need 0 <= start < end within the archive")
        window.controls.setVisible(False)
        window.resize(width, height)
        window.show()
        app.processEvents()
        frames = int(math.floor((end - start)*fps/speed)) + 1
        output = Path(output)
        output.parent.mkdir(parents=True, exist_ok=True)
        command = [ffmpeg, "-y", "-loglevel", "error", "-f", "rawvideo", "-pix_fmt", "rgb24",
                   "-s", f"{width}x{height}", "-r", f"{fps:g}", "-i", "-",
                   "-c:v", "libx264", "-pix_fmt", "yuv420p", "-crf", str(crf),
                   "-movflags", "+faststart", str(output)]
        encoder = subprocess.Popen(command, stdin=subprocess.PIPE)
        try:
            for k in range(frames):
                window.physical_time = min(end, start + k*speed/fps)
                window._render()
                app.processEvents()
                image = window.grab().toImage().convertToFormat(
                    QtGui.QImage.Format.Format_RGB888)
                if (image.width(), image.height()) != (width, height):
                    image = image.scaled(width, height)
                line = image.bytesPerLine()
                data = bytes(image.constBits())
                if line != 3*width:
                    data = b"".join(data[y*line:y*line + 3*width] for y in range(height))
                encoder.stdin.write(data)
                if progress is not None:
                    progress(k + 1, frames)
        finally:
            encoder.stdin.close()
            code = encoder.wait()
        if code:
            raise RuntimeError(f"ffmpeg exited with status {code}")
        return output
    finally:
        window.close()
