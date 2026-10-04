from __future__ import annotations

import math
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont


SCALE = 2
WIDTH, HEIGHT = 2400, 1800
LEAD_PULSES = 2
BASE_N = 33
N = BASE_N + LEAD_PULSES
WINDOW_LEN = 9
THRESHOLD = 2.0

NAVY = (7, 31, 79)
MAGENTA = (196, 0, 70)
MAGENTA_DARK = (139, 0, 48)
GRAY = (185, 190, 199)
GRAY_DARK = (74, 79, 87)
BLUE = (42, 132, 211)
BLUE_DARK = (11, 63, 124)
TERM_BLUE = (78, 161, 217)
TERM_BLUE_DARK = (11, 95, 153)
GUIDE = (183, 195, 216)
CYAN_FILL = (191, 239, 255, 70)
RED = (208, 0, 69)


def font(size: int, bold: bool = False) -> ImageFont.FreeTypeFont:
    candidates = ["Helvetica-Bold.ttf", "helveticab.ttf", "arialbd.ttf"] if bold else []
    candidates.extend(["Helvetica.ttf", "helvetica.ttf", "arial.ttf", "DejaVuSans.ttf"])
    for name in candidates:
        try:
            return ImageFont.truetype(name, size * SCALE)
        except OSError:
            pass
    return ImageFont.load_default()


# Rendered at 2x resolution; these map to about 20 pt in a 300 dpi export.
F_AXIS = font(40)
F_Y = font(42, bold=True)
F_LEGEND = font(40)


def draw_text(draw: ImageDraw.ImageDraw, xy: tuple[float, float], text: str, fnt, fill=NAVY):
    draw.text((xy[0] * SCALE, xy[1] * SCALE), text, font=fnt, fill=fill)


def draw_rotated_text(img: Image.Image, xy: tuple[float, float], text: str, fnt, fill=NAVY):
    tmp = Image.new("RGBA", (520 * SCALE, 70 * SCALE), (255, 255, 255, 0))
    tdraw = ImageDraw.Draw(tmp)
    tdraw.text((0, 0), text, font=fnt, fill=fill)
    tmp = tmp.rotate(90, expand=True)
    img.alpha_composite(tmp, (int(xy[0] * SCALE), int(xy[1] * SCALE)))


def line(draw: ImageDraw.ImageDraw, xy, fill=NAVY, width=2, dash: tuple[int, int] | None = None):
    xy = tuple(int(round(v * SCALE)) for v in xy)
    if dash is None:
        draw.line(xy, fill=fill, width=width * SCALE)
        return
    x1, y1, x2, y2 = xy
    dx, dy = x2 - x1, y2 - y1
    length = math.hypot(dx, dy)
    if length <= 0:
        return
    ux, uy = dx / length, dy / length
    pos = 0.0
    on, off = dash[0] * SCALE, dash[1] * SCALE
    while pos < length:
        end = min(pos + on, length)
        draw.line(
            (
                int(x1 + ux * pos),
                int(y1 + uy * pos),
                int(x1 + ux * end),
                int(y1 + uy * end),
            ),
            fill=fill,
            width=width * SCALE,
        )
        pos += on + off


def rect(draw: ImageDraw.ImageDraw, box, fill=None, outline=None, width=1, radius=0):
    box = tuple(int(round(v * SCALE)) for v in box)
    if radius:
        draw.rounded_rectangle(box, radius=radius * SCALE, fill=fill, outline=outline, width=width * SCALE)
    else:
        draw.rectangle(box, fill=fill, outline=outline, width=width * SCALE)


def axis(draw: ImageDraw.ImageDraw, left: float, base: float, right: float, top: float):
    line(draw, (left, base, right, base), NAVY, 2)
    line(draw, (left, base, left, top), NAVY, 2)
    draw.polygon(
        [
            (right * SCALE, base * SCALE),
            ((right - 10) * SCALE, (base - 6) * SCALE),
            ((right - 10) * SCALE, (base + 6) * SCALE),
        ],
        fill=NAVY,
    )
    draw.polygon(
        [
            (left * SCALE, top * SCALE),
            ((left - 6) * SCALE, (top + 10) * SCALE),
            ((left + 6) * SCALE, (top + 10) * SCALE),
        ],
        fill=NAVY,
    )


def ring_signal() -> list[float]:
    values = []
    for i in range(1, N + 1):
        j = i - LEAD_PULSES
        if j <= 0:
            values.append(0.0)
            continue
        lobe1 = 1.08 * math.exp(-0.5 * ((j - 11) / 2.25) ** 2)
        lobe2 = 1.12 * math.exp(-0.5 * ((j - 23) / 2.35) ** 2)
        values.append(lobe1 + lobe2)
    return values


def render(out_paths: list[Path]):
    img = Image.new("RGBA", (WIDTH * SCALE, HEIGHT * SCALE), "white")
    overlay = Image.new("RGBA", img.size, (255, 255, 255, 0))
    draw = ImageDraw.Draw(img)
    odraw = ImageDraw.Draw(overlay)

    signal = ring_signal()
    noise = [0.17 + 0.012 * math.sin((i + 1) * 0.72) for i in range(N)]
    weights = [s / (s + n) if s + n > 0 else 0.0 for s, n in zip(signal, noise)]
    numerator = [w * s for w, s in zip(weights, signal)]
    variance = [w * w * (s + n) for w, s, n in zip(weights, signal, noise)]

    statistic: list[float] = []
    first_cross = N - 1
    found = False
    for k in range(N):
        start = max(0, k - WINDOW_LEN + 1)
        num_sum = sum(numerator[start : k + 1])
        var_sum = sum(variance[start : k + 1])
        value = num_sum / math.sqrt(var_sum) if var_sum > 0 else 0.0
        statistic.append(value)
        if not found and value >= THRESHOLD:
            first_cross = k
            found = True

    win_start = max(0, first_cross - WINDOW_LEN + 1)

    left, right = 245.0, 2220.0
    rows = [125.0, 515.0, 905.0, 1295.0]
    plot_h = 245.0
    bar_w = 18.0
    plot_w = right - left

    def xpos(i: int) -> float:
        return left + i * plot_w / (N - 1)

    for row in rows:
        axis(draw, left, row + plot_h, right, row + 10)

    for i in range(0, N, 3):
        x = xpos(i)
        line(draw, (x, rows[0] + 10, x, rows[-1] + plot_h), GUIDE, 1, dash=(5, 9))
    x_cross = xpos(first_cross)
    line(draw, (x_cross, rows[0] + 10, x_cross, rows[-1] + plot_h), RED, 2, dash=(8, 8))

    wx0 = xpos(win_start) - bar_w
    wx1 = xpos(first_cross) + bar_w
    for row in rows[:3]:
        rect(odraw, (wx0, row + 15, wx1, row + plot_h), fill=CYAN_FILL, outline=(0, 164, 223, 180), width=2)

    # Row 1
    base = rows[0] + plot_h
    scale1 = 155.0 / max(s + n for s, n in zip(signal, noise))
    for i, (s, n) in enumerate(zip(signal, noise)):
        x = xpos(i) - bar_w / 2
        nh, sh = n * scale1, s * scale1
        rect(draw, (x, base - nh, x + bar_w, base), fill=GRAY, outline=GRAY_DARK, width=1)
        rect(draw, (x, base - nh - sh, x + bar_w, base - nh), fill=MAGENTA, outline=MAGENTA_DARK, width=1)

    # Row 2
    base = rows[1] + plot_h
    for i, w in enumerate(weights):
        x = xpos(i) - bar_w / 2
        h = w * 165.0
        rect(draw, (x, base - h, x + bar_w, base), fill=BLUE, outline=BLUE_DARK, width=1)

    # Row 3, overlaid equal-width terms
    base = rows[2] + plot_h
    scale3 = 160.0 / max(numerator)
    for i, value in enumerate(numerator):
        x = xpos(i) - bar_w / 2
        h = value * scale3
        alpha = 255 if win_start <= i <= first_cross else 65
        box = (x, base - h, x + bar_w, base)
        rect(odraw, box, fill=(*TERM_BLUE, alpha), outline=(*TERM_BLUE_DARK, alpha), width=1)
        rect(odraw, box, fill=(*MAGENTA, int(alpha * 0.45)), outline=(*MAGENTA_DARK, alpha), width=1)

    # Row 4
    base = rows[3] + plot_h
    scale4 = 165.0 / max(max(statistic), THRESHOLD * 1.15)
    thresh_y = base - THRESHOLD * scale4
    line(draw, (left, thresh_y, right, thresh_y), RED, 2, dash=(12, 9))
    points = [(xpos(i), base - value * scale4) for i, value in enumerate(statistic)]
    draw.line([(x * SCALE, y * SCALE) for x, y in points], fill=MAGENTA, width=5 * SCALE)
    for i in range(0, N, 2):
        x, y = points[i]
        draw.ellipse(
            ((x - 5) * SCALE, (y - 5) * SCALE, (x + 5) * SCALE, (y + 5) * SCALE),
            fill=MAGENTA,
        )
    x, y = points[first_cross]
    draw.ellipse(
        ((x - 9) * SCALE, (y - 9) * SCALE, (x + 9) * SCALE, (y + 9) * SCALE),
        fill="white",
        outline=MAGENTA,
        width=4 * SCALE,
    )

    ylabels = ["Photons", "Weight", ("Pulse", "contribution"), "Statistic"]
    xlabels = ["Pulse index", "Pulse index", "Pulse index", "Window index"]
    ylabel_x = left + 18
    for row, yl, xl in zip(rows, ylabels, xlabels):
        if isinstance(yl, tuple):
            draw_text(draw, (ylabel_x, row + 14), yl[0], F_Y)
            draw_text(draw, (ylabel_x, row + 58), yl[1], F_Y)
        else:
            draw_text(draw, (ylabel_x, row + 18), yl, F_Y)
        draw_text(draw, (right - 180, row + plot_h + 16), xl, F_AXIS)

    # Legends only; no formulas or explanatory paragraphs.
    legend_x = 1760
    rect(draw, (legend_x, rows[0] - 12, legend_x + 25, rows[0] + 13), fill=MAGENTA, outline=MAGENTA_DARK, width=1)
    draw_text(draw, (legend_x + 38, rows[0] - 18), "Signal", F_LEGEND)
    rect(draw, (legend_x + 180, rows[0] - 12, legend_x + 205, rows[0] + 13), fill=GRAY, outline=GRAY_DARK, width=1)
    draw_text(draw, (legend_x + 218, rows[0] - 18), "Noise", F_LEGEND)

    line(draw, (legend_x, rows[3], legend_x + 90, rows[3]), RED, 2, dash=(10, 8))
    draw_text(draw, (legend_x + 112, rows[3] - 22), "Threshold", F_LEGEND)

    img = Image.alpha_composite(img, overlay)
    for out_path in out_paths:
        out_path.parent.mkdir(parents=True, exist_ok=True)
        img.convert("RGB").save(out_path)
        print(out_path)


if __name__ == "__main__":
    asset_dir = Path(__file__).resolve().parent
    render(
        [
            asset_dir / "gpt_ring_pulse_window_snr_decision.png",
            asset_dir / "aligned_ring_pulse_window_snr_decision.png",
            asset_dir / "ring_pulse_window_snr_decision_compact.png",
        ]
    )
