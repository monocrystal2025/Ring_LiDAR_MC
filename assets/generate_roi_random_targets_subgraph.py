from __future__ import annotations

import math
import random
from pathlib import Path

from PIL import Image, ImageDraw


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "assets" / "random_targets_roi_subgraph.png"

W, H = 1200, 1200
SCALE = 3
CX = W / 2
CY = H / 2
R = 430


def sxy(p: tuple[float, float]) -> tuple[int, int]:
    return int(round(p[0] * SCALE)), int(round(p[1] * SCALE))


def draw_circle(draw: ImageDraw.ImageDraw, radius: float, fill, outline, width=2):
    box = [
        int((CX - radius) * SCALE),
        int((CY - radius) * SCALE),
        int((CX + radius) * SCALE),
        int((CY + radius) * SCALE),
    ]
    if fill is not None:
        draw.ellipse(box, fill=fill)
    draw.ellipse(box, outline=outline, width=width * SCALE)


def draw_polyline(draw: ImageDraw.ImageDraw, pts, fill, width=2):
    draw.line([sxy(p) for p in pts], fill=fill, width=width * SCALE, joint="curve")


def draw_arrow(
    draw: ImageDraw.ImageDraw,
    p0: tuple[float, float],
    direction: tuple[float, float],
    fill=(16, 112, 158, 235),
    length=68,
    width=5,
    head=16,
):
    dx, dy = direction
    norm = math.hypot(dx, dy)
    if norm < 1e-6:
        return
    ux, uy = dx / norm, dy / norm
    x0, y0 = p0
    x1, y1 = x0 + length * ux, y0 + length * uy
    draw.line([sxy((x0, y0)), sxy((x1, y1))], fill=fill, width=width * SCALE)
    ang = math.atan2(uy, ux)
    left = (
        x1 - head * math.cos(ang - math.pi / 6),
        y1 - head * math.sin(ang - math.pi / 6),
    )
    right = (
        x1 - head * math.cos(ang + math.pi / 6),
        y1 - head * math.sin(ang + math.pi / 6),
    )
    draw.polygon([sxy((x1, y1)), sxy(left), sxy(right)], fill=fill)


def roi_point(theta: float) -> tuple[float, float]:
    return math.cos(theta), math.sin(theta)


def inward_velocity(p: tuple[float, float]) -> tuple[float, float]:
    theta_v = 2 * math.pi * random.random()
    v = (math.cos(theta_v), math.sin(theta_v))
    if p[0] * v[0] + p[1] * v[1] > 0:
        v = (-v[0], -v[1])
    return v


def main():
    random.seed(53)
    img = Image.new("RGBA", (W * SCALE, H * SCALE), (255, 255, 255, 255))
    draw = ImageDraw.Draw(img, "RGBA")

    # ROI plane: a clean circular region matching a = z*tand(jiaodu).
    draw_circle(draw, R, (255, 255, 255, 255), (54, 91, 121, 245), 6)

    targets = []
    for _ in range(72):
        theta = 2 * math.pi * random.random()
        p = roi_point(theta)
        v = inward_velocity(p)
        targets.append((theta, p, v))

    targets.sort(key=lambda item: item[0])

    for _, p, v in targets:
        px = CX + R * p[0]
        py = CY + R * p[1]
        draw_arrow(draw, (px, py), (v[0], v[1]))

    for _, p, _ in targets:
        px = CX + R * p[0]
        py = CY + R * p[1]
        r = 9
        draw.ellipse(
            [
                int((px - r) * SCALE),
                int((py - r) * SCALE),
                int((px + r) * SCALE),
                int((py + r) * SCALE),
            ],
            fill=(230, 76, 62, 245),
            outline=(125, 32, 32, 245),
            width=1 * SCALE,
        )

    img = img.resize((W, H), Image.Resampling.LANCZOS)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    img.save(OUT)
    print(OUT)


if __name__ == "__main__":
    main()
