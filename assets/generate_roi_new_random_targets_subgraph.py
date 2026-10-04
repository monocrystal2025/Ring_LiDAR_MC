from __future__ import annotations

import math
import random
from pathlib import Path

from PIL import Image, ImageDraw


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "assets" / "random_targets_roi_new_subgraph.png"

W, H = 1600, 1000
SCALE = 3

CX = 800
BASE_Y = 900
SX = 720
SY = 175
SZ = 760

JIAODU_DEG = 30
ALPHA = math.radians(JIAODU_DEG)
S_MIN = 0.08
S_MAX = 1.0


def sxy(p: tuple[float, float]) -> tuple[int, int]:
    return int(round(p[0] * SCALE)), int(round(p[1] * SCALE))


def project(x: float, y: float, z: float) -> tuple[float, float]:
    return CX + SX * x, BASE_Y + SY * y - SZ * z


def draw_polyline(draw: ImageDraw.ImageDraw, pts, fill, width=2):
    draw.line([sxy(p) for p in pts], fill=fill, width=width * SCALE, joint="curve")


def split_visible_segments(samples, is_front: bool):
    segments = []
    current = []
    for point, y_depth in samples:
        front = y_depth >= 0
        if front == is_front:
            current.append(point)
        elif current:
            segments.append(current)
            current = []
    if current:
        segments.append(current)
    first_matches = (samples[0][1] >= 0) == is_front
    last_matches = (samples[-1][1] >= 0) == is_front
    if len(segments) > 1 and first_matches and last_matches:
        segments[0] = segments[-1] + segments[0]
        segments.pop()
    return [segment for segment in segments if len(segment) > 1]


def cone_point(slant_range: float, theta: float) -> tuple[float, float, float]:
    return (
        slant_range * math.sin(ALPHA) * math.cos(theta),
        slant_range * math.sin(ALPHA) * math.sin(theta),
        slant_range * math.cos(ALPHA),
    )


def random_inward_velocity(theta: float) -> tuple[float, float, float]:
    theta_v = 2 * math.pi * random.random()
    u_v = 2 * random.random() - 1
    phi_v = math.acos(u_v)
    v = [
        math.sin(phi_v) * math.cos(theta_v),
        math.sin(phi_v) * math.sin(theta_v),
        math.cos(phi_v),
    ]
    normal = [math.cos(theta), math.sin(theta), -math.tan(ALPHA)]
    n_norm = math.sqrt(sum(value * value for value in normal))
    normal = [value / n_norm for value in normal]
    dot = sum(v[i] * normal[i] for i in range(3))
    if dot >= 0:
        v = [v[i] - 2 * dot * normal[i] for i in range(3)]
    return v[0], v[1], v[2]


def draw_arrow(
    draw: ImageDraw.ImageDraw,
    p0: tuple[float, float],
    direction: tuple[float, float],
    fill=(16, 112, 158, 235),
    length=54,
    width=4,
    head=13,
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


def main():
    random.seed(86)
    img = Image.new("RGBA", (W * SCALE, H * SCALE), (255, 255, 255, 255))
    draw = ImageDraw.Draw(img, "RGBA")

    # Back-side slant lines.
    for theta in [7 * math.pi / 6, 4 * math.pi / 3, 3 * math.pi / 2, 5 * math.pi / 3, 11 * math.pi / 6]:
        pts = [project(*cone_point(S_MIN + (S_MAX - S_MIN) * i / 140, theta)) for i in range(141)]
        draw_polyline(draw, pts, (64, 115, 153, 70), 2)

    # Rings along slant range; rear half first, front half later.
    ring_slants = [S_MIN, 0.28, 0.48, 0.68, 0.84, S_MAX]
    for slant in ring_slants:
        samples = []
        for i in range(241):
            theta = 2 * math.pi * i / 240
            point = cone_point(slant, theta)
            samples.append((project(*point), point[1]))
        width = 4 if slant in (S_MIN, S_MAX) else 2
        for pts in split_visible_segments(samples, is_front=False):
            draw_polyline(draw, pts, (64, 115, 153, 65), width)

    # Front-side slant lines.
    for theta in [0, math.pi / 6, math.pi / 3, math.pi / 2, 2 * math.pi / 3, 5 * math.pi / 6, math.pi]:
        pts = [project(*cone_point(S_MIN + (S_MAX - S_MIN) * i / 140, theta)) for i in range(141)]
        draw_polyline(draw, pts, (64, 115, 153, 125), 2)

    for slant in ring_slants:
        samples = []
        for i in range(241):
            theta = 2 * math.pi * i / 240
            point = cone_point(slant, theta)
            samples.append((project(*point), point[1]))
        width = 5 if slant in (S_MIN, S_MAX) else 2
        alpha = 235 if slant == S_MAX else 150
        if slant == S_MIN:
            alpha = 170
        for pts in split_visible_segments(samples, is_front=True):
            draw_polyline(draw, pts, (55, 92, 122, alpha), width)

    targets = []
    for _ in range(84):
        theta = 2 * math.pi * random.random()
        slant = S_MIN + (S_MAX - S_MIN) * random.random()
        p = cone_point(slant, theta)
        v = random_inward_velocity(theta)
        targets.append((p, v, theta, slant))

    targets.sort(key=lambda item: item[0][1])

    for p, v, _, _ in targets:
        px, py = project(*p)
        q = (p[0] + 0.11 * v[0], p[1] + 0.11 * v[1], p[2] + 0.11 * v[2])
        qx, qy = project(*q)
        draw_arrow(draw, (px, py), (qx - px, qy - py))

    for p, _, _, slant in targets:
        px, py = project(*p)
        depth = (p[1] / (S_MAX * math.sin(ALPHA)) + 1) / 2
        r = 6.8 + 2.5 * depth + 1.4 * slant
        fill = (230, 76, 62, int(205 + 45 * depth))
        outline = (125, 32, 32, 240)
        draw.ellipse(
            [
                int((px - r) * SCALE),
                int((py - r) * SCALE),
                int((px + r) * SCALE),
                int((py + r) * SCALE),
            ],
            fill=fill,
            outline=outline,
            width=1 * SCALE,
        )

    img = img.resize((W, H), Image.Resampling.LANCZOS)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    img.save(OUT)
    print(OUT)


if __name__ == "__main__":
    main()
