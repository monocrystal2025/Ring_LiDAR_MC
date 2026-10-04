from __future__ import annotations

import math
import random
from pathlib import Path

from PIL import Image, ImageDraw


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "assets" / "random_targets_subgraph.png"

W, H = 1600, 1000
SCALE = 3

CX = 800
BASE_Y = 760
RX = 520
RY = 170
RZ = 560


def sxy(p: tuple[float, float]) -> tuple[int, int]:
    return int(round(p[0] * SCALE)), int(round(p[1] * SCALE))


def project(x: float, y: float, z: float) -> tuple[float, float]:
    # Upright orthographic projection: x is horizontal, z is vertical, y is depth.
    # There is no x-y shear, so latitude rings stay level with the ground ellipse.
    return CX + RX * x, BASE_Y + RY * y - RZ * z


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


def hemisphere_point(theta: float, u: float) -> tuple[float, float, float]:
    phi = math.acos(u)
    return (
        math.sin(phi) * math.cos(theta),
        math.sin(phi) * math.sin(theta),
        math.cos(phi),
    )


def random_inward_velocity(p: tuple[float, float, float]) -> tuple[float, float, float]:
    theta_v = 2 * math.pi * random.random()
    u_v = 2 * random.random() - 1
    phi_v = math.acos(u_v)
    v = (
        math.sin(phi_v) * math.cos(theta_v),
        math.sin(phi_v) * math.sin(theta_v),
        math.cos(phi_v),
    )
    if p[0] * v[0] + p[1] * v[1] + p[2] * v[2] > 0:
        v = (-v[0], -v[1], -v[2])
    return v


def main():
    random.seed(42)
    img = Image.new("RGBA", (W * SCALE, H * SCALE), (255, 255, 255, 255))
    draw = ImageDraw.Draw(img, "RGBA")

    # A very light base plane keeps the dome visually upright without adding labels.
    draw.ellipse(
        [
            (CX - RX - 4) * SCALE,
            (BASE_Y - RY - 4) * SCALE,
            (CX + RX + 4) * SCALE,
            (BASE_Y + RY + 4) * SCALE,
        ],
        fill=(255, 255, 255, 255),
        outline=(70, 108, 137, 220),
        width=4 * SCALE,
    )

    # Back meridians first, so the rear half of each latitude ring reads as
    # attached to the hemisphere surface rather than floating in space.
    for theta in [7 * math.pi / 6, 4 * math.pi / 3, 3 * math.pi / 2, 5 * math.pi / 3, 11 * math.pi / 6]:
        pts = []
        for i in range(121):
            phi = (math.pi / 2) * i / 120
            x = math.sin(phi) * math.cos(theta)
            y = math.sin(phi) * math.sin(theta)
            z = math.cos(phi)
            pts.append(project(x, y, z))
        draw_polyline(draw, pts, (64, 115, 153, 70), 2)

    # Latitude rings: rear half is lighter and drawn underneath; front half is
    # stronger and drawn later. This keeps every ring visually on the dome.
    for phi, alpha, line_width in [
        (math.pi / 8, 110, 2),
        (2 * math.pi / 8, 120, 2),
        (3 * math.pi / 8, 130, 2),
        (4 * math.pi / 8, 185, 3),
    ]:
        samples = []
        for i in range(241):
            theta = 2 * math.pi * i / 240
            x = math.sin(phi) * math.cos(theta)
            y = math.sin(phi) * math.sin(theta)
            z = math.cos(phi)
            samples.append((project(x, y, z), y))
        for pts in split_visible_segments(samples, is_front=False):
            draw_polyline(draw, pts, (64, 115, 153, max(45, alpha - 70)), line_width)

    # Front meridians drawn symmetrically across the upright dome.
    for theta in [0, math.pi / 6, math.pi / 3, math.pi / 2, 2 * math.pi / 3, 5 * math.pi / 6, math.pi]:
        pts = []
        for i in range(121):
            phi = (math.pi / 2) * i / 120
            x = math.sin(phi) * math.cos(theta)
            y = math.sin(phi) * math.sin(theta)
            z = math.cos(phi)
            pts.append(project(x, y, z))
        draw_polyline(draw, pts, (64, 115, 153, 125), 2)

    for phi, alpha, line_width in [
        (math.pi / 8, 130, 2),
        (2 * math.pi / 8, 145, 2),
        (3 * math.pi / 8, 155, 2),
        (4 * math.pi / 8, 185, 3),
    ]:
        samples = []
        for i in range(241):
            theta = 2 * math.pi * i / 240
            x = math.sin(phi) * math.cos(theta)
            y = math.sin(phi) * math.sin(theta)
            z = math.cos(phi)
            samples.append((project(x, y, z), y))
        for pts in split_visible_segments(samples, is_front=True):
            draw_polyline(draw, pts, (64, 115, 153, alpha), line_width)

    # Hidden/back half of the base ellipse is lighter; front half is stronger.
    back = []
    front = []
    for i in range(181):
        t = math.pi + math.pi * i / 180
        back.append((CX + RX * math.cos(t), BASE_Y + RY * math.sin(t)))
        t2 = math.pi * i / 180
        front.append((CX + RX * math.cos(t2), BASE_Y + RY * math.sin(t2)))
    draw_polyline(draw, back, (70, 108, 137, 95), 2)
    draw_polyline(draw, front, (55, 92, 122, 235), 5)

    targets = []
    for _ in range(78):
        p = hemisphere_point(2 * math.pi * random.random(), random.random())
        v = random_inward_velocity(p)
        targets.append((p, v))

    # Draw farther targets first so nearer points and arrows remain readable.
    targets.sort(key=lambda item: item[0][1])
    for p, v in targets:
        px, py = project(*p)
        p_next = (p[0] + 0.14 * v[0], p[1] + 0.14 * v[1], p[2] + 0.14 * v[2])
        qx, qy = project(*p_next)
        draw_arrow(draw, (px, py), (qx - px, qy - py))

    for p, _ in targets:
        px, py = project(*p)
        depth = (p[1] + 1) / 2
        r = 7.0 + 2.4 * depth
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
