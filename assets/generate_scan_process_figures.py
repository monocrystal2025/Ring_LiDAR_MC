from __future__ import annotations

import math
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw


ROOT = Path(__file__).resolve().parents[1]
ASSET_DIR = ROOT / "assets"

WHITE = (255, 255, 255, 255)
BOUNDARY = (54, 91, 121, 245)
PATH = (16, 112, 158, 210)
PATH_LIGHT = (16, 112, 158, 95)
SPOT = (16, 112, 158, 235)

SCALE = 3


def sxy(p: tuple[float, float]) -> tuple[int, int]:
    return int(round(p[0] * SCALE)), int(round(p[1] * SCALE))


def make_canvas(width: int, height: int):
    img = Image.new("RGBA", (width * SCALE, height * SCALE), WHITE)
    return img, ImageDraw.Draw(img, "RGBA")


def save_image(img: Image.Image, path: Path, size: tuple[int, int]):
    img = img.resize(size, Image.Resampling.LANCZOS)
    path.parent.mkdir(parents=True, exist_ok=True)
    img.save(path)
    print(path)


def draw_polyline(draw: ImageDraw.ImageDraw, pts, fill, width=2):
    if len(pts) > 1:
        draw.line([sxy(p) for p in pts], fill=fill, width=width * SCALE, joint="curve")


def draw_circle_outline(draw: ImageDraw.ImageDraw, cx, cy, r, fill=BOUNDARY, width=5):
    draw.ellipse(
        [
            int((cx - r) * SCALE),
            int((cy - r) * SCALE),
            int((cx + r) * SCALE),
            int((cy + r) * SCALE),
        ],
        outline=fill,
        width=width * SCALE,
    )


def draw_spots(draw: ImageDraw.ImageDraw, pts, radius=4.0, fill=SPOT):
    for x, y in pts:
        draw.ellipse(
            [
                int((x - radius) * SCALE),
                int((y - radius) * SCALE),
                int((x + radius) * SCALE),
                int((y + radius) * SCALE),
            ],
            fill=fill,
        )


def subsample(arr: np.ndarray, max_count: int) -> np.ndarray:
    if len(arr) <= max_count:
        return arr
    idx = np.linspace(0, len(arr) - 1, max_count).round().astype(int)
    return arr[idx]


def lidarpath_like_mc(f: float, full_divergence: float, omega: float, scan_radius: float) -> np.ndarray:
    # Same geometric relationship as path_init.m + lidarpath.m:
    # D = 2*R*tan(fasan_D/2), then a spherical Archimedean spiral in phi/theta.
    pitch = 2 * scan_radius * math.tan(full_divergence / 2)
    endpoint_value = pitch * math.sqrt(scan_radius**2 - pitch**2 / 4) / scan_radius**2
    endpoint_value = min(max(endpoint_value, 0), 1)
    phi_lower = math.acos(endpoint_value) + 0.5 * math.asin(endpoint_value)
    phi_upper = 0.5 * math.asin(endpoint_value)
    k = pitch / (2 * math.pi * scan_radius)

    phi_grid = np.linspace(phi_lower, phi_upper, 9000)
    integrand = np.sqrt(k**2 + np.sin(phi_grid) ** 2)
    dphi = -np.diff(phi_grid)
    q_grid = np.concatenate([[0], np.cumsum(0.5 * (integrand[:-1] + integrand[1:]) * dphi / k)])

    q_samples = np.arange(0, q_grid[-1], omega / f)
    phi = np.interp(q_samples, q_grid, phi_grid)
    theta = (phi_lower - phi) / k
    beam_vec = np.column_stack(
        [np.sin(phi) * np.cos(theta), np.sin(phi) * np.sin(theta), np.cos(phi)]
    )

    # MC_func.m prepends/appends circular scans at the first and last beam radius.
    all_parts = []
    for base in [beam_vec[0]]:
        theta0 = math.atan2(base[1], base[0])
        r0 = math.hypot(base[0], base[1])
        dtheta = omega / f / r0
        theta_ring = np.arange(theta0, theta0 + 2 * math.pi, dtheta)
        all_parts.append(
            np.column_stack([r0 * np.cos(theta_ring), r0 * np.sin(theta_ring), np.full_like(theta_ring, base[2])])
        )
    all_parts.append(beam_vec)
    for base in [beam_vec[-1]]:
        theta0 = math.atan2(base[1], base[0])
        r0 = math.hypot(base[0], base[1])
        dtheta = omega / f / r0
        theta_ring = np.arange(theta0, theta0 + 2 * math.pi, dtheta)
        all_parts.append(
            np.column_stack([r0 * np.cos(theta_ring), r0 * np.sin(theta_ring), np.full_like(theta_ring, base[2])])
        )
    return np.vstack(all_parts)


def generate_roi_spiral_path(r_end: float, pitch: float, step_size: float) -> np.ndarray:
    # Port of MC_func_ROI.m generate_spiral_path.
    a = pitch / (2 * math.pi)
    theta_max = r_end / a
    length = 0.5 * a * (theta_max * math.sqrt(theta_max**2 + 1) + math.log(theta_max + math.sqrt(theta_max**2 + 1)))
    path = np.zeros((int(math.ceil(length / step_size)) + 2, 2))
    point_count = 1
    theta = 0.0
    while True:
        ds_dtheta = a * math.sqrt(theta**2 + 1)
        dtheta = step_size / ds_dtheta
        theta_new = theta + dtheta
        r_new = a * theta_new
        if r_new >= r_end:
            theta_end = r_end / a
            path[point_count] = [r_end * math.cos(theta_end), r_end * math.sin(theta_end)]
            point_count += 1
            break
        path[point_count] = [r_new * math.cos(theta_new), r_new * math.sin(theta_new)]
        point_count += 1
        theta = theta_new
    return path[:point_count]


def generate_roi_raster_path(diameter: float, step_size: float, roi_radius: float) -> np.ndarray:
    # Port of the active MC_func_ROI.m generate_path(D,d,a).
    if diameter >= 2 * roi_radius:
        y = np.arange(-roi_radius - step_size, roi_radius + step_size + 1e-9, step_size)
        return np.column_stack([np.zeros_like(y), y])

    step_count_single = math.ceil(2 * roi_radius / step_size)
    path = [[-roi_radius + diameter / 2, -roi_radius]]
    columns = math.ceil(2 * roi_radius / diameter)
    for i in range(1, columns + 1):
        x0 = -roi_radius + diameter / 2 + diameter * (i - 1)
        if i % 2 == 1:
            y0 = -roi_radius
            y = np.linspace(y0, -y0, math.ceil(-2 * y0 / step_size))
        else:
            y0 = roi_radius
            y = np.linspace(y0, -y0, math.ceil(2 * y0 / step_size))
        x = np.full(step_count_single, x0)
        path.extend(np.column_stack([x, y]).tolist())
        if i != columns:
            x_conn = np.arange(x0, x0 + ((2 * roi_radius - diameter) / (columns - 1)) + 1e-9, step_size)
            y_conn = np.full_like(x_conn, -y0)
            path.extend(np.column_stack([x_conn, y_conn]).tolist())
    return np.asarray(path)


def project_hemisphere(point: np.ndarray) -> tuple[float, float]:
    # Same upright view style used by the EA random-target figure.
    cx, base_y = 800, 760
    rx, ry, rz = 520, 170, 560
    x, y, z = point
    return cx + rx * x, base_y + ry * y - rz * z


def split_samples_by_front(samples: list[tuple[tuple[float, float], float]], is_front: bool):
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


def draw_hemisphere_grid(draw: ImageDraw.ImageDraw):
    cx, base_y = 800, 760
    rx, ry = 520, 170

    # Back meridians first, then rear/front latitude arcs, then front meridians.
    for theta in [7 * math.pi / 6, 4 * math.pi / 3, 3 * math.pi / 2, 5 * math.pi / 3, 11 * math.pi / 6]:
        pts = []
        for i in range(121):
            phi = (math.pi / 2) * i / 120
            p = np.array([math.sin(phi) * math.cos(theta), math.sin(phi) * math.sin(theta), math.cos(phi)])
            pts.append(project_hemisphere(p))
        draw_polyline(draw, pts, (64, 115, 153, 8), 1)

    for phi, alpha, line_width in [
        (math.pi / 8, 105, 2),
        (2 * math.pi / 8, 112, 2),
        (3 * math.pi / 8, 122, 2),
        (4 * math.pi / 8, 170, 3),
    ]:
        samples = []
        for i in range(241):
            theta = 2 * math.pi * i / 240
            p = np.array([math.sin(phi) * math.cos(theta), math.sin(phi) * math.sin(theta), math.cos(phi)])
            samples.append((project_hemisphere(p), p[1]))
        for pts in split_samples_by_front(samples, is_front=False):
            draw_polyline(draw, pts, (64, 115, 153, max(7, int(alpha * 0.10))), 1)

    for theta in [0, math.pi / 6, math.pi / 3, math.pi / 2, 2 * math.pi / 3, 5 * math.pi / 6, math.pi]:
        pts = []
        for i in range(121):
            phi = (math.pi / 2) * i / 120
            p = np.array([math.sin(phi) * math.cos(theta), math.sin(phi) * math.sin(theta), math.cos(phi)])
            pts.append(project_hemisphere(p))
        draw_polyline(draw, pts, (64, 115, 153, 14), 1)

    for phi, alpha, line_width in [
        (math.pi / 8, 122, 2),
        (2 * math.pi / 8, 138, 2),
        (3 * math.pi / 8, 148, 2),
        (4 * math.pi / 8, 185, 3),
    ]:
        samples = []
        for i in range(241):
            theta = 2 * math.pi * i / 240
            p = np.array([math.sin(phi) * math.cos(theta), math.sin(phi) * math.sin(theta), math.cos(phi)])
            samples.append((project_hemisphere(p), p[1]))
        for pts in split_samples_by_front(samples, is_front=True):
            draw_polyline(draw, pts, (64, 115, 153, int(alpha * 0.14)), 1)

    back = []
    front = []
    for i in range(181):
        t = math.pi + math.pi * i / 180
        back.append((cx + rx * math.cos(t), base_y + ry * math.sin(t)))
        t2 = math.pi * i / 180
        front.append((cx + rx * math.cos(t2), base_y + ry * math.sin(t2)))
    draw_polyline(draw, back, (70, 108, 137, 22), 1)
    draw_polyline(draw, front, (55, 92, 122, 235), 5)


def draw_path_on_hemisphere(
    draw: ImageDraw.ImageDraw,
    beam_vec: np.ndarray,
    include_back: bool,
    include_front: bool,
):
    line_samples = subsample(beam_vec, 2400)
    back = []
    front = []
    for i in range(len(line_samples) - 1):
        segment = line_samples[i : i + 2]
        pts = [project_hemisphere(segment[0]), project_hemisphere(segment[1])]
        if float(segment[:, 1].mean()) >= 0:
            front.append(pts)
        else:
            back.append(pts)

    if include_back:
        for pts in back:
            draw_polyline(draw, pts, (16, 112, 158, 78), 2)
    if include_front:
        for pts in front:
            draw_polyline(draw, pts, (16, 112, 158, 225), 3)

    spot_samples = subsample(beam_vec, 560)
    for p in spot_samples:
        x, y = project_hemisphere(p)
        is_front = p[1] >= 0
        if is_front and not include_front:
            continue
        if (not is_front) and not include_back:
            continue
        radius = 3.8 if is_front else 3.1
        fill = (16, 112, 158, 235 if is_front else 125)
        draw.ellipse(
            [
                int((x - radius) * SCALE),
                int((y - radius) * SCALE),
                int((x + radius) * SCALE),
                int((y + radius) * SCALE),
            ],
            fill=fill,
        )


def draw_ea_angular_spots():
    width, height = 1600, 1000
    img, draw = make_canvas(width, height)

    # Representative value from MC_EA.m fasan_D_LIST; path construction follows
    # MC_func.m exactly: leading circle + lidarpath spherical spiral + trailing circle.
    beam_vec = lidarpath_like_mc(f=5e3, full_divergence=275e-3, omega=2 * math.pi, scan_radius=2000)

    # Draw rear trajectory first, then the hemisphere, then the foreground trajectory.
    draw_path_on_hemisphere(draw, beam_vec, include_back=True, include_front=False)
    draw_hemisphere_grid(draw)
    draw_path_on_hemisphere(draw, beam_vec, include_back=False, include_front=True)

    save_image(img, ASSET_DIR / "ea_scan_angular_spots.png", (width, height))


def draw_roi_spiral_scan():
    width, height = 1200, 1200
    img, draw = make_canvas(width, height)
    cx = cy = width / 2

    # MC_ROI.m relationships: a=z*tand(jiaodu), D=2*z*tan(fasan_D/2),
    # step_size=omega*z/f. fasan_D is chosen from the code's sweep range
    # but intentionally not the densest case.
    z = 500.0
    roi_radius_m = z * math.tan(math.radians(25))
    diameter_m = 2 * z * math.tan(80e-3 / 2)
    step_size = 2 * math.pi * z / 5e3
    r_end = roi_radius_m + 0.5 * diameter_m
    path = generate_roi_spiral_path(r_end, diameter_m, step_size)

    scale = 430 / r_end
    roi_r = roi_radius_m * scale
    pts = np.column_stack([cx + scale * path[:, 0], cy + scale * path[:, 1]])

    draw_circle_outline(draw, cx, cy, roi_r, BOUNDARY, 5)
    draw_polyline(draw, subsample(pts, 1600).tolist(), PATH_LIGHT, 2)
    draw_spots(draw, subsample(pts, 430).tolist(), radius=3.6, fill=SPOT)
    save_image(img, ASSET_DIR / "roi_equidistant_spiral_scan.png", (width, height))


def draw_roi_raster_scan():
    width, height = 1200, 1200
    img, draw = make_canvas(width, height)
    cx = cy = width / 2

    z = 500.0
    roi_radius_m = z * math.tan(math.radians(25))
    diameter_m = 2 * z * math.tan(80e-3 / 2)
    step_size = 2 * math.pi * z / 5e3
    path = generate_roi_raster_path(diameter_m, step_size, roi_radius_m)

    scale = 430 / roi_radius_m
    pts = np.column_stack([cx + scale * path[:, 0], cy + scale * path[:, 1]])

    draw_circle_outline(draw, cx, cy, 430, BOUNDARY, 5)
    draw_polyline(draw, subsample(pts, 1400).tolist(), PATH_LIGHT, 2)
    draw_spots(draw, subsample(pts, 420).tolist(), radius=3.6, fill=SPOT)
    save_image(img, ASSET_DIR / "roi_raster_scan.png", (width, height))


def main():
    draw_ea_angular_spots()
    draw_roi_spiral_scan()
    draw_roi_raster_scan()


if __name__ == "__main__":
    main()
