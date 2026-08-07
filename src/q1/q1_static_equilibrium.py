"""Reproducible static-equilibrium calculation for CUMCM 2016 A, question 1.

The equations and conditional pipe-buoyancy input are documented in
model_cards/q1.md and docs/open_questions.md.  Run from any working directory:
    python src/q1/q1_static_equilibrium.py
"""

from __future__ import annotations

import csv
import math
from pathlib import Path
from xml.sax.saxutils import escape


ROOT = Path(__file__).resolve().parents[2]
OUTPUT = ROOT / "outputs" / "q1"
TABLES = OUTPUT / "tables"
FIGURES = OUTPUT / "figures"
LOGS = OUTPUT / "logs"


def build_parameters() -> dict[str, float | str]:
    """Read II-chain mass per length and construct the fixed Q1 parameter set."""
    chain_csv = ROOT / "problem" / "original_problem" / "附表-锚链型号和参数表.csv"
    with chain_csv.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    selected = next(row for row in rows if row["型号"] == "II")

    p: dict[str, float | str] = {
        "rho": 1025.0,
        "g": 9.8,
        "float_radius": 1.0,
        "float_height": 2.0,
        "float_diameter": 2.0,
        "float_mass": 1000.0,
        "pipe_length": 1.0,
        "pipe_radius": 0.025,
        "pipe_mass": 10.0,
        "barrel_length": 1.0,
        "barrel_radius": 0.15,
        "barrel_mass": 100.0,
        "ball_mass": 1200.0,
        "steel_density": 7850.0,
        "chain_length": 22.05,
        "chain_type": "II",
        "chain_mass_per_length": float(selected["单位长度的质量(kg/m)"]),
    }
    # Conditional input OQ-1: fully displacing cylinder using the stated outer diameter.
    p["pipe_displaced_volume"] = math.pi * float(p["pipe_radius"]) ** 2 * float(p["pipe_length"])
    p["barrel_displaced_volume"] = math.pi * float(p["barrel_radius"]) ** 2 * float(p["barrel_length"])
    p["ball_displaced_volume"] = float(p["ball_mass"]) / float(p["steel_density"])
    p["pipe_effective_weight"] = float(p["g"]) * (
        float(p["pipe_mass"]) - float(p["rho"]) * float(p["pipe_displaced_volume"])
    )
    p["barrel_effective_weight"] = float(p["g"]) * (
        float(p["barrel_mass"]) - float(p["rho"]) * float(p["barrel_displaced_volume"])
    )
    p["ball_effective_weight"] = float(p["g"]) * (
        float(p["ball_mass"]) - float(p["rho"]) * float(p["ball_displaced_volume"])
    )
    p["chain_effective_weight_per_length"] = float(p["g"]) * float(p["chain_mass_per_length"]) * (
        1.0 - float(p["rho"]) / float(p["steel_density"])
    )
    assert all(float(p[key]) > 0 for key in (
        "pipe_effective_weight", "barrel_effective_weight", "ball_effective_weight",
        "chain_effective_weight_per_length",
    ))
    return p


def chain_geometry(horizontal: float, vertical_top: float, q: float, total_length: float) -> dict[str, float | str]:
    """Catenary geometry; switches correctly between grounded and fully suspended chain."""
    assert horizontal > 0 and vertical_top > 0 and q > 0 and total_length > 0
    vertical_capacity = q * total_length
    if vertical_top <= vertical_capacity:
        suspended_length = vertical_top / q
        bed_length = total_length - suspended_length
        anchor_vertical = 0.0
        vertical_span = (math.hypot(horizontal, vertical_top) - horizontal) / q
        suspended_horizontal = horizontal / q * math.asinh(vertical_top / horizontal)
        horizontal_span = bed_length + suspended_horizontal
        mode = "partly_grounded"
    else:
        suspended_length = total_length
        bed_length = 0.0
        anchor_vertical = vertical_top - vertical_capacity
        vertical_span = (math.hypot(horizontal, vertical_top) - math.hypot(horizontal, anchor_vertical)) / q
        horizontal_span = horizontal / q * (
            math.asinh(vertical_top / horizontal) - math.asinh(anchor_vertical / horizontal)
        )
        mode = "fully_suspended"
    return {
        "mode": mode,
        "suspended_length": suspended_length,
        "bed_length": bed_length,
        "vertical_span": vertical_span,
        "horizontal_span": horizontal_span,
        "anchor_vertical": anchor_vertical,
        "anchor_angle_rad": math.atan(anchor_vertical / horizontal),
        "arc_length_residual": abs(suspended_length - (vertical_top - anchor_vertical) / q),
    }


def geometry_at_draft(p: dict[str, float | str], wind_speed: float, draft: float) -> dict[str, float | str | list[float] | dict[str, float | str]]:
    """Evaluate the force recursion and geometric closure at one draft value."""
    assert 0.0 < draft < float(p["float_height"])
    horizontal = 0.625 * (float(p["float_height"]) - draft) * float(p["float_diameter"]) * wind_speed ** 2
    vertical_float = float(p["rho"]) * float(p["g"]) * math.pi * float(p["float_radius"]) ** 2 * draft - float(p["float_mass"]) * float(p["g"])

    pipe_angles: list[float] = []
    vertical_top = vertical_float
    pipe_weight = float(p["pipe_effective_weight"])
    for _ in range(4):
        denominator = vertical_top - pipe_weight / 2.0
        assert denominator > 0, "Pipe average vertical tension must be positive."
        pipe_angles.append(math.atan(horizontal / denominator))
        vertical_top -= pipe_weight

    barrel_vertical_top = vertical_top
    barrel_denominator = barrel_vertical_top - float(p["barrel_effective_weight"]) / 2.0
    assert barrel_denominator > 0, "Barrel average vertical tension must be positive."
    barrel_angle = math.atan(horizontal / barrel_denominator)
    chain_vertical_top = barrel_vertical_top - float(p["barrel_effective_weight"]) - float(p["ball_effective_weight"])
    assert chain_vertical_top > 0, "Chain-top vertical tension must be positive."
    chain = chain_geometry(horizontal, chain_vertical_top, float(p["chain_effective_weight_per_length"]), float(p["chain_length"]))

    water_depth = draft + float(p["pipe_length"]) * sum(math.cos(angle) for angle in pipe_angles) + float(p["barrel_length"]) * math.cos(barrel_angle) + float(chain["vertical_span"])
    swim_radius = float(chain["horizontal_span"]) + float(p["pipe_length"]) * sum(math.sin(angle) for angle in pipe_angles) + float(p["barrel_length"]) * math.sin(barrel_angle)
    return {
        "wind_speed": wind_speed,
        "draft": draft,
        "horizontal_tension": horizontal,
        "vertical_at_float": vertical_float,
        "pipe_angle_rad": pipe_angles,
        "pipe_angle_deg": [math.degrees(angle) for angle in pipe_angles],
        "barrel_angle_rad": barrel_angle,
        "barrel_angle_deg": math.degrees(barrel_angle),
        "chain_anchor_angle_rad": float(chain["anchor_angle_rad"]),
        "chain_anchor_angle_deg": math.degrees(float(chain["anchor_angle_rad"])),
        "chain": chain,
        "water_depth": water_depth,
        "swim_radius": swim_radius,
    }


def solve_condition(p: dict[str, float | str], wind_speed: float, target_depth: float) -> dict[str, float | str | list[float] | dict[str, float | str]]:
    """Bisection solution of H(delta; v)=target_depth on the physical draft interval."""
    lower = (
        float(p["float_mass"]) * float(p["g"]) + 4 * float(p["pipe_effective_weight"])
        + float(p["barrel_effective_weight"]) + float(p["ball_effective_weight"]) + 1e-6
    ) / (float(p["rho"]) * float(p["g"]) * math.pi * float(p["float_radius"]) ** 2)
    upper = float(p["float_height"]) - 1e-6
    assert lower < upper
    lo, hi = lower, upper
    assert geometry_at_draft(p, wind_speed, lo)["water_depth"] < target_depth
    assert geometry_at_draft(p, wind_speed, hi)["water_depth"] > target_depth
    for _ in range(100):
        mid = (lo + hi) / 2.0
        if float(geometry_at_draft(p, wind_speed, mid)["water_depth"]) < target_depth:
            lo = mid
        else:
            hi = mid
    result = geometry_at_draft(p, wind_speed, (lo + hi) / 2.0)
    result["depth_residual"] = float(result["water_depth"]) - target_depth
    assert abs(float(result["depth_residual"])) < 1e-10
    assert float(result["chain"]["arc_length_residual"]) < 1e-10  # type: ignore[index]
    return result


def write_csv(filename: Path, header: list[str], rows: list[list[object]]) -> None:
    with filename.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        writer.writerows(rows)


def write_latex(filename: Path, column_spec: str, header: str, rows: list[str]) -> None:
    filename.write_text(
        "% Generated by src/q1/q1_static_equilibrium.py\n"
        f"\\begin{{tabular}}{{{column_spec}}}\n\\toprule\n{header} \\\\\n\\midrule\n"
        + "\n".join(rows)
        + "\n\\bottomrule\n\\end{tabular}\n",
        encoding="utf-8",
    )


def write_tables_and_figures(p: dict[str, float | str], results: list[dict[str, object]], target_depth: float) -> None:
    for folder in (TABLES, FIGURES, LOGS):
        folder.mkdir(parents=True, exist_ok=True)

    parameter_rows = [
        ["sea_water_density", f"{float(p['rho']):.9f}", "kg/m^3"],
        ["gravity", f"{float(p['g']):.9f}", "m/s^2"],
        ["II_chain_total_length", f"{float(p['chain_length']):.9f}", "m"],
        ["II_chain_mass_per_length", f"{float(p['chain_mass_per_length']):.9f}", "kg/m"],
        ["pipe_effective_weight", f"{float(p['pipe_effective_weight']):.9f}", "N"],
        ["barrel_effective_weight", f"{float(p['barrel_effective_weight']):.9f}", "N"],
        ["ball_effective_weight", f"{float(p['ball_effective_weight']):.9f}", "N"],
        ["chain_effective_weight_per_length", f"{float(p['chain_effective_weight_per_length']):.9f}", "N/m"],
    ]
    write_csv(TABLES / "q1_input_parameters.csv", ["parameter", "value", "unit"], parameter_rows)
    write_latex(TABLES / "q1_input_parameters.tex", "lrr", "Parameter & Value & Unit", [f"{a} & {b} & {c} \\\\" for a, b, c in parameter_rows])

    r12, r24 = results
    metrics = [
        "第1节钢管与竖直方向夹角（度）", "第2节钢管与竖直方向夹角（度）",
        "第3节钢管与竖直方向夹角（度）", "第4节钢管与竖直方向夹角（度）",
        "钢桶与竖直方向夹角（度）", "锚端锚链与水平方向夹角（度）",
        "浮标的吃水深度（m）", "浮标的游动半径（m）",
    ]
    values12 = list(r12["pipe_angle_deg"]) + [r12["barrel_angle_deg"], r12["chain_anchor_angle_deg"], r12["draft"], r12["swim_radius"]]
    values24 = list(r24["pipe_angle_deg"]) + [r24["barrel_angle_deg"], r24["chain_anchor_angle_deg"], r24["draft"], r24["swim_radius"]]
    summary_rows = [[metric, f"{float(v12):.6f}", f"{float(v24):.6f}"] for metric, v12, v24 in zip(metrics, values12, values24)]
    write_csv(TABLES / "q1_summary.csv", ["metric", "wind_12_mps", "wind_24_mps"], summary_rows)
    write_latex(TABLES / "q1_summary.tex", "lrr", "指标 & 12 m/s & 24 m/s", [f"{a} & {b} & {c} \\\\" for a, b, c in summary_rows])

    minimum_draft = (
        float(p["float_mass"]) * float(p["g"]) + 4 * float(p["pipe_effective_weight"])
        + float(p["barrel_effective_weight"]) + float(p["ball_effective_weight"])
    ) / (float(p["rho"]) * float(p["g"]) * math.pi * float(p["float_radius"]) ** 2)
    drafts = [minimum_draft + 1e-4 + 0.001 * i for i in range(int((1.20 - minimum_draft - 1e-4) / 0.001) + 1)]
    curve_rows: list[list[object]] = []
    depth12, depth24 = [], []
    for draft in drafts:
        h12 = float(geometry_at_draft(p, 12.0, draft)["water_depth"])
        h24 = float(geometry_at_draft(p, 24.0, draft)["water_depth"])
        depth12.append(h12)
        depth24.append(h24)
        curve_rows.append([f"{draft:.6f}", f"{h12:.6f}", f"{h24:.6f}"])
    write_csv(TABLES / "q1_depth_draft_curve.csv", ["draft_m", "depth_at_12_mps_m", "depth_at_24_mps_m"], curve_rows)
    write_latex(TABLES / "q1_depth_draft_curve.tex", "rrr", "$\\delta$ (m) & $H_{12}$ (m) & $H_{24}$ (m)", [f"{a} & {b} & {c} \\\\" for a, b, c in curve_rows])

    write_svg_chart(
        FIGURES / "q1_depth_draft_curve.svg",
        [(drafts, depth12, "#1f77b4", "v = 12 m/s"), (drafts, depth24, "#d62728", "v = 24 m/s")],
        "Draft δ (m)", "Computed water depth H(δ) (m)", "Depth--draft curves under static wind loading",
        reference_y=(target_depth, "H = 18 m"),
        markers=[(float(r12["draft"]), target_depth, "#1f77b4", "circle"), (float(r24["draft"]), target_depth, "#d62728", "square")],
    )

    max_x = 0.0
    chain_series = []
    for result, label in zip(results, ["v = 12 m/s", "v = 24 m/s"]):
        x, z = chain_profile(float(result["horizontal_tension"]), float(result["chain"]["anchor_vertical"]) + float(p["chain_effective_weight_per_length"]) * float(result["chain"]["suspended_length"]), float(p["chain_effective_weight_per_length"]), float(p["chain_length"]))  # type: ignore[index]
        color = "#1f77b4" if label.startswith("v = 12") else "#d62728"
        chain_series.append((x, z, color, label))
        max_x = max(max_x, max(x))
    write_svg_chart(
        FIGURES / "q1_chain_shapes.svg", chain_series,
        "Horizontal distance from anchor (m)", "Height above seabed (m)", "Computed catenary-chain shapes",
        reference_y=(0.0, "Seabed"), x_limits=(0.0, max_x + 1.0), y_limits=(-0.25, 14.0),
    )

    lines = [
        "Q1 static equilibrium run",
        f"Target water depth: {target_depth:.6f} m",
        "Pipe buoyancy condition: fully displacing outer cylinder, radius 0.025 m",
    ]
    for result in results:
        chain = result["chain"]  # type: ignore[assignment]
        lines += [
            f"\nWind {float(result['wind_speed']):.0f} m/s",
            f"draft={float(result['draft']):.12f}, H={float(result['water_depth']):.12f}, residual={float(result['depth_residual']):.3e}",
            f"chain_mode={chain['mode']}, suspended_length={float(chain['suspended_length']):.12f}, bed_length={float(chain['bed_length']):.12f}",
            f"arc_length_residual={float(chain['arc_length_residual']):.3e}",
        ]
    (LOGS / "q1_run_summary.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")


def chain_profile(horizontal: float, vertical_top: float, q: float, total_length: float) -> tuple[list[float], list[float]]:
    """Return chain centerline coordinates, including any grounded segment."""
    chain = chain_geometry(horizontal, vertical_top, q, total_length)
    bed_length = float(chain["bed_length"])
    anchor_vertical = float(chain["anchor_vertical"])
    if bed_length > 0:
        x = [bed_length * i / 39.0 for i in range(40)]
        z = [0.0] * 40
    else:
        x, z = [0.0], [0.0]
    vertical_top = float(vertical_top)
    for i in range(200):
        vertical = anchor_vertical + (vertical_top - anchor_vertical) * i / 199.0
        x.append(bed_length + horizontal / q * (math.asinh(vertical / horizontal) - math.asinh(anchor_vertical / horizontal)))
        z.append((math.hypot(horizontal, vertical) - math.hypot(horizontal, anchor_vertical)) / q)
    return x, z


def write_svg_chart(
    filename: Path,
    series: list[tuple[list[float], list[float], str, str]],
    x_label: str,
    y_label: str,
    title: str,
    reference_y: tuple[float, str] | None = None,
    markers: list[tuple[float, float, str, str]] | None = None,
    x_limits: tuple[float, float] | None = None,
    y_limits: tuple[float, float] | None = None,
) -> None:
    """Write a self-contained publication-sized SVG chart without external packages."""
    width, height = 960, 640
    left, right, top, bottom = 120, 45, 70, 100
    all_x = [value for xs, _, _, _ in series for value in xs]
    all_y = [value for _, ys, _, _ in series for value in ys]
    if reference_y is not None:
        all_y.append(reference_y[0])
    xmin, xmax = x_limits if x_limits else (min(all_x), max(all_x))
    ymin, ymax = y_limits if y_limits else (min(all_y), max(all_y))
    xpad = (xmax - xmin) * 0.04 or 1.0
    ypad = (ymax - ymin) * 0.08 or 1.0
    if x_limits is None:
        xmin, xmax = xmin - xpad, xmax + xpad
    if y_limits is None:
        ymin, ymax = ymin - ypad, ymax + ypad
    plot_width, plot_height = width - left - right, height - top - bottom
    sx = lambda x: left + (x - xmin) / (xmax - xmin) * plot_width
    sy = lambda y: top + (ymax - y) / (ymax - ymin) * plot_height
    lines = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        '<style>text{font-family:Arial,sans-serif;fill:#111}.tick{font-size:14px}.label{font-size:18px}.title{font-size:22px;font-weight:bold}.legend{font-size:15px}</style>',
        f'<text class="title" x="{width / 2:.1f}" y="36" text-anchor="middle">{escape(title)}</text>',
    ]
    for i in range(6):
        x = xmin + (xmax - xmin) * i / 5
        px = sx(x)
        lines += [f'<line x1="{px:.2f}" y1="{top}" x2="{px:.2f}" y2="{height - bottom}" stroke="#dddddd"/>',
                  f'<text class="tick" x="{px:.2f}" y="{height - bottom + 25}" text-anchor="middle">{x:.2f}</text>']
        y = ymin + (ymax - ymin) * i / 5
        py = sy(y)
        lines += [f'<line x1="{left}" y1="{py:.2f}" x2="{width - right}" y2="{py:.2f}" stroke="#dddddd"/>',
                  f'<text class="tick" x="{left - 12}" y="{py + 5:.2f}" text-anchor="end">{y:.2f}</text>']
    lines += [f'<line x1="{left}" y1="{height - bottom}" x2="{width - right}" y2="{height - bottom}" stroke="black" stroke-width="1.2"/>',
              f'<line x1="{left}" y1="{top}" x2="{left}" y2="{height - bottom}" stroke="black" stroke-width="1.2"/>',
              f'<text class="label" x="{width / 2:.1f}" y="{height - 35}" text-anchor="middle">{escape(x_label)}</text>',
              f'<text class="label" transform="translate(30,{height / 2:.1f}) rotate(-90)" text-anchor="middle">{escape(y_label)}</text>']
    for xs, ys, color, _ in series:
        points = " ".join(f"{sx(x):.2f},{sy(y):.2f}" for x, y in zip(xs, ys))
        lines.append(f'<polyline points="{points}" fill="none" stroke="{color}" stroke-width="2.4"/>')
    if reference_y is not None:
        py = sy(reference_y[0])
        lines.append(f'<line x1="{left}" y1="{py:.2f}" x2="{width - right}" y2="{py:.2f}" stroke="#111" stroke-width="1.4" stroke-dasharray="7,5"/>')
    for x, y, color, shape in markers or []:
        px, py = sx(x), sy(y)
        if shape == "circle":
            lines.append(f'<circle cx="{px:.2f}" cy="{py:.2f}" r="5" fill="white" stroke="{color}" stroke-width="2.5"/>')
        else:
            lines.append(f'<rect x="{px - 5:.2f}" y="{py - 5:.2f}" width="10" height="10" fill="white" stroke="{color}" stroke-width="2.5"/>')
    legend_entries = [(color, label, False) for _, _, color, label in series]
    if reference_y is not None:
        legend_entries.append(("#111", reference_y[1], True))
    for index, (color, label, dashed) in enumerate(legend_entries):
        y = top + 24 + index * 23
        dash = ' stroke-dasharray="7,5"' if dashed else ''
        lines += [f'<line x1="{left + 15}" y1="{y}" x2="{left + 45}" y2="{y}" stroke="{color}" stroke-width="2.4"{dash}/>',
                  f'<text class="legend" x="{left + 55}" y="{y + 5}">{escape(label)}</text>']
    lines.append('</svg>')
    filename.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    p = build_parameters()
    target_depth = 18.0
    results = [solve_condition(p, wind, target_depth) for wind in (12.0, 24.0)]
    for result in results:
        result["depth_residual"] = float(result["water_depth"]) - target_depth
    write_tables_and_figures(p, results, target_depth)


if __name__ == "__main__":
    main()
