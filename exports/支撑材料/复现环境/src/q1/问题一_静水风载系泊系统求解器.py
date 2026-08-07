"""问题一静水风载系泊系统的确定性数值求解器。

从仓库根目录运行：python src/q1/q1_static_equilibrium.py
仅使用标准库；不读取或改写 data/raw。模型参数和公式口径以 model_cards/q1.md
及既有 q1_static_model.m 为准。本程序实现相同的递推与悬链线闭合关系，输出均写入
data/processed 和 outputs/q1。
"""

from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
PROCESSED = ROOT / "data" / "processed"
TABLES = ROOT / "outputs" / "q1" / "tables"
FIGURES = ROOT / "outputs" / "q1" / "figures"
LOGS = ROOT / "outputs" / "q1" / "logs"


@dataclass(frozen=True)
class Constants:
    g: float = 9.8
    rho_water: float = 1025.0
    depth: float = 18.0
    buoy_radius: float = 1.0
    buoy_height: float = 2.0
    buoy_mass: float = 1000.0
    pipe_length: float = 1.0
    pipe_radius: float = 0.025
    pipe_mass: float = 10.0
    n_pipes: int = 4
    barrel_length: float = 1.0
    barrel_radius: float = 0.15
    barrel_mass: float = 100.0
    ball_mass: float = 1200.0
    steel_density: float = 7850.0
    chain_length: float = 22.05
    chain_mass_per_length: float = 7.0

    @property
    def pipe_effective_weight(self) -> float:
        return self.pipe_mass * self.g - cylinder_buoyancy(
            self.pipe_radius, self.pipe_length, self
        )

    @property
    def barrel_effective_weight(self) -> float:
        return self.barrel_mass * self.g - cylinder_buoyancy(
            self.barrel_radius, self.barrel_length, self
        )

    @property
    def ball_effective_weight(self) -> float:
        return self.ball_mass * self.g * (1.0 - self.rho_water / self.steel_density)

    @property
    def chain_effective_weight_per_length(self) -> float:
        return self.chain_mass_per_length * self.g * (
            1.0 - self.rho_water / self.steel_density
        )


@dataclass(frozen=True)
class ChainGeometry:
    vertical_drop: float
    horizontal_span: float
    suspended_length: float
    seabed_length: float
    lower_vertical_tension: float
    lower_angle_deg: float
    regime: str


@dataclass(frozen=True)
class State:
    draft: float
    wind_speed: float
    horizontal_tension: float
    buoy_vertical_tension: float
    chain_top_vertical: float
    pipe_angles: tuple[float, ...]
    barrel_angle: float
    rigid_vertical_drop: float
    rigid_horizontal_shift: float
    chain: ChainGeometry
    model_depth: float
    radius: float


def cylinder_buoyancy(radius: float, length: float, c: Constants) -> float:
    return c.rho_water * c.g * math.pi * radius**2 * length


def member_angle(horizontal_tension: float, top_vertical_tension: float, effective_weight: float) -> float:
    denominator = top_vertical_tension - effective_weight / 2.0
    if denominator <= 0.0:
        raise ValueError("构件顶端竖向张力不足，刚体受拉平衡假设失效。")
    return math.atan2(horizontal_tension, denominator)


def chain_geometry(horizontal_tension: float, top_vertical_tension: float, c: Constants) -> ChainGeometry:
    if horizontal_tension <= 0.0 or top_vertical_tension <= 0.0:
        raise ValueError("悬链线顶端张力分量必须为正。")
    q = c.chain_effective_weight_per_length
    touchdown_length = top_vertical_tension / q
    if touchdown_length <= c.chain_length:
        lower_vertical = 0.0
        suspended_length = touchdown_length
        seabed_length = c.chain_length - suspended_length
        regime = "touchdown_with_seabed_contact"
    else:
        lower_vertical = top_vertical_tension - q * c.chain_length
        suspended_length = c.chain_length
        seabed_length = 0.0
        regime = "fully_suspended"
    top_norm = math.hypot(horizontal_tension, top_vertical_tension)
    lower_norm = math.hypot(horizontal_tension, lower_vertical)
    vertical_drop = (top_norm - lower_norm) / q
    horizontal_span = horizontal_tension / q * (
        math.asinh(top_vertical_tension / horizontal_tension)
        - math.asinh(lower_vertical / horizontal_tension)
    )
    return ChainGeometry(
        vertical_drop=vertical_drop,
        horizontal_span=horizontal_span,
        suspended_length=suspended_length,
        seabed_length=seabed_length,
        lower_vertical_tension=lower_vertical,
        lower_angle_deg=math.degrees(math.atan2(lower_vertical, horizontal_tension)),
        regime=regime,
    )


def evaluate_state(draft: float, wind_speed: float, c: Constants) -> State:
    if not 0.0 <= draft <= c.buoy_height:
        raise ValueError("浮标吃水必须在 [0, 2] m 内。")
    horizontal_tension = 0.625 * (2.0 * (c.buoy_height - draft)) * wind_speed**2
    vertical_top = (
        c.rho_water * c.g * math.pi * c.buoy_radius**2 * draft - c.buoy_mass * c.g
    )
    current_vertical = vertical_top
    pipe_angles: list[float] = []
    rigid_vertical_drop = 0.0
    rigid_horizontal_shift = 0.0
    for _ in range(c.n_pipes):
        angle = member_angle(horizontal_tension, current_vertical, c.pipe_effective_weight)
        pipe_angles.append(angle)
        rigid_vertical_drop += c.pipe_length * math.cos(angle)
        rigid_horizontal_shift += c.pipe_length * math.sin(angle)
        current_vertical -= c.pipe_effective_weight
    barrel_angle = member_angle(horizontal_tension, current_vertical, c.barrel_effective_weight)
    rigid_vertical_drop += c.barrel_length * math.cos(barrel_angle)
    rigid_horizontal_shift += c.barrel_length * math.sin(barrel_angle)
    chain_top_vertical = (
        current_vertical - c.barrel_effective_weight - c.ball_effective_weight
    )
    chain = chain_geometry(horizontal_tension, chain_top_vertical, c)
    model_depth = draft + rigid_vertical_drop + chain.vertical_drop
    radius = rigid_horizontal_shift + chain.horizontal_span + chain.seabed_length
    return State(
        draft=draft,
        wind_speed=wind_speed,
        horizontal_tension=horizontal_tension,
        buoy_vertical_tension=vertical_top,
        chain_top_vertical=chain_top_vertical,
        pipe_angles=tuple(pipe_angles),
        barrel_angle=barrel_angle,
        rigid_vertical_drop=rigid_vertical_drop,
        rigid_horizontal_shift=rigid_horizontal_shift,
        chain=chain,
        model_depth=model_depth,
        radius=radius,
    )


def solve_draft(wind_speed: float, c: Constants) -> State:
    low, high = 0.670, 0.850
    f_low = evaluate_state(low, wind_speed, c).model_depth - c.depth
    f_high = evaluate_state(high, wind_speed, c).model_depth - c.depth
    if not (f_low < 0.0 and f_high > 0.0):
        raise ValueError(f"吃水求根区间未包围解：f(low)={f_low:.6g}, f(high)={f_high:.6g}。")
    for _ in range(80):
        midpoint = (low + high) / 2.0
        if evaluate_state(midpoint, wind_speed, c).model_depth < c.depth:
            low = midpoint
        else:
            high = midpoint
    return evaluate_state((low + high) / 2.0, wind_speed, c)


def write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def state_row(state: State) -> dict[str, object]:
    angles = [math.degrees(angle) for angle in state.pipe_angles]
    return {
        "wind_speed_mps": f"{state.wind_speed:.0f}",
        "pipe_1_angle_deg": f"{angles[0]:.14g}",
        "pipe_2_angle_deg": f"{angles[1]:.14g}",
        "pipe_3_angle_deg": f"{angles[2]:.14g}",
        "pipe_4_angle_deg": f"{angles[3]:.14g}",
        "barrel_angle_deg": f"{math.degrees(state.barrel_angle):.14g}",
        "last_chain_angle_deg": f"{state.chain.lower_angle_deg:.14g}",
        "buoy_draft_m": f"{state.draft:.14g}",
        "swing_radius_m": f"{state.radius:.14g}",
        "suspended_chain_length_m": f"{state.chain.suspended_length:.14g}",
        "seabed_chain_length_m": f"{state.chain.seabed_length:.14g}",
        "depth_residual_m": f"{state.model_depth - 18.0:.14g}",
        "chain_regime": state.chain.regime,
    }


def write_summary_tex(path: Path, states: list[State], wind_speeds: list[int]) -> None:
    labels = [
        "第1节钢管与海平面竖直方向夹角（度）",
        "第2节钢管与海平面竖直方向夹角（度）",
        "第3节钢管与海平面竖直方向夹角（度）",
        "第4节钢管与海平面竖直方向夹角（度）",
        "钢桶与海平面竖直方向夹角（度）",
        "最后一条锚链与海平面水平方向夹角（度）",
        "浮标的吃水深度（m）",
        "浮标的游动半径（m）（以锚所在位置为圆心）",
    ]

    def values(state: State) -> list[float]:
        return [
            *(math.degrees(angle) for angle in state.pipe_angles),
            math.degrees(state.barrel_angle),
            state.chain.lower_angle_deg,
            state.draft,
            state.radius,
        ]

    left, right = values(states[0]), values(states[1])
    lines = [
        "% 由 src/q1/q1_static_equilibrium.py 自动生成；请勿手工修改。",
        "\\begin{tabular}{lrr}",
        "\\toprule",
        f"指标 & {wind_speeds[0]} m/s & {wind_speeds[1]} m/s \\\\ ",
        "\\midrule",
    ]
    lines.extend(f"{label} & {a:.3f} & {b:.3f} \\\\ " for label, a, b in zip(labels, left, right))
    lines.extend(["\\bottomrule", "\\end{tabular}", ""])
    path.write_text("\n".join(lines), encoding="utf-8")


def chain_rows(state: State, c: Constants) -> list[dict[str, object]]:
    chain = state.chain
    rows: list[dict[str, object]] = []
    n_suspended = 241
    top_x = state.radius - state.rigid_horizontal_shift
    top_z = -state.draft - state.rigid_vertical_drop
    for index in range(n_suspended):
        arc = chain.suspended_length * index / (n_suspended - 1)
        local_vertical = state.chain_top_vertical - c.chain_effective_weight_per_length * arc
        x_from_top = state.horizontal_tension / c.chain_effective_weight_per_length * (
            math.asinh(state.chain_top_vertical / state.horizontal_tension)
            - math.asinh(local_vertical / state.horizontal_tension)
        )
        z_drop = (
            math.hypot(state.horizontal_tension, state.chain_top_vertical)
            - math.hypot(state.horizontal_tension, local_vertical)
        ) / c.chain_effective_weight_per_length
        rows.append(
            {
                "segment": "suspended",
                "arc_length_m": f"{arc:.14g}",
                "x_m": f"{top_x - x_from_top:.14g}",
                "z_m": f"{top_z - z_drop:.14g}",
            }
        )
    if chain.seabed_length > 0.0:
        for index in range(1, 41):
            fraction = index / 40.0
            rows.append(
                {
                    "segment": "seabed",
                    "arc_length_m": f"{chain.suspended_length + chain.seabed_length * fraction:.14g}",
                    "x_m": f"{chain.seabed_length * (1.0 - fraction):.14g}",
                    "z_m": f"{-c.depth:.14g}",
                }
            )
    return rows


def svg_polyline(points: list[tuple[float, float]], color: str, width: float) -> str:
    encoded = " ".join(f"{x:.2f},{y:.2f}" for x, y in points)
    return f'<polyline fill="none" stroke="{color}" stroke-width="{width}" points="{encoded}"/>'


def write_depth_svg(path: Path, states: list[State], c: Constants) -> None:
    width, height, margin = 960, 600, 90
    x0, x1, y0, y1 = 0.670, 0.850, 6.0, 30.0

    def project(draft: float, depth: float) -> tuple[float, float]:
        x = margin + (draft - x0) / (x1 - x0) * (width - 2 * margin)
        y = height - margin - (depth - y0) / (y1 - y0) * (height - 2 * margin)
        return x, y

    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        f'<line x1="{margin}" y1="{height-margin}" x2="{width-margin}" y2="{height-margin}" stroke="black"/>',
        f'<line x1="{margin}" y1="{margin}" x2="{margin}" y2="{height-margin}" stroke="black"/>',
    ]
    for speed, color in ((12, "#1f77b4"), (24, "#c44e52")):
        points = [project(0.670 + 0.001 * index, evaluate_state(0.670 + 0.001 * index, speed, c).model_depth) for index in range(181)]
        lines.append(svg_polyline(points, color, 3.0))
    a, b = project(x0, c.depth), project(x1, c.depth)
    lines.append(f'<line x1="{a[0]:.2f}" y1="{a[1]:.2f}" x2="{b[0]:.2f}" y2="{b[1]:.2f}" stroke="#222" stroke-dasharray="8 5"/>')
    for state in states:
        x, y = project(state.draft, state.model_depth)
        lines.append(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="5" fill="black"/>')
    lines.extend([
        '<text x="480" y="40" text-anchor="middle" font-family="sans-serif" font-size="22">Model depth versus buoy draft</text>',
        '<text x="480" y="580" text-anchor="middle" font-family="sans-serif" font-size="18">Buoy draft (m)</text>',
        '<text x="28" y="300" text-anchor="middle" font-family="sans-serif" font-size="18" transform="rotate(-90 28 300)">Model depth (m)</text>',
        '<text x="690" y="85" font-family="sans-serif" font-size="16" fill="#1f77b4">12 m/s</text>',
        '<text x="690" y="110" font-family="sans-serif" font-size="16" fill="#c44e52">24 m/s</text>',
        '</svg>',
    ])
    path.write_text("\n".join(lines), encoding="utf-8")


def write_chain_svg(path: Path, states: list[State], c: Constants) -> None:
    width, height, margin = 900, 850, 85
    x0, x1, y0, y1 = -0.8, 18.8, -18.8, 1.2

    def project(x_value: float, z_value: float) -> tuple[float, float]:
        x = margin + (x_value - x0) / (x1 - x0) * (width - 2 * margin)
        y = height - margin - (z_value - y0) / (y1 - y0) * (height - 2 * margin)
        return x, y

    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        f'<line x1="{margin}" y1="{height-margin}" x2="{width-margin}" y2="{height-margin}" stroke="black"/>',
        f'<line x1="{margin}" y1="{margin}" x2="{margin}" y2="{height-margin}" stroke="black"/>',
    ]
    for state, color in zip(states, ("#1f77b4", "#c44e52")):
        rows = chain_rows(state, c)
        points = [project(float(row["x_m"]), float(row["z_m"])) for row in rows]
        lines.append(svg_polyline(points, color, 3.0))
    left, right = project(x0, 0.0), project(x1, 0.0)
    bed_left, bed_right = project(x0, -c.depth), project(x1, -c.depth)
    lines.extend([
        f'<line x1="{left[0]:.2f}" y1="{left[1]:.2f}" x2="{right[0]:.2f}" y2="{right[1]:.2f}" stroke="#4d8eb5"/>',
        f'<line x1="{bed_left[0]:.2f}" y1="{bed_left[1]:.2f}" x2="{bed_right[0]:.2f}" y2="{bed_right[1]:.2f}" stroke="#8b6a4e"/>',
        '<text x="450" y="42" text-anchor="middle" font-family="sans-serif" font-size="22">Static chain shapes</text>',
        '<text x="450" y="832" text-anchor="middle" font-family="sans-serif" font-size="18">Horizontal coordinate from anchor (m)</text>',
        '<text x="28" y="425" text-anchor="middle" font-family="sans-serif" font-size="18" transform="rotate(-90 28 425)">Elevation (m)</text>',
        '<text x="610" y="100" font-family="sans-serif" font-size="16" fill="#1f77b4">12 m/s</text>',
        '<text x="610" y="125" font-family="sans-serif" font-size="16" fill="#c44e52">24 m/s</text>',
        '</svg>',
    ])
    path.write_text("\n".join(lines), encoding="utf-8")


def write_check_log(path: Path, states: list[State], c: Constants) -> None:
    lines = [
        "Generated by src/q1/q1_static_equilibrium.py",
        f"pipe_effective_weight_N={c.pipe_effective_weight:.12f}",
        f"barrel_effective_weight_N={c.barrel_effective_weight:.12f}",
        f"ball_effective_weight_N={c.ball_effective_weight:.12f}",
        f"chain_effective_weight_N_per_m={c.chain_effective_weight_per_length:.12f}",
    ]
    for state in states:
        depth_residual = abs(state.model_depth - c.depth)
        length_residual = abs(
            state.chain.suspended_length + state.chain.seabed_length - c.chain_length
        )
        if depth_residual > 1e-10 or length_residual > 1e-10:
            raise AssertionError("水深或锚链长度闭合核验失败。")
        if not 0.0 <= state.draft <= c.buoy_height:
            raise AssertionError("浮标吃水越界。")
        lines.extend([
            "",
            f"wind_speed_mps={state.wind_speed:.0f}",
            f"draft_in_bounds={int(0.0 <= state.draft <= c.buoy_height)}",
            f"depth_residual_abs_m={depth_residual:.3e}",
            f"chain_length_residual_abs_m={length_residual:.3e}",
            f"anchor_angle_within_16deg={int(state.chain.lower_angle_deg <= 16.0)}",
            f"chain_regime={state.chain.regime}",
        ])
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    for directory in (PROCESSED, TABLES, FIGURES, LOGS):
        directory.mkdir(parents=True, exist_ok=True)
    c = Constants()
    wind_speeds = [12, 24]
    states = [solve_draft(speed, c) for speed in wind_speeds]

    curve_rows: list[dict[str, object]] = []
    for speed in wind_speeds:
        for index in range(181):
            draft = 0.670 + 0.001 * index
            state = evaluate_state(draft, speed, c)
            curve_rows.append({
                "wind_speed_mps": speed,
                "draft_m": f"{draft:.3f}",
                "model_depth_m": f"{state.model_depth:.14g}",
            })
    write_csv(
        PROCESSED / "q1_depth_vs_draft.csv",
        ["wind_speed_mps", "draft_m", "model_depth_m"],
        curve_rows,
    )
    write_csv(
        PROCESSED / "q1_depth_draft_solutions.csv",
        ["wind_speed_mps", "draft_m", "model_depth_m", "target_depth_m"],
        [{
            "wind_speed_mps": f"{state.wind_speed:.0f}",
            "draft_m": f"{state.draft:.14g}",
            "model_depth_m": f"{state.model_depth:.14g}",
            "target_depth_m": f"{c.depth:.14g}",
        } for state in states],
    )
    fieldnames = [
        "wind_speed_mps", "pipe_1_angle_deg", "pipe_2_angle_deg", "pipe_3_angle_deg",
        "pipe_4_angle_deg", "barrel_angle_deg", "last_chain_angle_deg", "buoy_draft_m",
        "swing_radius_m", "suspended_chain_length_m", "seabed_chain_length_m",
        "depth_residual_m", "chain_regime",
    ]
    write_csv(TABLES / "q1_static_results.csv", fieldnames, [state_row(state) for state in states])
    write_summary_tex(TABLES / "q1_static_results.tex", states, wind_speeds)
    for state in states:
        write_csv(
            PROCESSED / f"q1_chain_shape_{state.wind_speed:.0f}ms.csv",
            ["segment", "arc_length_m", "x_m", "z_m"],
            chain_rows(state, c),
        )
    write_depth_svg(FIGURES / "q1_depth_vs_draft.svg", states, c)
    write_chain_svg(FIGURES / "q1_chain_shapes.svg", states, c)
    write_check_log(LOGS / "q1_static_checks.txt", states, c)


if __name__ == "__main__":
    main()
