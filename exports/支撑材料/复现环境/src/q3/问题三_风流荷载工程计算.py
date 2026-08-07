"""问题三工程计算：风流荷载下的分段系泊静力平衡。

运行：使用已安装 Matplotlib 的 Python 执行
      python src/q3/q3_current_engineering.py

本脚本依赖 Matplotlib 和标准库，是供建模复核的工程实现，不修改论文或模型卡。它将建模手给出的
递推式应用到 4 节钢管、钢桶、重物球和离散锚链。角度约定：刚体构件
相对竖直方向计，锚链相对水平海床计。所有 CSV 以 UTF-8 写出。
"""

from __future__ import annotations

import csv
import io
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from matplotlib.ticker import MultipleLocator


ROOT = Path(__file__).resolve().parents[2]
PROCESSED = ROOT / "data" / "processed"
TABLES = ROOT / "outputs" / "q3" / "tables"
FIGURES = ROOT / "outputs" / "q3" / "figures"
LOGS = ROOT / "outputs" / "q3" / "logs"

CHINESE_FONT = FontProperties(fname="C:/Windows/Fonts/simsun.ttc")
plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times New Roman", "SimSun"],
    "mathtext.fontset": "stix",
    "axes.unicode_minus": False,
    "axes.linewidth": 0.6,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
})


@dataclass(frozen=True)
class ChainType:
    name: str
    link_length_m: float
    mass_per_length_kg_m: float


CHAIN_TYPES = (
    ChainType("I", 0.078, 3.2),
    ChainType("II", 0.105, 7.0),
    ChainType("III", 0.120, 12.5),
    ChainType("IV", 0.150, 19.5),
    ChainType("V", 0.180, 28.12),
)


@dataclass(frozen=True)
class Constants:
    g: float = 9.8
    rho_water: float = 1025.0
    rho_steel: float = 7850.0
    buoy_radius_m: float = 1.0
    buoy_height_m: float = 2.0
    buoy_mass_kg: float = 1000.0
    pipe_length_m: float = 1.0
    pipe_diameter_m: float = 0.05
    pipe_mass_kg: float = 10.0
    n_pipes: int = 4
    barrel_length_m: float = 1.0
    barrel_diameter_m: float = 0.30
    barrel_mass_kg: float = 100.0
    chain_diameter_m: float = 0.0175

    @property
    def pipe_effective_weight_n(self) -> float:
        return self.pipe_mass_kg * self.g - self.rho_water * self.g * math.pi * (self.pipe_diameter_m / 2) ** 2 * self.pipe_length_m

    @property
    def barrel_effective_weight_n(self) -> float:
        return self.barrel_mass_kg * self.g - self.rho_water * self.g * math.pi * (self.barrel_diameter_m / 2) ** 2 * self.barrel_length_m


@dataclass(frozen=True)
class Scenario:
    depth_m: float
    wind_speed_mps: float
    current_speed_mps: float
    ball_mass_kg: float
    chain: ChainType
    link_count: int


@dataclass(frozen=True)
class MemberState:
    angle_vertical_rad: float
    current_force_n: float
    horizontal_bottom_n: float
    vertical_bottom_n: float


@dataclass(frozen=True)
class ChainState:
    horizontal_span_m: float
    vertical_drop_m: float
    suspended_length_m: float
    seabed_length_m: float
    anchor_angle_deg: float
    lower_horizontal_tension_n: float
    lower_vertical_tension_n: float
    coordinates_from_top: tuple[tuple[float, float], ...]
    regime: str


@dataclass(frozen=True)
class State:
    scenario: Scenario
    draft_m: float
    pipe_angles_rad: tuple[float, ...]
    barrel_angle_rad: float
    buoy_current_force_n: float
    ball_current_force_n: float
    horizontal_tension_at_buoy_n: float
    chain_top_horizontal_tension_n: float
    chain_top_vertical_tension_n: float
    rigid_horizontal_shift_m: float
    rigid_vertical_drop_m: float
    chain: ChainState
    model_depth_m: float
    swing_radius_m: float


def cylinder_buoyancy(diameter_m: float, length_m: float, c: Constants) -> float:
    return c.rho_water * c.g * math.pi * (diameter_m / 2) ** 2 * length_m


def ball_effective_weight(ball_mass_kg: float, c: Constants) -> float:
    return ball_mass_kg * c.g * (1.0 - c.rho_water / c.rho_steel)


def ball_diameter(ball_mass_kg: float, c: Constants) -> float:
    return (6.0 * ball_mass_kg / (math.pi * c.rho_steel)) ** (1.0 / 3.0)


def bisect_root(function: Callable[[float], float], low: float, high: float) -> float:
    f_low = function(low)
    f_high = function(high)
    if f_low == 0:
        return low
    if f_high == 0:
        return high
    if f_low * f_high > 0:
        raise ValueError("角度方程未在给定区间内包围根")
    for _ in range(90):
        middle = (low + high) / 2.0
        f_middle = function(middle)
        if f_low * f_middle <= 0:
            high = middle
        else:
            low, f_low = middle, f_middle
    return (low + high) / 2.0


def solve_rigid_member(
    horizontal_top_n: float,
    vertical_top_n: float,
    effective_weight_n: float,
    diameter_m: float,
    length_m: float,
    current_speed_mps: float,
) -> MemberState:
    """以构件相对竖直的夹角求解手稿中的中点力矩平衡。"""
    vertical_mid_n = vertical_top_n - effective_weight_n / 2.0
    if vertical_mid_n <= 0:
        raise ValueError("刚体构件中点的竖向张力非正，受拉平衡不成立")
    coefficient = 374.0 * diameter_m * length_m * current_speed_mps**2

    def residual(angle: float) -> float:
        force = coefficient * math.cos(angle)
        return math.tan(angle) - (horizontal_top_n + force / 2.0) / vertical_mid_n

    angle = bisect_root(residual, 0.0, math.pi / 2.0 - 1e-10)
    current_force_n = coefficient * math.cos(angle)
    return MemberState(
        angle_vertical_rad=angle,
        current_force_n=current_force_n,
        horizontal_bottom_n=horizontal_top_n + current_force_n,
        vertical_bottom_n=vertical_top_n - effective_weight_n,
    )


def solve_chain(
    horizontal_top_n: float,
    vertical_top_n: float,
    chain: ChainType,
    link_count: int,
    current_speed_mps: float,
    c: Constants,
) -> ChainState:
    """离散链节递推；接触海床后按既有静水口径平铺且不再施加流阻。"""
    if horizontal_top_n <= 0 or vertical_top_n <= 0:
        raise ValueError("锚链顶端张力分量必须为正")
    effective_weight_per_m = chain.mass_per_length_kg_m * c.g * (1.0 - c.rho_water / c.rho_steel)
    horizontal_n = horizontal_top_n
    vertical_n = vertical_top_n
    x_from_top = 0.0
    z_drop = 0.0
    suspended_m = 0.0
    seabed_m = 0.0
    points: list[tuple[float, float]] = [(0.0, 0.0)]
    fully_laid = False
    for _ in range(link_count):
        length_m = chain.link_length_m
        if fully_laid:
            x_from_top += length_m
            seabed_m += length_m
            points.append((x_from_top, z_drop))
            continue
        submerged_length_m = min(length_m, vertical_n / effective_weight_per_m)
        vertical_mid_n = vertical_n - effective_weight_per_m * submerged_length_m / 2.0
        if submerged_length_m <= 0 or vertical_mid_n < -1e-9:
            raise ValueError("锚链离散递推遇到非正竖向张力")
        coefficient = 374.0 * c.chain_diameter_m * submerged_length_m * current_speed_mps**2

        def residual(angle: float) -> float:
            force = coefficient * math.sin(angle)
            return math.tan(angle) - vertical_mid_n / (horizontal_n + force / 2.0)

        angle_from_horizontal = bisect_root(residual, 0.0, math.pi / 2.0 - 1e-10)
        current_force_n = coefficient * math.sin(angle_from_horizontal)
        x_from_top += submerged_length_m * math.cos(angle_from_horizontal)
        z_drop += submerged_length_m * math.sin(angle_from_horizontal)
        horizontal_n += current_force_n
        vertical_n -= effective_weight_per_m * submerged_length_m
        suspended_m += submerged_length_m
        if submerged_length_m < length_m - 1e-12:
            laid_piece_m = length_m - submerged_length_m
            x_from_top += laid_piece_m
            seabed_m += laid_piece_m
            vertical_n = 0.0
            fully_laid = True
        points.append((x_from_top, z_drop))

    total_length_m = chain.link_length_m * link_count
    if abs(suspended_m + seabed_m - total_length_m) > 1e-9:
        raise AssertionError("锚链长度守恒校验失败")
    if seabed_m > 1e-10:
        anchor_angle_deg = 0.0
        regime = "touchdown_with_seabed_contact"
    else:
        anchor_angle_deg = math.degrees(math.atan2(vertical_n, horizontal_n))
        regime = "fully_suspended"
    return ChainState(
        horizontal_span_m=x_from_top - seabed_m,
        vertical_drop_m=z_drop,
        suspended_length_m=suspended_m,
        seabed_length_m=seabed_m,
        anchor_angle_deg=anchor_angle_deg,
        lower_horizontal_tension_n=horizontal_n,
        lower_vertical_tension_n=vertical_n,
        coordinates_from_top=tuple(points),
        regime=regime,
    )


def evaluate_state(draft_m: float, scenario: Scenario, c: Constants) -> State:
    if not 0.0 < draft_m < c.buoy_height_m:
        raise ValueError("浮标吃水必须位于 (0, 2) m")
    buoyancy_n = c.rho_water * c.g * math.pi * c.buoy_radius_m**2 * draft_m
    vertical_top_n = buoyancy_n - c.buoy_mass_kg * c.g
    if vertical_top_n <= 0:
        raise ValueError("浮标提供的竖向张力非正")
    wind_force_n = 0.625 * (2.0 * (c.buoy_height_m - draft_m)) * scenario.wind_speed_mps**2
    buoy_current_force_n = 374.0 * (2.0 * draft_m) * scenario.current_speed_mps**2
    horizontal_n = wind_force_n + buoy_current_force_n
    rigid_vertical_m = 0.0
    rigid_horizontal_m = 0.0
    pipe_angles: list[float] = []
    for _ in range(c.n_pipes):
        member = solve_rigid_member(
            horizontal_n, vertical_top_n, c.pipe_effective_weight_n,
            c.pipe_diameter_m, c.pipe_length_m, scenario.current_speed_mps,
        )
        pipe_angles.append(member.angle_vertical_rad)
        rigid_vertical_m += c.pipe_length_m * math.cos(member.angle_vertical_rad)
        rigid_horizontal_m += c.pipe_length_m * math.sin(member.angle_vertical_rad)
        horizontal_n, vertical_top_n = member.horizontal_bottom_n, member.vertical_bottom_n
    barrel = solve_rigid_member(
        horizontal_n, vertical_top_n, c.barrel_effective_weight_n,
        c.barrel_diameter_m, c.barrel_length_m, scenario.current_speed_mps,
    )
    rigid_vertical_m += c.barrel_length_m * math.cos(barrel.angle_vertical_rad)
    rigid_horizontal_m += c.barrel_length_m * math.sin(barrel.angle_vertical_rad)
    ball_diameter_m = ball_diameter(scenario.ball_mass_kg, c)
    ball_current_force_n = 374.0 * (math.pi * ball_diameter_m**2 / 4.0) * scenario.current_speed_mps**2
    chain_top_horizontal_n = barrel.horizontal_bottom_n + ball_current_force_n
    chain_top_vertical_n = barrel.vertical_bottom_n - ball_effective_weight(scenario.ball_mass_kg, c)
    chain_state = solve_chain(
        chain_top_horizontal_n, chain_top_vertical_n, scenario.chain, scenario.link_count,
        scenario.current_speed_mps, c,
    )
    model_depth_m = draft_m + rigid_vertical_m + chain_state.vertical_drop_m
    swing_radius_m = rigid_horizontal_m + chain_state.horizontal_span_m + chain_state.seabed_length_m
    return State(
        scenario=scenario,
        draft_m=draft_m,
        pipe_angles_rad=tuple(pipe_angles),
        barrel_angle_rad=barrel.angle_vertical_rad,
        buoy_current_force_n=buoy_current_force_n,
        ball_current_force_n=ball_current_force_n,
        horizontal_tension_at_buoy_n=wind_force_n + buoy_current_force_n,
        chain_top_horizontal_tension_n=chain_top_horizontal_n,
        chain_top_vertical_tension_n=chain_top_vertical_n,
        rigid_horizontal_shift_m=rigid_horizontal_m,
        rigid_vertical_drop_m=rigid_vertical_m,
        chain=chain_state,
        model_depth_m=model_depth_m,
        swing_radius_m=swing_radius_m,
    )


def solve_draft(scenario: Scenario, c: Constants) -> State:
    last_draft: float | None = None
    last_residual: float | None = None
    for index in range(1, 1900):
        draft_m = index / 1000.0
        try:
            residual = evaluate_state(draft_m, scenario, c).model_depth_m - scenario.depth_m
        except ValueError:
            continue
        if last_residual is not None and last_residual * residual <= 0:
            low, high = last_draft, draft_m
            assert low is not None
            for _ in range(90):
                middle = (low + high) / 2.0
                middle_residual = evaluate_state(middle, scenario, c).model_depth_m - scenario.depth_m
                if last_residual * middle_residual <= 0:
                    high = middle
                else:
                    low, last_residual = middle, middle_residual
            return evaluate_state((low + high) / 2.0, scenario, c)
        last_draft, last_residual = draft_m, residual
    raise ValueError("在 0<δ<2 m 的可行域中未找到 H(δ)=目标水深 的根")


def state_row(state: State) -> dict[str, object]:
    s = state.scenario
    return {
        "depth_m": s.depth_m,
        "wind_speed_mps": s.wind_speed_mps,
        "current_speed_mps": s.current_speed_mps,
        "chain_model": s.chain.name,
        "link_count": s.link_count,
        "chain_length_m": s.chain.link_length_m * s.link_count,
        "ball_mass_kg": s.ball_mass_kg,
        "barrel_angle_deg": math.degrees(state.barrel_angle_rad),
        "last_chain_angle_deg": state.chain.anchor_angle_deg,
        "buoy_draft_m": state.draft_m,
        "swing_radius_m": state.swing_radius_m,
        "pipe_1_angle_deg": math.degrees(state.pipe_angles_rad[0]),
        "pipe_2_angle_deg": math.degrees(state.pipe_angles_rad[1]),
        "pipe_3_angle_deg": math.degrees(state.pipe_angles_rad[2]),
        "pipe_4_angle_deg": math.degrees(state.pipe_angles_rad[3]),
        "suspended_chain_length_m": state.chain.suspended_length_m,
        "seabed_chain_length_m": state.chain.seabed_length_m,
        "depth_residual_m": state.model_depth_m - s.depth_m,
        "chain_regime": state.chain.regime,
        "barrel_angle_within_5deg": int(math.degrees(state.barrel_angle_rad) <= 5.0),
        "anchor_angle_within_16deg": int(state.chain.anchor_angle_deg <= 16.0),
    }


def solve_row(scenario: Scenario, c: Constants) -> dict[str, object]:
    try:
        row = state_row(solve_draft(scenario, c))
        row["status"] = "ok"
        row["message"] = ""
        return row
    except ValueError as exc:
        return {
            "depth_m": scenario.depth_m, "wind_speed_mps": scenario.wind_speed_mps,
            "current_speed_mps": scenario.current_speed_mps, "chain_model": scenario.chain.name,
            "link_count": scenario.link_count,
            "chain_length_m": scenario.chain.link_length_m * scenario.link_count,
            "ball_mass_kg": scenario.ball_mass_kg, "status": "infeasible", "message": str(exc),
        }


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    keys: list[str] = []
    for row in rows:
        for key in row:
            if key not in keys:
                keys.append(key)
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=keys, extrasaction="ignore")
    writer.writeheader()
    writer.writerows(rows)
    content = buffer.getvalue()
    write_text_if_unchanged_locked(path, content)


def write_text_if_unchanged_locked(path: Path, content: str) -> None:
    try:
        path.write_text(content, encoding="utf-8", newline="")
    except PermissionError:
        # 某些本地预览程序会短暂锁住 CSV。若计算结果未变，则保留原文件；
        # 若内容不同则显式失败，避免把更新误报为成功。
        if not path.exists() or path.read_text(encoding="utf-8") != content:
            raise


def tex_number(row: dict[str, object], key: str, precision: int = 3) -> str:
    return f"{float(row[key]):.{precision}f}"


def write_tex_table(path: Path, alignment: str, header: list[str], body: list[list[str]]) -> None:
    lines = [
        "% 由 src/q3/q3_current_engineering.py 自动生成，请勿手工修改。",
        rf"\begin{{tabular}}{{{alignment}}}",
        r"\toprule",
        " & ".join(header) + r" \\",
        r"\midrule",
    ]
    lines.extend(" & ".join(row) + r" \\" for row in body)
    lines.extend([r"\bottomrule", r"\end{tabular}", ""])
    write_text_if_unchanged_locked(path, "\n".join(lines))


def write_tex_tables(
    static_rows: list[dict[str, object]], current_chain_rows: list[dict[str, object]],
    link_rows: list[dict[str, object]], wind_rows: list[dict[str, object]],
    current_rows: list[dict[str, object]], ball_rows: list[dict[str, object]],
    depth_rows: list[dict[str, object]],
) -> None:
    metrics_header = [r"$\theta_{\mathrm d}/(^\circ)$", r"$\varphi_{\mathrm a}/(^\circ)$", r"$\delta/\mathrm m$", r"$r/\mathrm m$"]
    metrics = lambda row: [tex_number(row, key) for key in ("barrel_angle_deg", "last_chain_angle_deg", "buoy_draft_m", "swing_radius_m")]
    write_tex_table(
        TABLES / "q3_static_chain_model_relation.tex", "crrrrrrr",
        ["型号", r"$\ell/\mathrm{mm}$", r"$\mu/(\mathrm{kg\,m^{-1}})$", r"$L/\mathrm m$", *metrics_header],
        [[str(row["chain_model"]), str(int(round(float(row["link_length_m"]) * 1000))), tex_number(row, "mass_per_length_kg_m", 2), tex_number(row, "chain_length_m", 2), *metrics(row)] for row in static_rows],
    )
    write_tex_table(
        TABLES / "q3_current_speed_chain_model_relation.tex", "ccrrrr",
        [r"$v_{\mathrm c}/(\mathrm{m\,s^{-1}})$", "型号", *metrics_header],
        [[tex_number(row, "current_speed_mps", 1), str(row["chain_model"]), *metrics(row)] for row in current_chain_rows],
    )
    write_tex_table(
        TABLES / "q3_link_count_effect.tex", "rrrrrr",
        ["链节数", r"$L/\mathrm m$", *metrics_header],
        [[str(int(float(row["link_count"]))), tex_number(row, "chain_length_m", 2), *metrics(row)] for row in link_rows],
    )
    write_tex_table(
        TABLES / "q3_wind_speed_effect.tex", "rrrrr",
        [r"$v_{\mathrm w}/(\mathrm{m\,s^{-1}})$", *metrics_header],
        [[tex_number(row, "wind_speed_mps", 0), *metrics(row)] for row in wind_rows],
    )
    write_tex_table(
        TABLES / "q3_current_speed_effect.tex", "rrrrr",
        [r"$v_{\mathrm c}/(\mathrm{m\,s^{-1}})$", *metrics_header],
        [[tex_number(row, "current_speed_mps", 1), *metrics(row)] for row in current_rows],
    )
    write_tex_table(
        TABLES / "q3_ball_mass_effect.tex", "rrrrr",
        [r"$m_{\mathrm s}/\mathrm{kg}$", *metrics_header],
        [[tex_number(row, "ball_mass_kg", 0), *metrics(row)] for row in ball_rows],
    )
    write_tex_table(
        TABLES / "q3_depth_effect.tex", "rrrrr",
        [r"$H/\mathrm m$", *metrics_header],
        [[tex_number(row, "depth_m", 0), *metrics(row)] for row in depth_rows],
    )


def make_depth_curve(c: Constants) -> tuple[list[dict[str, object]], dict[float, list[tuple[float, float]]]]:
    rows: list[dict[str, object]] = []
    curves: dict[float, list[tuple[float, float]]] = {}
    for current_speed in (0.0, 0.5, 1.0, 1.5):
        scenario = Scenario(18.0, 36.0, current_speed, 1200.0, CHAIN_TYPES[1], 210)
        curve: list[tuple[float, float]] = []
        for index in range(450, 1201):
            draft_m = index / 1000.0
            try:
                state = evaluate_state(draft_m, scenario, c)
            except ValueError:
                continue
            rows.append({
                "wind_speed_mps": 36.0, "current_speed_mps": current_speed,
                "chain_model": "II", "link_count": 210, "ball_mass_kg": 1200.0,
                "draft_m": draft_m, "model_depth_m": state.model_depth_m,
            })
            curve.append((draft_m, state.model_depth_m))
        curves[current_speed] = curve
    return rows, curves


def plot_depth_curve(curves: dict[float, list[tuple[float, float]]]) -> None:
    styles = {
        0.0: ("#222222", "-", "o"),
        0.5: ("#1f4e79", "--", "s"),
        1.0: ("#9a5a12", "-.", "^"),
        1.5: ("#8a2525", ":", "D"),
    }
    figure, axis = plt.subplots(figsize=(6.85, 4.25))
    for current_speed, values in curves.items():
        x_values, y_values = zip(*values)
        color, line_style, marker = styles[current_speed]
        axis.plot(
            x_values, y_values, color=color, linestyle=line_style, marker=marker,
            markevery=50, markersize=3.5, markerfacecolor="white", markeredgewidth=0.75,
            linewidth=1.8, label=f"v_c = {current_speed:g} m/s",
        )
    axis.set_xlim(0.68, 0.81)
    axis.set_ylim(15.5, 20.5)
    axis.xaxis.set_major_locator(MultipleLocator(0.02))
    axis.xaxis.set_minor_locator(MultipleLocator(0.01))
    axis.yaxis.set_major_locator(MultipleLocator(1.0))
    axis.yaxis.set_minor_locator(MultipleLocator(0.5))
    axis.grid(which="major", color="#d0d0d0", linewidth=0.45)
    axis.tick_params(which="major", length=3.5, labelsize=8)
    axis.tick_params(which="minor", length=2)
    axis.set_xlabel("浮标吃水 $\\delta$ / m", fontproperties=CHINESE_FONT, fontsize=9)
    axis.set_ylabel("海水深度 $H$ / m", fontproperties=CHINESE_FONT, fontsize=9)
    legend = axis.legend(loc="upper left", frameon=True, facecolor="white", edgecolor="#555555", fontsize=8)
    for text in legend.get_texts():
        text.set_fontproperties(CHINESE_FONT)
    figure.subplots_adjust(left=0.13, right=0.98, bottom=0.16, top=0.98)
    for suffix in ("pdf", "svg", "png"):
        figure.savefig(FIGURES / f"q3_depth_draft_curve.{suffix}", dpi=600, bbox_inches="tight")
    plt.close(figure)


def anchored_chain_coordinates(state: State) -> list[tuple[float, float]]:
    chain = state.chain
    return [
        (chain.seabed_length_m + chain.horizontal_span_m - x_from_top, -state.draft_m - state.rigid_vertical_drop_m - z_drop)
        for x_from_top, z_drop in chain.coordinates_from_top
    ]


def write_chain_coordinates(state: State, filename: str) -> None:
    coordinates = anchored_chain_coordinates(state)
    rows = [
        {"node_index": index, "x_m": x_value, "z_m": z_value}
        for index, (x_value, z_value) in enumerate(coordinates)
    ]
    write_csv(PROCESSED / filename, rows)


def plot_chain_shapes(states: dict[float, State]) -> None:
    figure, axes = plt.subplots(1, 2, figsize=(6.85, 3.55))
    colors = {16.0: "#1f4e79", 20.0: "#8a2525"}
    panels = {16.0: "(a)", 20.0: "(b)"}
    for axis, depth_m in zip(axes, sorted(states)):
        state = states[depth_m]
        x_values, z_values = zip(*anchored_chain_coordinates(state))
        color = colors[depth_m]
        axis.plot(x_values, z_values, color=color, linewidth=1.8, marker="o", markevery=30,
                  markersize=2.8, markerfacecolor="white", markeredgewidth=0.65)
        axis.axhline(0.0, color="#266a96", linewidth=0.8)
        axis.axhline(-depth_m, color="#6d5035", linestyle="--", linewidth=0.8)
        axis.plot(0.0, -depth_m, marker="v", color="#222222", markersize=4.8)
        axis.plot(x_values[0], z_values[0], marker="s", color="#222222", markersize=4.2)
        axis.text(0.04, 0.94, f"{panels[depth_m]} $H={depth_m:g}$ m", transform=axis.transAxes,
                  va="top", fontsize=8.5)
        axis.text(0.04, 0.84, "海平面", transform=axis.transAxes, color="#266a96",
                  fontproperties=CHINESE_FONT, fontsize=8)
        axis.text(0.04, 0.07, "海床", transform=axis.transAxes, color="#6d5035",
                  fontproperties=CHINESE_FONT, fontsize=8)
        axis.annotate("锚点", xy=(0.0, -depth_m), xytext=(5, 7), textcoords="offset points",
                      fontproperties=CHINESE_FONT, fontsize=7.5)
        axis.annotate("钢桶--锚链节点", xy=(x_values[0], z_values[0]), xytext=(-2, 7),
                      textcoords="offset points", ha="right", fontproperties=CHINESE_FONT, fontsize=7.2)
        axis.set_xlim(-0.5, max(x_values) + 1.0)
        axis.set_ylim(-depth_m - 1.0, 1.0)
        axis.xaxis.set_major_locator(MultipleLocator(5.0))
        axis.xaxis.set_minor_locator(MultipleLocator(2.5))
        axis.yaxis.set_major_locator(MultipleLocator(4.0))
        axis.yaxis.set_minor_locator(MultipleLocator(2.0))
        axis.grid(which="major", color="#d5d5d5", linewidth=0.4)
        axis.tick_params(which="major", length=3.2, labelsize=7.5)
        axis.tick_params(which="minor", length=1.8)
        axis.set_xlabel("水平坐标 $x$ / m", fontproperties=CHINESE_FONT, fontsize=8.5)
    axes[0].set_ylabel("高程 $z$ / m", fontproperties=CHINESE_FONT, fontsize=8.5)
    figure.subplots_adjust(left=0.095, right=0.99, bottom=0.18, top=0.97, wspace=0.28)
    for suffix in ("pdf", "svg", "png"):
        figure.savefig(FIGURES / f"q3_chain_shapes_depth_16_20.{suffix}", dpi=600, bbox_inches="tight")
    plt.close(figure)


def main() -> None:
    for directory in (PROCESSED, TABLES, FIGURES, LOGS):
        directory.mkdir(parents=True, exist_ok=True)
    c = Constants()
    base = dict(depth_m=18.0, wind_speed_mps=36.0, current_speed_mps=1.5, ball_mass_kg=1200.0, chain=CHAIN_TYPES[1], link_count=210)

    static_rows: list[dict[str, object]] = []
    for chain in CHAIN_TYPES:
        row = solve_row(Scenario(18.0, 36.0, 0.0, 1200.0, chain, 210), c)
        row.update({"link_length_m": chain.link_length_m, "mass_per_length_kg_m": chain.mass_per_length_kg_m, "chain_mass_kg": chain.link_length_m * 210 * chain.mass_per_length_kg_m})
        static_rows.append(row)
    write_csv(TABLES / "q3_static_chain_model_relation.csv", static_rows)

    current_chain_rows = [
        solve_row(Scenario(18.0, 36.0, current_speed, 1200.0, chain, 210), c)
        for current_speed in (0.0, 0.5, 1.0, 1.5)
        for chain in CHAIN_TYPES
    ]
    write_csv(TABLES / "q3_current_speed_chain_model_relation.csv", current_chain_rows)

    link_rows = [
        solve_row(Scenario(18.0, 36.0, 1.5, 1200.0, CHAIN_TYPES[1], link_count), c)
        for link_count in (180, 200, 210, 220, 240)
    ]
    write_csv(TABLES / "q3_link_count_effect.csv", link_rows)

    wind_rows = [
        solve_row(Scenario(18.0, wind_speed, 1.5, 1200.0, CHAIN_TYPES[1], 210), c)
        for wind_speed in (12.0, 24.0, 36.0)
    ]
    write_csv(TABLES / "q3_wind_speed_effect.csv", wind_rows)

    current_rows = [
        solve_row(Scenario(18.0, 36.0, current_speed, 1200.0, CHAIN_TYPES[1], 210), c)
        for current_speed in (0.0, 0.5, 1.0, 1.5)
    ]
    write_csv(TABLES / "q3_current_speed_effect.csv", current_rows)

    ball_rows = [
        solve_row(Scenario(18.0, 36.0, 1.5, ball_mass, CHAIN_TYPES[1], 210), c)
        for ball_mass in (800.0, 1000.0, 1200.0, 1500.0, 1800.0, 2100.0, 2400.0, 2700.0, 3000.0)
    ]
    write_csv(TABLES / "q3_ball_mass_effect.csv", ball_rows)

    depth_rows = [
        solve_row(Scenario(depth_m, 36.0, 1.5, 1200.0, CHAIN_TYPES[1], 210), c)
        for depth_m in (16.0, 17.0, 18.0, 19.0, 20.0)
    ]
    write_csv(TABLES / "q3_depth_effect.csv", depth_rows)
    write_tex_tables(static_rows, current_chain_rows, link_rows, wind_rows, current_rows, ball_rows, depth_rows)

    curve_rows, curves = make_depth_curve(c)
    write_csv(PROCESSED / "q3_depth_draft_curve.csv", curve_rows)
    plot_depth_curve(curves)

    shape_states: dict[float, State] = {}
    for depth_m in (16.0, 20.0):
        scenario = Scenario(depth_m, 36.0, 1.5, 1200.0, CHAIN_TYPES[1], 210)
        state = solve_draft(scenario, c)
        shape_states[depth_m] = state
        write_chain_coordinates(state, f"q3_chain_shape_depth_{int(depth_m)}m.csv")
    plot_chain_shapes(shape_states)

    checks = [
        "Generated by src/q3/q3_current_engineering.py",
        "Engineering baseline for one-factor sweeps: H=18 m, v_w=36 m/s, v_c=1.5 m/s, chain II, 210 links, m_s=1200 kg.",
        "Static chain-model table changes only chain type; current-chain table changes chain type and v_c.",
        "Water-current force directions are aligned with wind; the ball uses S=pi*D_s^2/4.",
        "After touchdown, spare chain is laid on the flat seabed and receives no additional current drag in this engineering discretization.",
        "No paper section or model card is written by this script.",
    ]
    all_rows = static_rows + current_chain_rows + link_rows + wind_rows + current_rows + ball_rows + depth_rows
    checked = 0
    for row in all_rows:
        if row.get("status") == "ok":
            checked += 1
            if abs(float(row["depth_residual_m"])) > 1e-8:
                raise AssertionError("水深闭合残差超限")
    checks.append(f"successful_scenarios={checked}")
    (LOGS / "q3_engineering_checks.txt").write_text("\n".join(checks) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
