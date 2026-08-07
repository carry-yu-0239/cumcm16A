# 问题三工程计算说明（未写入论文）

## 范围与追溯

本工程实现位于 `src/q3/q3_current_engineering.py`，对应链条为：题面与既有假设 → 建模手给出的风力、水流力、构件递推与力矩式 → 数值求解脚本 → `data/processed/q3_*` 和 `outputs/q3/`。本次未修改 `paper/sections/08_problem3.tex`、模型卡或结果总账。

本次使用了带 Matplotlib 的工作区 Python 运行时；图形同时写出 PDF、PNG 与 UTF-8 SVG，正文优先引用矢量 PDF。

已沿用项目已有口径：海水密度 (1025\ \mathrm{kg/m^3})、(g=9.8\ \mathrm{m/s^2})、钢材密度 (7850\ \mathrm{kg/m^3})、重物球为实心钢球、锚链链材直径 (17.5\ \mathrm{mm})。钢管和钢桶的夹角从竖直方向量起，锚链末端夹角从水平海床量起。

## 递推实现

浮标的水平载荷取 (F_W+F_{cb})，其中 (F_{cb}=374(2\delta)v_c^2)。对长度 (L_i)、直径 (D_i) 的钢管或钢桶，按手稿采用

\[
F_{ci}=374D_iL_i\cos\varphi_i\,v_c^2,
\quad H_{i+1}=H_i+F_{ci},
\quad V_{i+1}=V_i-(G_i-F_i).
\]

构件中点的力矩平衡给出隐式角度方程

\[
\tan\varphi_i=\frac{H_i+F_{ci}/2}{V_i-(G_i-F_i)/2},
\]

每个构件以二分法求解。重物球的水流投影面积取 (S_s=\pi D_s^2/4)，其水平水流力作为钢桶—锚链节点的集中水平荷载；这与手稿中的 (F_{cs}) 符号对应。

锚链按链节离散。悬空链节相对水平的夹角由同样的中点平衡递推，链节水流投影采用题面和手稿的 (D_cL_c\sin\alpha) 形式。若竖向张力降到零，余链在平坦海床上平铺；由于题面没有海床摩擦和卧底链流阻参数，平铺段不再施加水流阻力。这是数值离散约定，不是新的正式模型结论。

## 输出与扫描口径

`outputs/q3/tables/` 中有七张 UTF-8 CSV 表，覆盖：静水时链型、流速—链型、链节数、风速（12/24/36 m/s）、流速、重物球质量与水深的影响。除正在改变的自变量外，单因素扫描采用工程基准：

\[
H=18\ \mathrm{m},\quad v_w=36\ \mathrm{m/s},\quad v_c=1.5\ \mathrm{m/s},
\quad \text{II 型锚链、210 节、}m_s=1200\ \mathrm{kg}.
\]

这组基准仅用于让表格可复现，并不等于问题三的最终设计选择。结果中的 `infeasible` 表示在 (0<\delta<2\ \mathrm{m}) 内未找到满足指定水深的静力闭合根，不能被解读为工程上绝对不可布放。

| 需求表格 | 输出文件 | 变化的自变量 |
|---|---|---|
| 海水静止时锚链型号与参数的关系 | `q3_static_chain_model_relation.csv` | 链型 I--V |
| 海水流速与锚链型号的关系 | `q3_current_speed_chain_model_relation.csv` | (v_c=0,0.5,1.0,1.5\ \mathrm{m/s}) 与链型 |
| 不同链条节数与各参数的关系 | `q3_link_count_effect.csv` | II 型链节数 180、200、210、220、240 |
| 不同风速对各指标的影响 | `q3_wind_speed_effect.csv` | (v_w=12,24,36\ \mathrm{m/s}) |
| 水流速度对各参数的影响 | `q3_current_speed_effect.csv` | (v_c=0,0.5,1.0,1.5\ \mathrm{m/s}) |
| 重物球质量与各参数的关系 | `q3_ball_mass_effect.csv` | (m_s=800\) 至 (3000\ \mathrm{kg}) |
| 海水深度与各参数的关系 | `q3_depth_effect.csv` | (H=16,17,18,19,20\ \mathrm{m}) |

每张结果表均含钢桶夹角、最后一条锚链角、浮标吃水、游动半径，以及四节钢管角度、悬空/卧底链长、闭合残差与两个角度约束标记。

`data/processed/q3_depth_draft_curve.csv` 与 `outputs/q3/figures/q3_depth_draft_curve.{png,svg}` 给出以吃水为自变量、以海水深度为因变量的曲线。`q3_chain_shapes_depth_16_20.{png,svg}` 及对应 CSV 给出 (H=16\) m 和 (H=20\) m 的锚链形状，使用相同的工程基准风流与构型。

## 验证层级

- 已完成：脚本执行、每个成功场景的水深闭合残差检查、锚链总长度守恒检查、CSV UTF-8 写出。
- 已完成：PNG/SVG 图形导出；字体优先为宋体等中文衬线字体，数学符号按 LaTex 数学字体显示。
- 未完成：不同深度实测流速分布、海床摩擦、链环精细投影系数与浮标姿态的物理验证；也未做 MATLAB 复现、论文编译或人工图形验收。
