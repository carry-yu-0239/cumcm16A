%% CUMCM 2016 A题 - 问题二完整程序
% 说明：本程序严格沿用问题一的静力递推思路。
% 与问题一相比只改变两点：
%   1. 海面风速取 36 m/s；
%   2. 将重物球质量 m_s 作为外层搜索变量。
%
% 问题一模型口径：
%   - 浮标平浮，先由浮标水平、竖直受力平衡确定顶部张力分量；
%   - 四节钢管和钢桶均按刚性圆柱逐节进行受力、力矩平衡；
%   - II 型锚链按链环长度 0.105 m 离散为 210 节，逐链环递推；
%   - 以浮标吃水深度 delta 为内层未知量，用二分法满足 H(delta)=18 m；
%   - 外层搜索重物球质量，使 theta_d <= 5 deg 且 phi_a <= 16 deg。
%
% 兼容：MATLAB R2025b；不需要额外工具箱。
% 运行：直接在 MATLAB 中运行本文件即可。
% 输出：处理数据写入 data/processed，结果写入 outputs/q2。

clear; clc; close all;

%% 0. 输出目录：全部写入仓库约定目录
script_dir = fileparts(mfilename('fullpath'));
root_dir = fullfile(script_dir, '..', '..');
processed_dir = fullfile(root_dir, 'data', 'processed');
tables_dir = fullfile(root_dir, 'outputs', 'q2', 'tables');
figures_dir = fullfile(root_dir, 'outputs', 'q2', 'figures');
logs_dir = fullfile(root_dir, 'outputs', 'q2', 'logs');
output_dirs = {processed_dir, tables_dir, figures_dir, logs_dir};
for k = 1:numel(output_dirs)
    if ~exist(output_dirs{k}, 'dir')
        mkdir(output_dirs{k});
    end
end

%% 1. 参数：与问题一保持一致
p.g = 9.8;                        % 重力加速度 m/s^2
p.rho = 1025;                     % 海水密度 kg/m^3
p.depth = 18;                     % 水深 H, m

% 浮标
p.buoy_radius = 1.0;              % m
p.buoy_height = 2.0;              % m
p.buoy_mass = 1000;               % kg

% 四节钢管：完全密闭，按外轮廓排水体积计算浮力
p.pipe_length = 1.0;              % m
p.pipe_radius = 0.025;            % m
p.pipe_mass = 10;                 % kg
p.n_pipes = 4;

% 密封钢桶
p.barrel_length = 1.0;            % m
p.barrel_radius = 0.15;           % m
p.barrel_mass = 100;              % kg

% 钢材、重物球
p.rho_steel = 7850;               % kg/m^3

% II 型锚链
p.chain_length = 22.05;           % m
p.chain_link_length = 0.105;      % m
p.chain_mass_per_length = 7.0;    % kg/m
p.n_chain_links = round(p.chain_length / p.chain_link_length);
if abs(p.n_chain_links * p.chain_link_length - p.chain_length) > 1e-10
    error('锚链总长不能由给定链环长度整除，请检查参数。');
end

% 水中有效重力
p.pipe_effective_weight = p.pipe_mass*p.g - ...
    p.rho*p.g*pi*p.pipe_radius^2*p.pipe_length;
p.barrel_effective_weight = p.barrel_mass*p.g - ...
    p.rho*p.g*pi*p.barrel_radius^2*p.barrel_length;
p.chain_effective_weight_per_length = p.chain_mass_per_length*p.g * ...
    (1 - p.rho/p.rho_steel);

wind_speed = 36;                  % 问题二海面风速 m/s
original_mass = 1200;             % 问题一原配重 kg
selected_mass = 2250;             % 最终推荐配重 kg（留出小幅约束裕量）

%% 2. 原 1200 kg 方案
state_original = solve_draft(wind_speed, original_mass, p);

%% 3. 按 10 kg 步长搜索最小可行质量
mass_grid = (1200:10:3000)';
nm = numel(mass_grid);
barrel_angle_deg = nan(nm,1);
anchor_angle_deg = nan(nm,1);
draft_m = nan(nm,1);
radius_m = nan(nm,1);
feasible = false(nm,1);

for k = 1:nm
    s = solve_draft(wind_speed, mass_grid(k), p);
    barrel_angle_deg(k) = rad2deg(s.barrel_angle);
    anchor_angle_deg(k) = s.chain.anchor_angle_deg;
    draft_m(k) = s.draft;
    radius_m(k) = s.radius;
    feasible(k) = barrel_angle_deg(k) <= 5 && anchor_angle_deg(k) <= 16;
end

idx = find(feasible, 1, 'first');
if isempty(idx)
    error('在 1200-3000 kg 搜索区间内未找到满足两个角度约束的方案。');
end
min_grid_mass = mass_grid(idx);
state_grid = solve_draft(wind_speed, min_grid_mass, p);

mass_scan = table(mass_grid, barrel_angle_deg, anchor_angle_deg, ...
    draft_m, radius_m, feasible, ...
    'VariableNames', {'ball_mass_kg','barrel_angle_deg','anchor_angle_deg', ...
    'buoy_draft_m','swing_radius_m','feasible'});
writetable(mass_scan, fullfile(processed_dir, 'q2_mass_scan.csv'));

%% 4. 二分搜索理论临界质量
critical_mass = solve_critical_mass(wind_speed, p, 1200, 3000);
state_critical = solve_draft(wind_speed, critical_mass, p);

%% 5. 最终 2250 kg 方案
state_selected = solve_draft(wind_speed, selected_mass, p);

%% 6. 输出汇总表
scenario = ["原1200kg方案"; "理论临界方案"; "10kg步长最小可行方案"; "最终2250kg方案"];
ball_mass_kg = [original_mass; critical_mass; min_grid_mass; selected_mass];
states = {state_original; state_critical; state_grid; state_selected};
result_table = make_result_table(scenario, ball_mass_kg, states);
writetable(result_table, fullfile(tables_dir, 'q2_results.csv'));
write_result_table_tex(fullfile(tables_dir, 'q2_results.tex'), result_table);

%% 7. 输出锚链坐标
chain_1200 = chain_coordinates_from_anchor(state_original, p);
chain_2250 = chain_coordinates_from_anchor(state_selected, p);
writetable(chain_1200, fullfile(processed_dir, 'q2_chain_1200kg.csv'));
writetable(chain_2250, fullfile(processed_dir, 'q2_chain_2250kg.csv'));

%% 8. 绘图
plot_chain_shapes(chain_1200, chain_2250, p, figures_dir);
plot_mass_constraints(mass_scan, min_grid_mass, selected_mass, figures_dir);

%% 9. 输出文本摘要与自检日志
summary_path = fullfile(logs_dir, 'q2_summary.txt');
fid = fopen(summary_path, 'w', 'n', 'UTF-8');
if fid < 0
    error('无法创建输出摘要文件。');
end
write_state_text(fid, '原 1200 kg 方案', state_original, original_mass);
fprintf(fid, '\n理论临界质量 = %.6f kg\n', critical_mass);
fprintf(fid, '10 kg 步长最小可行质量 = %.0f kg\n', min_grid_mass);
write_state_text(fid, '最终 2250 kg 方案', state_selected, selected_mass);
fclose(fid);

%% 10. 命令行结果
fprintf('\n============================================\n');
fprintf('问题二计算完成（严格沿用问题一递推模型）\n');
fprintf('============================================\n');
print_state('原 1200 kg 方案', state_original, original_mass);
fprintf('理论临界质量：%.6f kg\n', critical_mass);
fprintf('10 kg 步长最小可行质量：%.0f kg\n', min_grid_mass);
print_state('最终 2250 kg 方案', state_selected, selected_mass);
fprintf('输出目录：%s\n', root_dir);

%% ======================== Local functions ========================
function state = solve_draft(wind_speed, ball_mass, p)
% 对给定风速和重物球质量，以 delta 为未知量求 H(delta)=18 m。

    draft_grid = linspace(0.30, 1.60, 800);
    prev_valid = false;
    low = NaN; high = NaN;

    for k = 1:numel(draft_grid)
        delta = draft_grid(k);
        try
            s = evaluate_state(delta, wind_speed, ball_mass, p);
            f = s.model_depth - p.depth;
        catch
            continue;
        end

        if prev_valid && prev_f*f <= 0
            low = prev_delta;
            high = delta;
            break;
        end
        prev_valid = true;
        prev_delta = delta;
        prev_f = f;
    end

    if isnan(low)
        error('未找到 H(delta)=18 m 的吃水根区间。');
    end

    f_low = evaluate_state(low, wind_speed, ball_mass, p).model_depth - p.depth;
    for k = 1:100
        mid = 0.5*(low + high);
        s_mid = evaluate_state(mid, wind_speed, ball_mass, p);
        f_mid = s_mid.model_depth - p.depth;
        if f_low*f_mid <= 0
            high = mid;
        else
            low = mid;
            f_low = f_mid;
        end
    end

    state = evaluate_state(0.5*(low + high), wind_speed, ball_mass, p);
end

function state = evaluate_state(delta, wind_speed, ball_mass, p)
% 沿用问题一：浮标 -> 4节钢管 -> 钢桶/重物球 -> 210节锚链。

    if delta <= 0 || delta >= p.buoy_height
        error('浮标吃水 delta 必须位于 (0,2) m。');
    end

    % (1) 浮标受力平衡
    Fb = p.rho*p.g*pi*p.buoy_radius^2*delta;
    Gb = p.buoy_mass*p.g;
    Fw = 0.625 * ((p.buoy_height-delta)*2*p.buoy_radius) * wind_speed^2;

    Th = Fw;                       % 系泊张力水平分量
    Tv = Fb - Gb;                  % 浮标下接构件顶部的竖向张力分量
    if Th <= 0 || Tv <= 0
        error('当前吃水下顶部张力分量不满足受拉平衡条件。');
    end

    % (2) 四节钢管：角度均按“与竖直方向夹角”定义
    pipe_angles = zeros(p.n_pipes,1);
    rigid_vertical = 0;
    rigid_horizontal = 0;
    current_Tv = Tv;

    for i = 1:p.n_pipes
        pipe_angles(i) = rigid_member_angle_from_vertical( ...
            Th, current_Tv, p.pipe_effective_weight);
        rigid_vertical = rigid_vertical + p.pipe_length*cos(pipe_angles(i));
        rigid_horizontal = rigid_horizontal + p.pipe_length*sin(pipe_angles(i));
        current_Tv = current_Tv - p.pipe_effective_weight;
    end

    % (3) 钢桶：仍按问题一的刚体受力、力矩平衡
    barrel_angle = rigid_member_angle_from_vertical( ...
        Th, current_Tv, p.barrel_effective_weight);
    rigid_vertical = rigid_vertical + p.barrel_length*cos(barrel_angle);
    rigid_horizontal = rigid_horizontal + p.barrel_length*sin(barrel_angle);
    current_Tv = current_Tv - p.barrel_effective_weight;

    % (4) 重物球作用在钢桶下端连接点，进入锚链顶部竖向张力
    ball_effective_weight = ball_mass*p.g*(1 - p.rho/p.rho_steel);
    chain_top_Tv = current_Tv - ball_effective_weight;
    if chain_top_Tv <= 0
        error('重物球过重或吃水过浅，锚链顶部竖向张力非正。');
    end

    % (5) II 型锚链逐链环递推；每节 0.105 m，共 210 节
    chain = solve_chain_links(Th, chain_top_Tv, p);

    state = struct();
    state.wind_speed = wind_speed;
    state.ball_mass = ball_mass;
    state.draft = delta;
    state.wind_force = Fw;
    state.horizontal_tension = Th;
    state.top_vertical_tension = Tv;
    state.pipe_angles = pipe_angles;
    state.barrel_angle = barrel_angle;
    state.rigid_vertical_drop = rigid_vertical;
    state.rigid_horizontal_shift = rigid_horizontal;
    state.chain_top_vertical_tension = chain_top_Tv;
    state.chain = chain;
    state.model_depth = delta + rigid_vertical + chain.vertical_drop;
    state.radius = rigid_horizontal + chain.horizontal_span;
end

function theta = rigid_member_angle_from_vertical(Th, Tv_top, W)
% 匀质刚性圆柱：由受力平衡和绕端点力矩平衡得到。
% theta 为构件轴线与竖直方向的夹角。
    Tv_mid = Tv_top - 0.5*W;
    if Tv_mid <= 0
        error('刚性构件中点等效竖向张力非正，当前受拉平衡失效。');
    end
    theta = atan2(Th, Tv_mid);
end

function chain = solve_chain_links(Th, Tv_top, p)
% 锚链按问题一的“逐节受力 + 力矩平衡”思想离散计算。
% 内部 alpha_j 定义为第 j 节链环与水平海床的夹角。

    n = p.n_chain_links;
    ell = p.chain_link_length;
    w = p.chain_effective_weight_per_length;

    alpha = zeros(n,1);
    dx = zeros(n,1);
    dz = zeros(n,1);
    x_from_top = zeros(n+1,1);
    z_from_top = zeros(n+1,1);

    Tv = Tv_top;
    suspended_length = 0;
    seabed_length = 0;

    for j = 1:n
        if Tv <= 1e-12
            % 已到达海床，余下链环全部卧底
            alpha(j) = 0;
            dx(j) = ell;
            dz(j) = 0;
            seabed_length = seabed_length + ell;
            Tv = 0;
        else
            Tv_bottom = Tv - w*ell;
            if Tv_bottom >= 0
                % 整节悬空：链环姿态由两端竖向张力平均值确定
                Tv_mid = 0.5*(Tv + Tv_bottom);
                alpha(j) = atan2(Tv_mid, Th);   % 与水平面夹角
                dx(j) = ell*cos(alpha(j));
                dz(j) = ell*sin(alpha(j));
                suspended_length = suspended_length + ell;
                Tv = Tv_bottom;
            else
                % 触地点位于本链环内部：悬空部分 + 卧底部分
                ell_s = Tv / w;
                ell_b = ell - ell_s;
                Tv_mid = 0.5*Tv;
                alpha(j) = atan2(Tv_mid, Th);
                dx(j) = ell_s*cos(alpha(j)) + ell_b;
                dz(j) = ell_s*sin(alpha(j));
                suspended_length = suspended_length + ell_s;
                seabed_length = seabed_length + ell_b;
                Tv = 0;
            end
        end

        x_from_top(j+1) = x_from_top(j) + dx(j);
        z_from_top(j+1) = z_from_top(j) + dz(j);
    end

    if seabed_length > 1e-10
        anchor_angle_deg = 0;
        regime = "touchdown_with_seabed_contact";
    else
        anchor_angle_deg = atan2d(Tv, Th);
        regime = "fully_suspended";
    end

    chain = struct();
    chain.link_angle_from_horizontal = alpha;
    chain.dx = dx;
    chain.dz = dz;
    chain.x_from_top = x_from_top;
    chain.z_from_top = z_from_top;
    chain.horizontal_span = x_from_top(end);
    chain.vertical_drop = z_from_top(end);
    chain.suspended_length = suspended_length;
    chain.seabed_length = seabed_length;
    chain.bottom_vertical_tension = Tv;
    chain.anchor_angle_deg = anchor_angle_deg;
    chain.regime = regime;
end

function tf = is_feasible_mass(ball_mass, wind_speed, p)
    s = solve_draft(wind_speed, ball_mass, p);
    tf = rad2deg(s.barrel_angle) <= 5 && s.chain.anchor_angle_deg <= 16;
end

function critical_mass = solve_critical_mass(wind_speed, p, low_mass, high_mass)
    if is_feasible_mass(low_mass, wind_speed, p)
        error('给定质量下界已经可行，不能用于临界质量二分。');
    end
    if ~is_feasible_mass(high_mass, wind_speed, p)
        error('给定质量上界仍不可行，请增大质量上界。');
    end

    for k = 1:80
        mid_mass = 0.5*(low_mass + high_mass);
        if is_feasible_mass(mid_mass, wind_speed, p)
            high_mass = mid_mass;
        else
            low_mass = mid_mass;
        end
    end
    critical_mass = 0.5*(low_mass + high_mass);
end

function T = make_result_table(scenario, ball_mass_kg, states)
    n = numel(states);
    theta1 = zeros(n,1); theta2 = zeros(n,1);
    theta3 = zeros(n,1); theta4 = zeros(n,1);
    theta_d = zeros(n,1); phi_a = zeros(n,1);
    delta = zeros(n,1); radius = zeros(n,1);
    suspended = zeros(n,1); seabed = zeros(n,1);
    depth_residual = zeros(n,1); regime = strings(n,1);

    for k = 1:n
        s = states{k};
        a = rad2deg(s.pipe_angles);
        theta1(k)=a(1); theta2(k)=a(2); theta3(k)=a(3); theta4(k)=a(4);
        theta_d(k)=rad2deg(s.barrel_angle);
        phi_a(k)=s.chain.anchor_angle_deg;
        delta(k)=s.draft;
        radius(k)=s.radius;
        suspended(k)=s.chain.suspended_length;
        seabed(k)=s.chain.seabed_length;
        depth_residual(k)=s.model_depth-18;
        regime(k)=s.chain.regime;
    end

    T = table(scenario, ball_mass_kg, theta1, theta2, theta3, theta4, ...
        theta_d, phi_a, delta, radius, suspended, seabed, depth_residual, regime, ...
        'VariableNames', {'scenario','ball_mass_kg','theta1_deg','theta2_deg', ...
        'theta3_deg','theta4_deg','theta_d_deg','phi_a_deg','delta_m','radius_m', ...
        'suspended_chain_m','seabed_chain_m','depth_residual_m','chain_regime'});
end

function write_result_table_tex(path, T)
% 与 q2_results.csv 同步输出的 LaTeX 汇总表，供论文复核使用。
    fid = fopen(path, 'w', 'n', 'UTF-8');
    if fid < 0
        error('无法创建问题二 LaTeX 结果表。');
    end
    fprintf(fid, '%% 由 src/q2/q2_complete.m 自动生成。\n');
    fprintf(fid, '\\begin{tabular}{lrrrr}\n');
    fprintf(fid, '\\toprule\n');
    fprintf(fid, '方案 & 重物球质量/kg & 钢桶倾角/(度) & 锚端夹角/(度) & 游动半径/m \\\\ \n');
    fprintf(fid, '\\midrule\n');
    for k = 1:height(T)
        fprintf(fid, '%s & %.6f & %.4f & %.4f & %.4f \\\\ \n', ...
            char(T.scenario(k)), T.ball_mass_kg(k), T.theta_d_deg(k), ...
            T.phi_a_deg(k), T.radius_m(k));
    end
    fprintf(fid, '\\bottomrule\n');
    fprintf(fid, '\\end{tabular}\n');
    fclose(fid);
end

function C = chain_coordinates_from_anchor(state, p)
% 将“从锚链顶端向下”的坐标转换为“锚点为原点、海床为 z=0”的坐标。
    xt = state.chain.x_from_top;
    zt = state.chain.z_from_top;
    x = state.chain.horizontal_span - xt;
    z = state.chain.vertical_drop - zt;
    link_index = (0:p.n_chain_links)';
    C = table(link_index, x, z, 'VariableNames', {'node_index','x_m','z_m'});
end

function plot_chain_shapes(C1, C2, p, out_dir)
    fig = figure('Visible','off','Color','w','Position',[100 100 1050 720]);
    hold on; grid on; box on;
    h1 = plot(C1.x_m, C1.z_m, 'LineWidth', 2.0);
    h2 = plot(C2.x_m, C2.z_m, 'LineWidth', 2.0);
    h3 = yline(0, '-', 'LineWidth', 1.0);
    h4 = yline(p.depth, '--', 'LineWidth', 1.0);
    plot(0,0,'ko','MarkerFaceColor','k','MarkerSize',6);
    xlabel('相对锚点的水平距离 / m');
    ylabel('距海床高度 / m');
    title('风速 36 m/s 时调节配重前后的锚链形状');
    legend([h1 h2 h3 h4], {'1200 kg','2250 kg','海床','海面'}, 'Location','best');
    exportgraphics(fig, fullfile(out_dir,'q2_chain_shapes.png'), 'Resolution',300);
    close(fig);
end

function plot_mass_constraints(T, min_grid_mass, selected_mass, out_dir)
    fig = figure('Visible','off','Color','w','Position',[100 100 1050 720]);
    hold on; grid on; box on;
    h1 = plot(T.ball_mass_kg, T.barrel_angle_deg, 'LineWidth',2.0);
    h2 = plot(T.ball_mass_kg, T.anchor_angle_deg, 'LineWidth',2.0);
    h3 = yline(5,'--','5 deg');
    h4 = yline(16,'--','16 deg');
    xline(min_grid_mass, ':', sprintf('%.0f kg',min_grid_mass));
    xline(selected_mass, ':', sprintf('%.0f kg',selected_mass));
    xlabel('重物球质量 m_s / kg');
    ylabel('角度 / deg');
    title('重物球质量对钢桶倾角和锚端夹角的影响');
    legend([h1 h2 h3 h4], {'钢桶倾角','锚端夹角','钢桶约束','锚端约束'}, 'Location','northeast');
    exportgraphics(fig, fullfile(out_dir,'q2_mass_constraints.png'), 'Resolution',300);
    close(fig);
end

function print_state(name, s, mass)
    a = rad2deg(s.pipe_angles);
    fprintf('\n%s\n', name);
    fprintf('m_s = %.6f kg\n', mass);
    fprintf('theta_1 = %.4f deg\n', a(1));
    fprintf('theta_2 = %.4f deg\n', a(2));
    fprintf('theta_3 = %.4f deg\n', a(3));
    fprintf('theta_4 = %.4f deg\n', a(4));
    fprintf('theta_d = %.4f deg\n', rad2deg(s.barrel_angle));
    fprintf('phi_a   = %.4f deg\n', s.chain.anchor_angle_deg);
    fprintf('delta   = %.4f m\n', s.draft);
    fprintf('r       = %.4f m\n', s.radius);
    fprintf('chain suspended = %.4f m, seabed = %.4f m\n', ...
        s.chain.suspended_length, s.chain.seabed_length);
    fprintf('|H_model-18| = %.3e m\n', abs(s.model_depth-18));
end

function write_state_text(fid, name, s, mass)
    a = rad2deg(s.pipe_angles);
    fprintf(fid, '\n[%s]\n', name);
    fprintf(fid, 'm_s = %.6f kg\n', mass);
    fprintf(fid, 'theta_1 = %.6f deg\n', a(1));
    fprintf(fid, 'theta_2 = %.6f deg\n', a(2));
    fprintf(fid, 'theta_3 = %.6f deg\n', a(3));
    fprintf(fid, 'theta_4 = %.6f deg\n', a(4));
    fprintf(fid, 'theta_d = %.6f deg\n', rad2deg(s.barrel_angle));
    fprintf(fid, 'phi_a = %.6f deg\n', s.chain.anchor_angle_deg);
    fprintf(fid, 'delta = %.6f m\n', s.draft);
    fprintf(fid, 'r = %.6f m\n', s.radius);
    fprintf(fid, 'depth residual = %.3e m\n', abs(s.model_depth-18));
    fprintf(fid, 'chain regime = %s\n', char(s.chain.regime));
end
