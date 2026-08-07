%% 问题一：静水风载下的浮标—刚体构件—悬链线静力模型
% 从仓库根目录运行：matlab -batch "run('src/q1/q1_static_model.m')"
% 仅使用相对路径；原始数据不作改写。
% 重物球通过无质量竖直短绳悬挂在钢桶—锚链连接点，是分支载荷而非主链构件。

clear; clc;

script_dir = fileparts(mfilename('fullpath'));
root_dir = fullfile(script_dir, '..', '..');
processed_dir = fullfile(root_dir, 'data', 'processed');
tables_dir = fullfile(root_dir, 'outputs', 'q1', 'tables');
figures_dir = fullfile(root_dir, 'outputs', 'q1', 'figures');
logs_dir = fullfile(root_dir, 'outputs', 'q1', 'logs');
output_dirs = {processed_dir, tables_dir, figures_dir, logs_dir};
for k = 1:numel(output_dirs)
    if ~exist(output_dirs{k}, 'dir')
        mkdir(output_dirs{k});
    end
end

% 题目与 paper/sections/03_assumptions.tex 中已确认的参数
c.g = 9.8;
c.rho_water = 1025;
c.depth = 18;
c.buoy_radius = 1;
c.buoy_height = 2;
c.buoy_mass = 1000;
c.pipe_length = 1;
c.pipe_radius = 0.025;
c.pipe_mass = 10;
c.n_pipes = 4;
c.barrel_length = 1;
c.barrel_radius = 0.15;
c.barrel_mass = 100;
c.ball_mass = 1200;
c.steel_density = 7850;
c.chain_length = 22.05;
c.chain_mass_per_length = 7;

% 密闭圆柱的外轮廓排水浮力；锚链由材料密度折算排水体积。
c.pipe_effective_weight = c.pipe_mass*c.g - cylinder_buoyancy(c.pipe_radius, c.pipe_length, c);
c.barrel_effective_weight = c.barrel_mass*c.g - cylinder_buoyancy(c.barrel_radius, c.barrel_length, c);
c.ball_effective_weight = c.ball_mass*c.g*(1-c.rho_water/c.steel_density);
c.chain_effective_weight_per_length = c.chain_mass_per_length*c.g*(1-c.rho_water/c.steel_density);

wind_speeds = [12, 24];
states = cell(size(wind_speeds));
for k = 1:numel(wind_speeds)
    states{k} = solve_draft(wind_speeds(k), c);
end

% 函数模型 H(delta)：分别输出两种风速的离散曲线数据。
curve_rows = table();
for k = 1:numel(wind_speeds)
    % 0.67 m 以下不足以使重物球支线与锚链同时保持张紧，故不属于本模型的可行域。
    draft_grid = (0.670:0.001:0.850)';
    h_grid = zeros(size(draft_grid));
    for j = 1:numel(draft_grid)
        state = evaluate_state(draft_grid(j), wind_speeds(k), c);
        h_grid(j) = state.model_depth;
    end
    curve_rows = [curve_rows; table(repmat(wind_speeds(k), numel(draft_grid), 1), draft_grid, h_grid, ...
        'VariableNames', {'wind_speed_mps', 'draft_m', 'model_depth_m'})]; %#ok<AGROW>
end
writetable(curve_rows, fullfile(processed_dir, 'q1_depth_vs_draft.csv'), 'Encoding', 'UTF-8');

% H(delta)=18 m 的精确求解点单独保存，供图形和结果复核直接引用。
solution_drafts = cellfun(@(s) s.draft, states)';
solution_depths = cellfun(@(s) s.model_depth, states)';
solution_rows = table(wind_speeds(:), solution_drafts, solution_depths, c.depth*ones(numel(wind_speeds), 1), ...
    'VariableNames', {'wind_speed_mps', 'draft_m', 'model_depth_m', 'target_depth_m'});
writetable(solution_rows, fullfile(processed_dir, 'q1_depth_draft_solutions.csv'), 'Encoding', 'UTF-8');

% 汇总指标表（CSV 和 LaTeX 同时输出）。
summary = make_summary_table(states, wind_speeds);
writetable(summary, fullfile(tables_dir, 'q1_static_results.csv'), 'Encoding', 'UTF-8');
write_summary_tex(fullfile(tables_dir, 'q1_static_results.tex'), states, wind_speeds);

% 锚链坐标为中间数据，便于复算曲线；坐标原点为锚点，z 轴向上。
chain_tables = cell(size(states));
for k = 1:numel(states)
    chain_tables{k} = chain_coordinates(states{k}, c);
    writetable(chain_tables{k}, fullfile(processed_dir, sprintf('q1_chain_shape_%dms.csv', wind_speeds(k))), 'Encoding', 'UTF-8');
end

plot_depth_curve(curve_rows, solution_rows, wind_speeds, c, figures_dir);
plot_chain_shapes(chain_tables, wind_speeds, c, figures_dir);
write_check_log(fullfile(logs_dir, 'q1_static_checks.txt'), states, c);

%% Local functions
function Fb = cylinder_buoyancy(radius, length, c)
    Fb = c.rho_water*c.g*pi*radius^2*length;
end

function state = evaluate_state(draft, wind_speed, c)
    if draft < 0 || draft > c.buoy_height
        error('浮标吃水必须在 [0, 2] m 内。');
    end
    horizontal_tension = 0.625 * (2*(c.buoy_height-draft)) * wind_speed^2;
    vertical_top = c.rho_water*c.g*pi*c.buoy_radius^2*draft - c.buoy_mass*c.g;
    pipe_angles = zeros(c.n_pipes, 1);
    rigid_vertical_drop = 0;
    rigid_horizontal_shift = 0;
    current_vertical = vertical_top;
    for i = 1:c.n_pipes
        pipe_angles(i) = member_angle(horizontal_tension, current_vertical, c.pipe_effective_weight);
        rigid_vertical_drop = rigid_vertical_drop + c.pipe_length*cos(pipe_angles(i));
        rigid_horizontal_shift = rigid_horizontal_shift + c.pipe_length*sin(pipe_angles(i));
        current_vertical = current_vertical - c.pipe_effective_weight;
    end
    barrel_angle = member_angle(horizontal_tension, current_vertical, c.barrel_effective_weight);
    rigid_vertical_drop = rigid_vertical_drop + c.barrel_length*cos(barrel_angle);
    rigid_horizontal_shift = rigid_horizontal_shift + c.barrel_length*sin(barrel_angle);

    % 钢桶下端节点：锚链张力与重物球支线拉力共同由钢桶承受。
    chain_top_vertical = current_vertical - c.barrel_effective_weight - c.ball_effective_weight;
    chain = chain_geometry(horizontal_tension, chain_top_vertical, c);
    state = struct( ...
        'draft', draft, 'wind_speed', wind_speed, ...
        'horizontal_tension', horizontal_tension, ...
        'buoy_vertical_tension', vertical_top, ...
        'chain_top_vertical', chain_top_vertical, ...
        'pipe_angles', pipe_angles, 'barrel_angle', barrel_angle, ...
        'rigid_vertical_drop', rigid_vertical_drop, ...
        'rigid_horizontal_shift', rigid_horizontal_shift, ...
        'chain', chain, ...
        'model_depth', draft + rigid_vertical_drop + chain.vertical_drop, ...
        'radius', rigid_horizontal_shift + chain.horizontal_span + chain.seabed_length);
end

function angle = member_angle(horizontal_tension, top_vertical_tension, effective_weight)
    denominator = top_vertical_tension - effective_weight/2;
    if denominator <= 0
        error('构件顶端竖向张力不足，刚体受拉平衡假设失效。');
    end
    % 对匀质刚体圆柱上端取矩：tan(phi)=H/(V_top-W/2)。
    angle = atan2(horizontal_tension, denominator);
end

function chain = chain_geometry(horizontal_tension, top_vertical_tension, c)
    if horizontal_tension <= 0 || top_vertical_tension <= 0
        error('悬链线顶端张力分量必须为正。');
    end
    required_length_to_touchdown = top_vertical_tension/c.chain_effective_weight_per_length;
    if required_length_to_touchdown <= c.chain_length
        lower_vertical_tension = 0;
        suspended_length = required_length_to_touchdown;
        seabed_length = c.chain_length-suspended_length;
        regime = 'touchdown_with_seabed_contact';
    else
        lower_vertical_tension = top_vertical_tension-c.chain_effective_weight_per_length*c.chain_length;
        suspended_length = c.chain_length;
        seabed_length = 0;
        regime = 'fully_suspended';
    end
    top_norm = hypot(horizontal_tension, top_vertical_tension);
    lower_norm = hypot(horizontal_tension, lower_vertical_tension);
    chain.vertical_drop = (top_norm-lower_norm)/c.chain_effective_weight_per_length;
    chain.horizontal_span = horizontal_tension/c.chain_effective_weight_per_length * ...
        (asinh(top_vertical_tension/horizontal_tension)-asinh(lower_vertical_tension/horizontal_tension));
    chain.suspended_length = suspended_length;
    chain.seabed_length = seabed_length;
    chain.lower_vertical_tension = lower_vertical_tension;
    chain.lower_angle_deg = atan2d(lower_vertical_tension, horizontal_tension);
    chain.regime = regime;
end

function state = solve_draft(wind_speed, c)
    low = 0.67; high = 0.85;
    f_low = evaluate_state(low, wind_speed, c).model_depth-c.depth;
    f_high = evaluate_state(high, wind_speed, c).model_depth-c.depth;
    if ~(f_low < 0 && f_high > 0)
        error('吃水求根区间未包围解：f(low)=%.6g, f(high)=%.6g。', f_low, f_high);
    end
    for k = 1:80
        mid = (low+high)/2;
        if evaluate_state(mid, wind_speed, c).model_depth < c.depth
            low = mid;
        else
            high = mid;
        end
    end
    state = evaluate_state((low+high)/2, wind_speed, c);
end

function summary = make_summary_table(states, wind_speeds)
    n = numel(states);
    values = zeros(n, 12);
    regime = strings(n, 1);
    for k = 1:n
        s = states{k};
        values(k,:) = [wind_speeds(k), rad2deg(s.pipe_angles)', rad2deg(s.barrel_angle), ...
            s.chain.lower_angle_deg, s.draft, s.radius, ...
            s.chain.suspended_length, s.chain.seabed_length, s.model_depth-18];
        regime(k) = string(s.chain.regime);
    end
    summary = array2table(values, 'VariableNames', {'wind_speed_mps', 'pipe_1_angle_deg', 'pipe_2_angle_deg', ...
        'pipe_3_angle_deg', 'pipe_4_angle_deg', 'barrel_angle_deg', 'last_chain_angle_deg', ...
        'buoy_draft_m', 'swing_radius_m', 'suspended_chain_length_m', 'seabed_chain_length_m', 'depth_residual_m'});
    summary.chain_regime = regime;
end

function write_summary_tex(path, states, wind_speeds)
    labels = {'第1节钢管与海平面竖直方向夹角（度）', '第2节钢管与海平面竖直方向夹角（度）', ...
        '第3节钢管与海平面竖直方向夹角（度）', '第4节钢管与海平面竖直方向夹角（度）', ...
        '钢桶与海平面竖直方向夹角（度）', '最后一条锚链与海平面水平方向夹角（度）', ...
        '浮标的吃水深度（m）', '浮标的游动半径（m）（以锚所在位置为圆心）'};
    f = fopen(path, 'w', 'n', 'UTF-8');
    fprintf(f, '%% 由 src/q1/q1_static_model.m 自动生成；请勿手工修改。\n');
    fprintf(f, '\\begin{tabular}{lrr}\n\\toprule\n指标 & %d m/s & %d m/s \\\\ \n\\midrule\n', wind_speeds(1), wind_speeds(2));
    for i = 1:numel(labels)
        a = result_value(states{1}, i); b = result_value(states{2}, i);
        fprintf(f, '%s & %.3f & %.3f \\\\ \n', labels{i}, a, b);
    end
    fprintf(f, '\\bottomrule\n\\end{tabular}\n');
    fclose(f);
end

function x = result_value(s, i)
    values = [rad2deg(s.pipe_angles)', rad2deg(s.barrel_angle), s.chain.lower_angle_deg, s.draft, s.radius];
    x = values(i);
end

function coords = chain_coordinates(state, c)
    chain = state.chain;
    n_suspended = 241;
    arc = linspace(0, chain.suspended_length, n_suspended)';
    h = state.horizontal_tension;
    vt = state.chain_top_vertical;
    local_v = vt-c.chain_effective_weight_per_length*arc;
    x_from_top = h/c.chain_effective_weight_per_length * (asinh(vt/h)-asinh(local_v/h));
    z_drop = (hypot(h,vt)-hypot(h,local_v))/c.chain_effective_weight_per_length;
    top_x = state.radius-state.rigid_horizontal_shift;
    top_z = -state.draft-state.rigid_vertical_drop;
    coords = table(repmat("suspended", n_suspended, 1), arc, top_x-x_from_top, top_z-z_drop, ...
        'VariableNames', {'segment', 'arc_length_m', 'x_m', 'z_m'});
    if chain.seabed_length > 0
        bed_arc = chain.suspended_length + linspace(chain.seabed_length/40, chain.seabed_length, 40)';
        bed_x = linspace(chain.seabed_length*(39/40), 0, 40)';
        bed = table(repmat("seabed", 40, 1), bed_arc, bed_x, -c.depth*ones(40,1), ...
            'VariableNames', {'segment', 'arc_length_m', 'x_m', 'z_m'});
        coords = [coords; bed];
    end
end

function plot_depth_curve(curve_rows, solution_rows, wind_speeds, c, figures_dir)
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100, 100, 1320, 840]);
    ax = axes(fig); set(ax, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'GridColor', [0.75 0.75 0.75], ...
        'FontName', 'Microsoft YaHei', 'FontSize', 11);
    hold on; grid on; box on;
    colors = [0.12, 0.47, 0.71; 0.77, 0.31, 0.32];
    h_curve = gobjects(size(wind_speeds));
    h_solution = gobjects(size(wind_speeds));
    for k = 1:numel(wind_speeds)
        rows = curve_rows.wind_speed_mps == wind_speeds(k);
        h_curve(k) = plot(curve_rows.draft_m(rows), curve_rows.model_depth_m(rows), ...
            'LineWidth', 2, 'Color', colors(k,:));
        h_solution(k) = plot(solution_rows.draft_m(k), solution_rows.model_depth_m(k), 'ko', ...
            'MarkerSize', 7, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
    end
    h_target = yline(c.depth, '--k', 'LineWidth', 1.2);
    text(solution_rows.draft_m(1)-0.0010, c.depth+0.75, ...
        sprintf('\\delta_{12}=%.3f m', solution_rows.draft_m(1)), ...
        'Interpreter', 'tex', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', ...
        'BackgroundColor', 'w', 'Margin', 2);
    text(solution_rows.draft_m(2)+0.0010, c.depth-0.85, ...
        sprintf('\\delta_{24}=%.3f m', solution_rows.draft_m(2)), ...
        'Interpreter', 'tex', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
        'BackgroundColor', 'w', 'Margin', 2);
    title('不同风速下浮标吃水与模型水深的关系', 'FontName', 'Microsoft YaHei', 'Interpreter', 'none');
    xlabel('浮标吃水 \delta（m）', 'FontName', 'Microsoft YaHei');
    ylabel('模型水深 H(\delta)（m）', 'FontName', 'Microsoft YaHei');
    lgd = legend([h_curve(1), h_curve(2), h_target, h_solution(1), h_solution(2)], ...
        {'风速 12 m/s', '风速 24 m/s', '目标水深 18 m', '12 m/s 的求解吃水', '24 m/s 的求解吃水'}, ...
        'Location', 'southeast');
    set(lgd, 'Color', 'w', 'TextColor', 'k', 'FontName', 'Microsoft YaHei');
    xlim([0.67, 0.85]);
    print(fig, fullfile(figures_dir, 'q1_depth_vs_draft.png'), '-dpng', '-r300');
    print(fig, fullfile(figures_dir, 'q1_depth_vs_draft.pdf'), '-dpdf', '-r300');
    close(fig);
end

function plot_chain_shapes(chain_tables, wind_speeds, c, figures_dir)
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100, 100, 1320, 1000]);
    ax = axes(fig); set(ax, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'GridColor', [0.75 0.75 0.75], ...
        'FontName', 'Microsoft YaHei', 'FontSize', 11);
    hold on; grid on; box on;
    colors = [0.12, 0.47, 0.71; 0.77, 0.31, 0.32];
    h_chain = gobjects(size(wind_speeds));
    for k = 1:numel(wind_speeds)
        h_chain(k) = plot(chain_tables{k}.x_m, chain_tables{k}.z_m, 'LineWidth', 2.2, 'Color', colors(k,:));
    end
    h_surface = yline(0, '-', 'Color', [0.3, 0.55, 0.75], 'LineWidth', 1.2);
    h_bed = yline(-c.depth, '-', 'Color', [0.55, 0.42, 0.34], 'LineWidth', 1.2);
    h_anchor = plot(0, -c.depth, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6);
    title('不同风速下锚链的静力平衡形状', 'FontName', 'Microsoft YaHei', 'Interpreter', 'none');
    xlabel('相对锚点的水平坐标（m）', 'FontName', 'Microsoft YaHei');
    ylabel('高程（m）', 'FontName', 'Microsoft YaHei');
    text(0.5, -1.0, '锚点为坐标原点，z 轴向上', 'FontName', 'Microsoft YaHei', ...
        'BackgroundColor', 'w', 'Margin', 2);
    axis equal;
    lgd = legend([h_chain(1), h_chain(2), h_surface, h_bed, h_anchor], ...
        {'风速 12 m/s', '风速 24 m/s', '静水面', '海床（-18 m）', '锚点'}, 'Location', 'none');
    set(lgd, 'Color', 'w', 'TextColor', 'k', 'FontName', 'Microsoft YaHei', ...
        'Units', 'normalized', 'Position', [0.14, 0.58, 0.28, 0.26]);
    print(fig, fullfile(figures_dir, 'q1_chain_shapes.png'), '-dpng', '-r300');
    print(fig, fullfile(figures_dir, 'q1_chain_shapes.pdf'), '-dpdf', '-r300');
    close(fig);
end

function write_check_log(path, states, c)
    f = fopen(path, 'w', 'n', 'UTF-8');
    fprintf(f, 'Generated by src/q1/q1_static_model.m\n');
    fprintf(f, 'pipe_effective_weight_N=%.12f\n', c.pipe_effective_weight);
    fprintf(f, 'barrel_effective_weight_N=%.12f\n', c.barrel_effective_weight);
    fprintf(f, 'ball_effective_weight_N=%.12f\n', c.ball_effective_weight);
    fprintf(f, 'chain_effective_weight_N_per_m=%.12f\n', c.chain_effective_weight_per_length);
    for k = 1:numel(states)
        s = states{k};
        fprintf(f, '\nwind_speed_mps=%.0f\n', s.wind_speed);
        fprintf(f, 'draft_in_bounds=%d\n', s.draft >= 0 && s.draft <= c.buoy_height);
        fprintf(f, 'depth_residual_abs_m=%.3e\n', abs(s.model_depth-c.depth));
        fprintf(f, 'chain_length_residual_abs_m=%.3e\n', abs(s.chain.suspended_length+s.chain.seabed_length-c.chain_length));
        fprintf(f, 'anchor_angle_within_16deg=%d\n', s.chain.lower_angle_deg <= 16);
        fprintf(f, 'chain_regime=%s\n', s.chain.regime);
    end
    fclose(f);
end
