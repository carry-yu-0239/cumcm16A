%% 问题二：极端风速下的配重调节
% 从仓库根目录运行：matlab -batch "run('src/q2/q2_static_model.m')"
% 模型与 q2_complete.m 保持一致；仅按仓库规范调整脚本和输出路径。

clear; clc; close all;

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

% 参数：与问题一保持一致
p.g = 9.8; p.rho = 1025; p.depth = 18;
p.buoy_radius = 1.0; p.buoy_height = 2.0; p.buoy_mass = 1000;
p.pipe_length = 1.0; p.pipe_radius = 0.025; p.pipe_mass = 10; p.n_pipes = 4;
p.barrel_length = 1.0; p.barrel_radius = 0.15; p.barrel_mass = 100;
p.rho_steel = 7850;
p.chain_length = 22.05; p.chain_link_length = 0.105; p.chain_mass_per_length = 7.0;
p.n_chain_links = round(p.chain_length / p.chain_link_length);
if abs(p.n_chain_links * p.chain_link_length - p.chain_length) > 1e-10
    error('锚链总长不能由给定链环长度整除，请检查参数。');
end
p.pipe_effective_weight = p.pipe_mass*p.g - p.rho*p.g*pi*p.pipe_radius^2*p.pipe_length;
p.barrel_effective_weight = p.barrel_mass*p.g - p.rho*p.g*pi*p.barrel_radius^2*p.barrel_length;
p.chain_effective_weight_per_length = p.chain_mass_per_length*p.g * (1-p.rho/p.rho_steel);

wind_speed = 36;
original_mass = 1200;
selected_mass = 2250;

% 原 1200 kg 方案
state_original = solve_draft(wind_speed, original_mass, p);

% 按 10 kg 步长搜索最小可行质量
mass_grid = (1200:10:3000)';
nm = numel(mass_grid);
barrel_angle_deg = nan(nm,1); anchor_angle_deg = nan(nm,1);
draft_m = nan(nm,1); radius_m = nan(nm,1); feasible = false(nm,1);
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
    error('在 1200--3000 kg 搜索区间内未找到满足两个角度约束的方案。');
end
min_grid_mass = mass_grid(idx);
state_grid = solve_draft(wind_speed, min_grid_mass, p);
mass_scan = table(mass_grid, barrel_angle_deg, anchor_angle_deg, draft_m, radius_m, feasible, ...
    'VariableNames', {'ball_mass_kg','barrel_angle_deg','anchor_angle_deg','buoy_draft_m','swing_radius_m','feasible'});
writetable(mass_scan, fullfile(processed_dir, 'q2_mass_scan.csv'), 'Encoding', 'UTF-8');
writetable(mass_scan, fullfile(tables_dir, 'q2_mass_scan.csv'), 'Encoding', 'UTF-8');
write_table_tex(fullfile(tables_dir, 'q2_mass_scan.tex'), mass_scan);

% 二分搜索理论临界质量与最终推荐方案
critical_mass = solve_critical_mass(wind_speed, p, 1200, 3000);
state_critical = solve_draft(wind_speed, critical_mass, p);
state_selected = solve_draft(wind_speed, selected_mass, p);

scenario = ["原1200kg方案"; "理论临界方案"; "10kg步长最小可行方案"; "最终2250kg方案"];
ball_mass_kg = [original_mass; critical_mass; min_grid_mass; selected_mass];
states = {state_original; state_critical; state_grid; state_selected};
result_table = make_result_table(scenario, ball_mass_kg, states, p.depth);
writetable(result_table, fullfile(tables_dir, 'q2_static_results.csv'), 'Encoding', 'UTF-8');
write_table_tex(fullfile(tables_dir, 'q2_static_results.tex'), result_table);

% 锚链坐标是可复算中间数据。
chain_1200 = chain_coordinates_from_anchor(state_original, p);
chain_2250 = chain_coordinates_from_anchor(state_selected, p);
writetable(chain_1200, fullfile(processed_dir, 'q2_chain_1200kg.csv'), 'Encoding', 'UTF-8');
writetable(chain_2250, fullfile(processed_dir, 'q2_chain_2250kg.csv'), 'Encoding', 'UTF-8');

plot_chain_shapes(chain_1200, chain_2250, p, figures_dir);
plot_mass_constraints(mass_scan, min_grid_mass, selected_mass, figures_dir);
write_summary(fullfile(logs_dir, 'q2_summary.txt'), state_original, state_selected, ...
    original_mass, critical_mass, min_grid_mass, selected_mass);
write_check_log(fullfile(logs_dir, 'q2_static_checks.txt'), states, p);

fprintf('\n问题二计算完成。\n');
fprintf('理论临界质量：%.6f kg\n', critical_mass);
fprintf('10 kg 步长最小可行质量：%.0f kg\n', min_grid_mass);
fprintf('输出目录：%s\n', fullfile(root_dir, 'outputs', 'q2'));

%% Local functions
function state = solve_draft(wind_speed, ball_mass, p)
    draft_grid = linspace(0.30, 1.60, 800);
    prev_valid = false; low = NaN; high = NaN;
    for k = 1:numel(draft_grid)
        delta = draft_grid(k);
        try
            s = evaluate_state(delta, wind_speed, ball_mass, p);
            f = s.model_depth-p.depth;
        catch
            continue;
        end
        if prev_valid && prev_f*f <= 0
            low = prev_delta; high = delta; break;
        end
        prev_valid = true; prev_delta = delta; prev_f = f;
    end
    if isnan(low)
        error('未找到 H(delta)=18 m 的吃水根区间。');
    end
    f_low = evaluate_state(low, wind_speed, ball_mass, p).model_depth-p.depth;
    for k = 1:100
        mid = 0.5*(low+high);
        f_mid = evaluate_state(mid, wind_speed, ball_mass, p).model_depth-p.depth;
        if f_low*f_mid <= 0
            high = mid;
        else
            low = mid; f_low = f_mid;
        end
    end
    state = evaluate_state(0.5*(low+high), wind_speed, ball_mass, p);
end

function state = evaluate_state(delta, wind_speed, ball_mass, p)
    if delta <= 0 || delta >= p.buoy_height
        error('浮标吃水 delta 必须位于 (0,2) m。');
    end
    Fb = p.rho*p.g*pi*p.buoy_radius^2*delta;
    Gb = p.buoy_mass*p.g;
    Fw = 0.625*((p.buoy_height-delta)*2*p.buoy_radius)*wind_speed^2;
    Th = Fw; Tv = Fb-Gb;
    if Th <= 0 || Tv <= 0
        error('当前吃水下顶部张力分量不满足受拉平衡条件。');
    end
    pipe_angles = zeros(p.n_pipes,1);
    rigid_vertical = 0; rigid_horizontal = 0; current_Tv = Tv;
    for i = 1:p.n_pipes
        pipe_angles(i) = rigid_member_angle_from_vertical(Th, current_Tv, p.pipe_effective_weight);
        rigid_vertical = rigid_vertical+p.pipe_length*cos(pipe_angles(i));
        rigid_horizontal = rigid_horizontal+p.pipe_length*sin(pipe_angles(i));
        current_Tv = current_Tv-p.pipe_effective_weight;
    end
    barrel_angle = rigid_member_angle_from_vertical(Th, current_Tv, p.barrel_effective_weight);
    rigid_vertical = rigid_vertical+p.barrel_length*cos(barrel_angle);
    rigid_horizontal = rigid_horizontal+p.barrel_length*sin(barrel_angle);
    current_Tv = current_Tv-p.barrel_effective_weight;
    ball_effective_weight = ball_mass*p.g*(1-p.rho/p.rho_steel);
    chain_top_Tv = current_Tv-ball_effective_weight;
    if chain_top_Tv <= 0
        error('重物球过重或吃水过浅，锚链顶部竖向张力非正。');
    end
    chain = solve_chain_links(Th, chain_top_Tv, p);
    state = struct('wind_speed', wind_speed, 'ball_mass', ball_mass, 'draft', delta, ...
        'wind_force', Fw, 'horizontal_tension', Th, 'top_vertical_tension', Tv, ...
        'pipe_angles', pipe_angles, 'barrel_angle', barrel_angle, ...
        'rigid_vertical_drop', rigid_vertical, 'rigid_horizontal_shift', rigid_horizontal, ...
        'chain_top_vertical_tension', chain_top_Tv, 'chain', chain, ...
        'model_depth', delta+rigid_vertical+chain.vertical_drop, ...
        'radius', rigid_horizontal+chain.horizontal_span);
end

function theta = rigid_member_angle_from_vertical(Th, Tv_top, W)
    Tv_mid = Tv_top-0.5*W;
    if Tv_mid <= 0
        error('刚性构件中点等效竖向张力非正，当前受拉平衡失效。');
    end
    theta = atan2(Th, Tv_mid);
end

function chain = solve_chain_links(Th, Tv_top, p)
    n = p.n_chain_links; ell = p.chain_link_length; w = p.chain_effective_weight_per_length;
    alpha = zeros(n,1); dx = zeros(n,1); dz = zeros(n,1);
    x_from_top = zeros(n+1,1); z_from_top = zeros(n+1,1);
    Tv = Tv_top; suspended_length = 0; seabed_length = 0;
    for j = 1:n
        if Tv <= 1e-12
            alpha(j) = 0; dx(j) = ell; dz(j) = 0; seabed_length = seabed_length+ell; Tv = 0;
        else
            Tv_bottom = Tv-w*ell;
            if Tv_bottom >= 0
                alpha(j) = atan2(0.5*(Tv+Tv_bottom), Th);
                dx(j) = ell*cos(alpha(j)); dz(j) = ell*sin(alpha(j));
                suspended_length = suspended_length+ell; Tv = Tv_bottom;
            else
                ell_s = Tv/w; ell_b = ell-ell_s; alpha(j) = atan2(0.5*Tv, Th);
                dx(j) = ell_s*cos(alpha(j))+ell_b; dz(j) = ell_s*sin(alpha(j));
                suspended_length = suspended_length+ell_s; seabed_length = seabed_length+ell_b; Tv = 0;
            end
        end
        x_from_top(j+1) = x_from_top(j)+dx(j);
        z_from_top(j+1) = z_from_top(j)+dz(j);
    end
    if seabed_length > 1e-10
        anchor_angle_deg = 0; regime = "touchdown_with_seabed_contact";
    else
        anchor_angle_deg = atan2d(Tv, Th); regime = "fully_suspended";
    end
    chain = struct('link_angle_from_horizontal', alpha, 'dx', dx, 'dz', dz, ...
        'x_from_top', x_from_top, 'z_from_top', z_from_top, ...
        'horizontal_span', x_from_top(end), 'vertical_drop', z_from_top(end), ...
        'suspended_length', suspended_length, 'seabed_length', seabed_length, ...
        'bottom_vertical_tension', Tv, 'anchor_angle_deg', anchor_angle_deg, 'regime', regime);
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
        mid_mass = 0.5*(low_mass+high_mass);
        if is_feasible_mass(mid_mass, wind_speed, p)
            high_mass = mid_mass;
        else
            low_mass = mid_mass;
        end
    end
    critical_mass = 0.5*(low_mass+high_mass);
end

function T = make_result_table(scenario, ball_mass_kg, states, depth)
    n = numel(states); values = zeros(n,11); regime = strings(n,1);
    for k = 1:n
        s = states{k}; a = rad2deg(s.pipe_angles);
        values(k,:) = [a', rad2deg(s.barrel_angle), s.chain.anchor_angle_deg, s.draft, ...
            s.radius, s.chain.suspended_length, s.chain.seabed_length, s.model_depth-depth];
        regime(k) = s.chain.regime;
    end
    T = array2table(values, 'VariableNames', {'theta1_deg','theta2_deg','theta3_deg','theta4_deg', ...
        'theta_d_deg','phi_a_deg','delta_m','radius_m','suspended_chain_m','seabed_chain_m','depth_residual_m'});
    T = addvars(T, scenario, ball_mass_kg, regime, 'Before', 1, ...
        'NewVariableNames', {'scenario','ball_mass_kg','chain_regime'});
end

function C = chain_coordinates_from_anchor(state, p)
    x = state.chain.horizontal_span-state.chain.x_from_top;
    z = state.chain.vertical_drop-state.chain.z_from_top;
    C = table((0:p.n_chain_links)', x, z, 'VariableNames', {'node_index','x_m','z_m'});
end

function plot_chain_shapes(C1, C2, p, figures_dir)
    fig = figure('Visible','off','Color','w','Position',[100 100 1050 720]);
    ax = axes(fig); set(ax, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', ...
        'GridColor', [0.75 0.75 0.75], 'FontName', 'Microsoft YaHei');
    hold on; grid on; box on;
    h1 = plot(C1.x_m, C1.z_m, 'LineWidth',2.0);
    h2 = plot(C2.x_m, C2.z_m, 'LineWidth',2.0);
    h3 = yline(0, '-', 'LineWidth',1.0); h4 = yline(p.depth, '--', 'LineWidth',1.0);
    plot(0,0,'ko','MarkerFaceColor','k','MarkerSize',6);
    xlabel('相对锚点的水平距离 / m'); ylabel('距海床高度 / m');
    title('风速 36 m/s 时调节配重前后的锚链形状');
    lgd = legend([h1 h2 h3 h4], {'1200 kg','2250 kg','海床','海面'}, 'Location','best');
    set(lgd, 'Color', 'w', 'TextColor', 'k', 'FontName', 'Microsoft YaHei');
    exportgraphics(fig, fullfile(figures_dir,'q2_chain_shapes.png'), 'Resolution',300);
    exportgraphics(fig, fullfile(figures_dir,'q2_chain_shapes.pdf'), 'Resolution',300);
    close(fig);
end

function plot_mass_constraints(T, min_grid_mass, selected_mass, figures_dir)
    fig = figure('Visible','off','Color','w','Position',[100 100 1050 720]);
    ax = axes(fig); set(ax, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', ...
        'GridColor', [0.75 0.75 0.75], 'FontName', 'Microsoft YaHei');
    hold on; grid on; box on;
    h1 = plot(T.ball_mass_kg, T.barrel_angle_deg, 'LineWidth',2.0);
    h2 = plot(T.ball_mass_kg, T.anchor_angle_deg, 'LineWidth',2.0);
    h3 = yline(5,'--','5 deg'); h4 = yline(16,'--','16 deg');
    xline(min_grid_mass, ':', sprintf('%.0f kg',min_grid_mass));
    xline(selected_mass, ':', sprintf('%.0f kg',selected_mass));
    xlabel('重物球质量 m_s / kg'); ylabel('角度 / deg');
    title('重物球质量对钢桶倾角和锚端夹角的影响');
    lgd = legend([h1 h2 h3 h4], {'钢桶倾角','锚端夹角','钢桶约束','锚端约束'}, 'Location','northeast');
    set(lgd, 'Color', 'w', 'TextColor', 'k', 'FontName', 'Microsoft YaHei');
    exportgraphics(fig, fullfile(figures_dir,'q2_mass_constraints.png'), 'Resolution',300);
    exportgraphics(fig, fullfile(figures_dir,'q2_mass_constraints.pdf'), 'Resolution',300);
    close(fig);
end

function write_summary(path, state_original, state_selected, original_mass, critical_mass, min_grid_mass, selected_mass)
    f = fopen(path, 'w', 'n', 'UTF-8');
    if f < 0, error('无法创建输出摘要文件。'); end
    write_state_text(f, '原 1200 kg 方案', state_original, original_mass);
    fprintf(f, '\n理论临界质量 = %.6f kg\n', critical_mass);
    fprintf(f, '10 kg 步长最小可行质量 = %.0f kg\n', min_grid_mass);
    write_state_text(f, '最终 2250 kg 方案', state_selected, selected_mass);
    fclose(f);
end

function write_check_log(path, states, p)
    f = fopen(path, 'w', 'n', 'UTF-8');
    fprintf(f, 'Generated by src/q2/q2_static_model.m\n');
    fprintf(f, 'pipe_effective_weight_N=%.12f\n', p.pipe_effective_weight);
    fprintf(f, 'barrel_effective_weight_N=%.12f\n', p.barrel_effective_weight);
    fprintf(f, 'chain_effective_weight_N_per_m=%.12f\n', p.chain_effective_weight_per_length);
    for k = 1:numel(states)
        s = states{k};
        fprintf(f, '\nscenario=%d\n', k);
        fprintf(f, 'draft_in_bounds=%d\n', s.draft > 0 && s.draft < p.buoy_height);
        fprintf(f, 'depth_residual_abs_m=%.3e\n', abs(s.model_depth-p.depth));
        fprintf(f, 'chain_length_residual_abs_m=%.3e\n', abs(s.chain.suspended_length+s.chain.seabed_length-p.chain_length));
        fprintf(f, 'barrel_angle_within_5deg=%d\n', rad2deg(s.barrel_angle) <= 5);
        fprintf(f, 'anchor_angle_within_16deg=%d\n', s.chain.anchor_angle_deg <= 16);
        fprintf(f, 'chain_regime=%s\n', char(s.chain.regime));
    end
    fclose(f);
end

function write_state_text(f, name, s, mass)
    a = rad2deg(s.pipe_angles);
    fprintf(f, '\n[%s]\n', name);
    fprintf(f, 'm_s = %.6f kg\n', mass);
    fprintf(f, 'theta_1 = %.6f deg\ntheta_2 = %.6f deg\ntheta_3 = %.6f deg\ntheta_4 = %.6f deg\n', a);
    fprintf(f, 'theta_d = %.6f deg\nphi_a = %.6f deg\ndelta = %.6f m\nr = %.6f m\n', ...
        rad2deg(s.barrel_angle), s.chain.anchor_angle_deg, s.draft, s.radius);
    fprintf(f, 'depth_residual = %.3e m\nchain_regime = %s\n', abs(s.model_depth-18), char(s.chain.regime));
end

function write_table_tex(path, T)
    f = fopen(path, 'w', 'n', 'UTF-8');
    if f < 0, error('无法创建 LaTeX 表格文件。'); end
    ncol = width(T); names = T.Properties.VariableNames;
    fprintf(f, '%% 由 src/q2/q2_static_model.m 自动生成；请勿手工修改。\n');
    fprintf(f, '\\begin{longtable}{%s}\n\\toprule\n', repmat('l', 1, ncol));
    fprintf(f, '%s', tex_escape(names{1}));
    for j = 2:ncol, fprintf(f, ' & %s', tex_escape(names{j})); end
    fprintf(f, ' \\\\\n\\midrule\n');
    for i = 1:height(T)
        for j = 1:ncol
            if j > 1, fprintf(f, ' & '); end
            fprintf(f, '%s', tex_escape(table_value_to_text(T{i,j})));
        end
        fprintf(f, ' \\\\\n');
    end
    fprintf(f, '\\bottomrule\n\\end{longtable}\n');
    fclose(f);
end

function text_value = table_value_to_text(value)
    if islogical(value), text_value = sprintf('%d', value);
    elseif isnumeric(value), text_value = sprintf('%.10g', value);
    elseif isstring(value), text_value = char(value);
    else, text_value = value;
    end
end

function text_value = tex_escape(text_value)
    text_value = strrep(char(text_value), '\', '\\textbackslash{}');
    text_value = strrep(text_value, '_', '\\_');
    text_value = strrep(text_value, '%', '\\%');
    text_value = strrep(text_value, '&', '\\&');
end
