function q1_static_equilibrium()
%Q1_STATIC_EQUILIBRIUM Reproducible working model for CUMCM 2016 A, Q1.
% Run this file from any working directory. Outputs are written under outputs/q1.

    root_dir = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    output_dir = fullfile(root_dir, 'outputs', 'q1');
    table_dir = fullfile(output_dir, 'tables');
    figure_dir = fullfile(output_dir, 'figures');
    log_dir = fullfile(output_dir, 'logs');
    make_directory(table_dir);
    make_directory(figure_dir);
    make_directory(log_dir);

    p = build_parameters(root_dir);
    winds = [12, 24];
    depth_target = 18;
    results = repmat(struct(), size(winds));
    for k = 1:numel(winds)
        results(k) = solve_condition(p, winds(k), depth_target);
    end

    write_parameter_table(p, table_dir);
    write_summary_table(results, table_dir);
    write_depth_curve(p, winds, depth_target, results, table_dir, figure_dir);
    write_chain_figure(p, results, figure_dir);
    write_log(p, results, depth_target, log_dir);
end

function p = build_parameters(root_dir)
    p.rho = 1025;
    p.g = 9.8;
    p.water_depth = 18;
    p.float_radius = 1.0;
    p.float_height = 2.0;
    p.float_diameter = 2.0;
    p.float_mass = 1000;
    p.pipe_length = 1.0;
    p.pipe_radius = 0.025;
    p.pipe_mass = 10;
    p.barrel_length = 1.0;
    p.barrel_radius = 0.15;
    p.barrel_mass = 100;
    p.ball_mass = 1200;
    p.steel_density = 7850;
    p.chain_length = 22.05;
    p.chain_type = "II";

    chain_file = fullfile(root_dir, 'problem', 'original_problem', ...
        '附表-锚链型号和参数表.csv');
    options = detectImportOptions(chain_file, 'Encoding', 'UTF-8');
    chain_data = readtable(chain_file, options);
    chain_type = string(chain_data{:, 1});
    row = find(chain_type == p.chain_type, 1);
    assert(~isempty(row), 'Q1:ChainType', 'II 型锚链未在附表中找到。');
    p.chain_mass_per_length = chain_data{row, 3};

    % Conditional assumption OQ-1: every steel pipe is a fully displacing outer cylinder.
    p.pipe_displaced_volume = pi * p.pipe_radius^2 * p.pipe_length;
    p.barrel_displaced_volume = pi * p.barrel_radius^2 * p.barrel_length;
    p.ball_displaced_volume = p.ball_mass / p.steel_density;
    p.pipe_effective_weight = p.g * (p.pipe_mass - p.rho * p.pipe_displaced_volume);
    p.barrel_effective_weight = p.g * (p.barrel_mass - p.rho * p.barrel_displaced_volume);
    p.ball_effective_weight = p.g * (p.ball_mass - p.rho * p.ball_displaced_volume);
    p.chain_effective_weight_per_length = p.g * p.chain_mass_per_length * ...
        (1 - p.rho / p.steel_density);
    assert(p.pipe_effective_weight > 0 && p.barrel_effective_weight > 0 && ...
        p.ball_effective_weight > 0 && p.chain_effective_weight_per_length > 0, ...
        'Q1:EffectiveWeight', '所有水中有效重力必须为正。');
end

function result = solve_condition(p, wind_speed, target_depth)
    residual = @(draft) geometry_at_draft(p, wind_speed, draft).water_depth - target_depth;
    lower = (p.float_mass * p.g + 4 * p.pipe_effective_weight + ...
        p.barrel_effective_weight + p.ball_effective_weight + 1e-6) / ...
        (p.rho * p.g * pi * p.float_radius^2);
    upper = p.float_height - 1e-6;
    assert(lower < upper, 'Q1:DraftBracket', '物理吃水区间为空。');
    assert(residual(lower) < 0 && residual(upper) > 0, ...
        'Q1:RootBracket', 'H(delta)=18 的求根区间未被包围。');
    draft = fzero(residual, [lower, upper]);
    result = geometry_at_draft(p, wind_speed, draft);
    result.depth_residual = result.water_depth - target_depth;
    assert(abs(result.depth_residual) < 1e-8, 'Q1:DepthResidual', ...
        '水深闭合残差超过容差。');
    assert(result.chain.arc_length_residual < 1e-8, 'Q1:ArcLength', ...
        '锚链弧长校验失败。');
end

function state = geometry_at_draft(p, wind_speed, draft)
    assert(draft > 0 && draft < p.float_height, 'Q1:DraftRange', ...
        '吃水必须在 (0, 2) m 内。');
    horizontal_tension = 0.625 * (p.float_height - draft) * ...
        p.float_diameter * wind_speed^2;
    vertical_at_float = p.rho * p.g * pi * p.float_radius^2 * draft - p.float_mass * p.g;

    pipe_angle = zeros(1, 4);
    vertical_top = vertical_at_float;
    for i = 1:4
        denominator = vertical_top - p.pipe_effective_weight / 2;
        assert(denominator > 0, 'Q1:PipeTension', '钢管平均竖向张力非正。');
        pipe_angle(i) = atan(horizontal_tension / denominator);
        vertical_top = vertical_top - p.pipe_effective_weight;
    end
    barrel_vertical_top = vertical_top;
    barrel_denominator = barrel_vertical_top - p.barrel_effective_weight / 2;
    assert(barrel_denominator > 0, 'Q1:BarrelTension', '钢桶平均竖向张力非正。');
    barrel_angle = atan(horizontal_tension / barrel_denominator);
    chain_vertical_top = barrel_vertical_top - p.barrel_effective_weight - p.ball_effective_weight;
    assert(chain_vertical_top > 0, 'Q1:ChainTension', '钢桶端锚链竖向张力非正。');

    chain = chain_geometry(horizontal_tension, chain_vertical_top, p.chain_effective_weight_per_length, p.chain_length);
    state.wind_speed = wind_speed;
    state.draft = draft;
    state.horizontal_tension = horizontal_tension;
    state.vertical_at_float = vertical_at_float;
    state.pipe_angle_rad = pipe_angle;
    state.pipe_angle_deg = rad2deg(pipe_angle);
    state.barrel_angle_rad = barrel_angle;
    state.barrel_angle_deg = rad2deg(barrel_angle);
    state.chain_anchor_angle_rad = chain.anchor_angle_rad;
    state.chain_anchor_angle_deg = rad2deg(chain.anchor_angle_rad);
    state.chain = chain;
    state.swim_radius = chain.horizontal_span + p.pipe_length * sum(sin(pipe_angle)) + ...
        p.barrel_length * sin(barrel_angle);
    state.water_depth = draft + p.pipe_length * sum(cos(pipe_angle)) + ...
        p.barrel_length * cos(barrel_angle) + chain.vertical_span;
end

function chain = chain_geometry(horizontal, vertical_top, q, total_length)
    assert(horizontal > 0 && vertical_top > 0 && q > 0 && total_length > 0, ...
        'Q1:ChainInput', '锚链几何参数必须为正。');
    vertical_capacity = q * total_length;
    if vertical_top <= vertical_capacity
        suspended_length = vertical_top / q;
        bed_length = total_length - suspended_length;
        vertical_span = (hypot(horizontal, vertical_top) - horizontal) / q;
        suspended_horizontal = horizontal / q * asinh(vertical_top / horizontal);
        horizontal_span = bed_length + suspended_horizontal;
        anchor_vertical = 0;
        mode = "partly_grounded";
    else
        suspended_length = total_length;
        bed_length = 0;
        anchor_vertical = vertical_top - vertical_capacity;
        vertical_span = (hypot(horizontal, vertical_top) - hypot(horizontal, anchor_vertical)) / q;
        horizontal_span = horizontal / q * (asinh(vertical_top / horizontal) - ...
            asinh(anchor_vertical / horizontal));
        mode = "fully_suspended";
    end
    anchor_angle_rad = atan(anchor_vertical / horizontal);
    chain.mode = mode;
    chain.suspended_length = suspended_length;
    chain.bed_length = bed_length;
    chain.vertical_span = vertical_span;
    chain.horizontal_span = horizontal_span;
    chain.anchor_vertical = anchor_vertical;
    chain.anchor_angle_rad = anchor_angle_rad;
    chain.arc_length_residual = abs(suspended_length - ...
        (vertical_top - anchor_vertical) / q);
end

function write_parameter_table(p, table_dir)
    values = [p.rho; p.g; p.chain_length; p.chain_mass_per_length; ...
        p.pipe_effective_weight; p.barrel_effective_weight; p.ball_effective_weight; ...
        p.chain_effective_weight_per_length];
    names = ["sea_water_density"; "gravity"; "II_chain_total_length"; ...
        "II_chain_mass_per_length"; "pipe_effective_weight"; "barrel_effective_weight"; ...
        "ball_effective_weight"; "chain_effective_weight_per_length"];
    units = ["kg/m^3"; "m/s^2"; "m"; "kg/m"; "N"; "N"; "N"; "N/m"];
    T = table(names, values, units, 'VariableNames', {"parameter", "value", "unit"});
    writetable(T, fullfile(table_dir, 'q1_input_parameters.csv'), 'Encoding', 'UTF-8');
    write_latex_parameters(fullfile(table_dir, 'q1_input_parameters.tex'), T);
end

function write_summary_table(results, table_dir)
    metrics = ["第1节钢管与竖直方向夹角（度）"; "第2节钢管与竖直方向夹角（度）"; ...
        "第3节钢管与竖直方向夹角（度）"; "第4节钢管与竖直方向夹角（度）"; ...
        "钢桶与竖直方向夹角（度）"; "锚端锚链与水平方向夹角（度）"; ...
        "浮标的吃水深度（m）"; "浮标的游动半径（m）"];
    values_12 = [results(1).pipe_angle_deg(:); results(1).barrel_angle_deg; ...
        results(1).chain_anchor_angle_deg; results(1).draft; results(1).swim_radius];
    values_24 = [results(2).pipe_angle_deg(:); results(2).barrel_angle_deg; ...
        results(2).chain_anchor_angle_deg; results(2).draft; results(2).swim_radius];
    T = table(metrics, values_12, values_24, 'VariableNames', {"metric", "wind_12_mps", "wind_24_mps"});
    writetable(T, fullfile(table_dir, 'q1_summary.csv'), 'Encoding', 'UTF-8');
    write_latex_summary(fullfile(table_dir, 'q1_summary.tex'), T);
end

function write_depth_curve(p, winds, target_depth, results, table_dir, figure_dir)
    minimum_draft = (p.float_mass * p.g + 4 * p.pipe_effective_weight + ...
        p.barrel_effective_weight + p.ball_effective_weight) / ...
        (p.rho * p.g * pi * p.float_radius^2);
    draft_grid = ((minimum_draft + 1e-4):0.001:1.20)';
    H12 = nan(size(draft_grid));
    H24 = nan(size(draft_grid));
    for i = 1:numel(draft_grid)
        H12(i) = geometry_at_draft(p, winds(1), draft_grid(i)).water_depth;
        H24(i) = geometry_at_draft(p, winds(2), draft_grid(i)).water_depth;
    end
    T = table(draft_grid, H12, H24, 'VariableNames', {"draft_m", "depth_at_12_mps_m", "depth_at_24_mps_m"});
    writetable(T, fullfile(table_dir, 'q1_depth_draft_curve.csv'), 'Encoding', 'UTF-8');
    write_latex_curve(fullfile(table_dir, 'q1_depth_draft_curve.tex'), T);

    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100, 100, 840, 560]);
    plot(draft_grid, H12, 'LineWidth', 1.8); hold on;
    plot(draft_grid, H24, 'LineWidth', 1.8);
    yline(target_depth, '--k', 'LineWidth', 1.2);
    plot(results(1).draft, target_depth, 'o', 'MarkerSize', 7, 'LineWidth', 1.4);
    plot(results(2).draft, target_depth, 's', 'MarkerSize', 7, 'LineWidth', 1.4);
    grid on; box on;
    xlabel('Draft \delta (m)'); ylabel('Computed water depth H(\delta) (m)');
    title('Depth--draft curves under static wind loading');
    legend('v = 12 m/s', 'v = 24 m/s', 'H = 18 m', 'Location', 'northwest');
    exportgraphics(fig, fullfile(figure_dir, 'q1_depth_draft_curve.png'), 'Resolution', 300);
    close(fig);
end

function write_chain_figure(p, results, figure_dir)
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100, 100, 840, 560]);
    colors = lines(numel(results));
    hold on;
    max_x = 0;
    for k = 1:numel(results)
        [x, z] = chain_profile(results(k).horizontal_tension, results(k).vertical_at_float - ...
            4 * p.pipe_effective_weight - p.barrel_effective_weight - p.ball_effective_weight, ...
            p.chain_effective_weight_per_length, p.chain_length);
        plot(x, z, 'LineWidth', 1.8, 'Color', colors(k, :));
        max_x = max(max_x, max(x));
    end
    yline(0, 'k-', 'LineWidth', 1.0);
    grid on; box on;
    xlim([0, max_x + 1]); ylim([-0.25, 14]);
    xlabel('Horizontal distance from anchor (m)'); ylabel('Height above seabed (m)');
    title('Computed catenary-chain shapes');
    legend('v = 12 m/s', 'v = 24 m/s', 'Seabed', 'Location', 'northwest');
    exportgraphics(fig, fullfile(figure_dir, 'q1_chain_shapes.png'), 'Resolution', 300);
    close(fig);
end

function [x, z] = chain_profile(horizontal, vertical_top, q, total_length)
    chain = chain_geometry(horizontal, vertical_top, q, total_length);
    if chain.bed_length > 0
        x_bed = linspace(0, chain.bed_length, 40);
        z_bed = zeros(size(x_bed));
    else
        x_bed = 0;
        z_bed = 0;
    end
    vertical = linspace(chain.anchor_vertical, vertical_top, 200);
    x_suspended = chain.bed_length + horizontal / q * ...
        (asinh(vertical / horizontal) - asinh(chain.anchor_vertical / horizontal));
    z_suspended = (hypot(horizontal, vertical) - hypot(horizontal, chain.anchor_vertical)) / q;
    x = [x_bed, x_suspended];
    z = [z_bed, z_suspended];
end

function write_log(p, results, target_depth, log_dir)
    file = fullfile(log_dir, 'q1_run_summary.txt');
    fid = fopen(file, 'w', 'n', 'UTF-8');
    assert(fid >= 0, 'Q1:Log', '无法创建日志文件。');
    cleaner = onCleanup(@() fclose(fid)); %#ok<NASGU>
    fprintf(fid, 'Q1 static equilibrium run\n');
    fprintf(fid, 'Target water depth: %.6f m\n', target_depth);
    fprintf(fid, 'Pipe buoyancy condition: fully displacing outer cylinder, radius %.6f m\n', p.pipe_radius);
    for k = 1:numel(results)
        r = results(k);
        fprintf(fid, '\nWind %.0f m/s\n', r.wind_speed);
        fprintf(fid, 'draft=%.12f, H=%.12f, residual=%.3e\n', r.draft, r.water_depth, r.depth_residual);
        fprintf(fid, 'chain_mode=%s, suspended_length=%.12f, bed_length=%.12f\n', ...
            r.chain.mode, r.chain.suspended_length, r.chain.bed_length);
        fprintf(fid, 'arc_length_residual=%.3e\n', r.chain.arc_length_residual);
    end
end

function write_latex_parameters(filename, T)
    fid = open_latex_file(filename);
    fprintf(fid, '%% Generated by src/q1/q1_static_equilibrium.m\\n');
    fprintf(fid, '\\begin{tabular}{lrr}\\toprule\\nParameter & Value & Unit \\\\ \\midrule\\n');
    for i = 1:height(T)
        fprintf(fid, '%s & %.9f & %s \\\\n', T.parameter(i), T.value(i), T.unit(i));
    end
    fprintf(fid, '\\bottomrule\\n\\end{tabular}\\n');
    fclose(fid);
end

function write_latex_summary(filename, T)
    fid = open_latex_file(filename);
    fprintf(fid, '%% Generated by src/q1/q1_static_equilibrium.m\\n');
    fprintf(fid, '\\begin{tabular}{lrr}\\toprule\\n指标 & 12 m/s & 24 m/s \\\\ \\midrule\\n');
    for i = 1:height(T)
        fprintf(fid, '%s & %.6f & %.6f \\\\n', T.metric(i), T.wind_12_mps(i), T.wind_24_mps(i));
    end
    fprintf(fid, '\\bottomrule\\n\\end{tabular}\\n');
    fclose(fid);
end

function write_latex_curve(filename, T)
    fid = open_latex_file(filename);
    fprintf(fid, '%% Generated by src/q1/q1_static_equilibrium.m\\n');
    fprintf(fid, '\\begin{tabular}{rrr}\\toprule\\n$\\delta$ (m) & $H_{12}$ (m) & $H_{24}$ (m) \\\\ \\midrule\\n');
    for i = 1:height(T)
        fprintf(fid, '%.3f & %.6f & %.6f \\\\n', T.draft_m(i), T.depth_at_12_mps_m(i), T.depth_at_24_mps_m(i));
    end
    fprintf(fid, '\\bottomrule\\n\\end{tabular}\\n');
    fclose(fid);
end

function fid = open_latex_file(filename)
    fid = fopen(filename, 'w', 'n', 'UTF-8');
    assert(fid >= 0, 'Q1:Latex', '无法创建 LaTeX 表格文件。');
end

function make_directory(folder)
    if ~isfolder(folder)
        mkdir(folder);
    end
end
