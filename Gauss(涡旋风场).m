%% 三维区域地形模型 - 旋转风场
% 加载区域地形数据
load topo;  % MATLAB内置全球地形数据

% 定义区域范围
lat_range = [0, 70];  % 纬度范围
lon_range = [60, 150]; % 经度范围

% 计算索引范围 (topo数据: 第一行90°N，最后一行90°S)
row_start = 90 - lat_range(2) + 1;
row_end = 90 - lat_range(1) + 1;
col_start = lon_range(1) + 1;
col_end = lon_range(2) + 1;

% 确保索引在有效范围内
row_start = max(1, min(row_start, 180));
row_end = max(1, min(row_end, 180));
col_start = max(1, min(col_start, 360));
col_end = max(1, min(col_end, 360));

% 提取区域地形
Z_china = topo(row_start:row_end, col_start:col_end);
[nrows, ncols] = size(Z_china);

% 创建坐标网格 (单位：米)
x_range = linspace(0, 1000000, ncols);  % 0-1000km
y_range = linspace(0, 1000000, nrows);  % 0-1000km

% 地形高度 (0~5000米)
Z_china = double(Z_china);
min_elev = min(Z_china(:));
max_elev = max(Z_china(:));

Z_china = (Z_china - min_elev) ./ (max_elev - min_elev) * 5000;

% 创建地形插值函数
[X_terrain, Y_terrain] = meshgrid(x_range, y_range);
terrain = @(x,y) interp2(X_terrain, Y_terrain, Z_china, x, y, 'linear', 0);

% 定义污染源参数
sourceX = [300000, 700000];   % 污染源X坐标 (m) 300km, 700km
sourceY = [300000, 700000];   % 污染源Y坐标 (m)
sourceZ = [5000, 5000];       % 污染源高度 (m)，确保高于地形最高点
Q = [1.0, 0.5];               % 排放速率 (kg/s)

% 模拟参数
T_end = 24*3600;        % 总时长24小时
dt = 300;               % 时间步长300秒 (5分钟)
numP = 5000;            % 粒子总数

% 初始化粒子
particles = zeros(numP, 3);
absorbed = false(numP, 1);  % 吸收状态标记
deposit_time = nan(numP, 1);
deposit_position = nan(numP, 3);
rebounded = false(numP, 1); % 反弹粒子标志

% 分配粒子到源
for i = 1:length(sourceX)
    idx = (1 + (i-1)*numP/length(sourceX)) : (i*numP/length(sourceX));
    particles(idx, :) = repmat([sourceX(i), sourceY(i), sourceZ(i)], length(idx), 1);
end

% 旋转风场参数
vortex_center = [500000, 500000];  % 涡旋中心500km处
vortex_strength = 5e6;             % 涡旋强度 (适配1000km区域)
max_wind_speed = 20.0;             % 最大风速 (m/s)
wind_direction = 1;                % 1=逆时针, -1=顺时针
Kx = 100; Ky = 100; Kz = 10;       % 扩散系数 (m²/s) - 适配大区域

%% 创建地形可视化
figure('Position', [100, 100, 1200, 800]);
h_terrain = surf(X_terrain, Y_terrain, Z_china, 'FaceAlpha', 0.6, 'EdgeColor', 'none');
colormap(jet);
colorbar;
hold on;

% 初始化粒子图
h_particles = scatter3(particles(:,1), particles(:,2), particles(:,3), 10, 'filled', 'MarkerFaceColor', 'b');

% 标记涡旋中心
max_terrain = max(Z_china(:));
plot3(vortex_center(1), vortex_center(2), max_terrain+1000, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');

% 图形设置
xlabel('X (m)'); ylabel('Y (m)'); zlabel('高程 (m)');
title('中国区域污染物扩散模拟 (旋转风场)');
view(35, 25);
axis tight;
grid on;
hold off;

%% 主模拟循环
for t = 0:dt:T_end
    active_idx = ~absorbed;
    numActive = sum(active_idx);

    if numActive > 0
        active_particles = particles(active_idx, :);
        
        % 计算旋转风场
        [u_local, v_local, w_local] = vortex_velocity(... 
            active_particles(:,1), active_particles(:,2), ...
            vortex_center, vortex_strength, max_wind_speed, wind_direction);

        % 湍流扩散位移
        dx = u_local * dt + sqrt(2*Kx*dt) * randn(numActive, 1);
        dy = v_local * dt + sqrt(2*Ky*dt) * randn(numActive, 1);
        dz = w_local * dt + sqrt(2*Kz*dt) * randn(numActive, 1);

        particles(active_idx, 1) = active_particles(:,1) + dx;
        particles(active_idx, 2) = active_particles(:,2) + dy;
        particles(active_idx, 3) = active_particles(:,3) + dz;

        % 检测地形碰撞
        for i = find(active_idx)'
            x = particles(i, 1);
            y = particles(i, 2);
            terrain_height = terrain(x, y);

            if particles(i, 3) < terrain_height
                if rand < 0.9 % 90%沉积
                    absorbed(i) = true;
                    deposit_time(i) = t;
                    deposit_position(i, :) = [x, y, terrain_height];
                else % 10%反弹
                    rebounded(i) = true;
                    particles(i,3) = terrain_height + 10; % 反弹10米高度
                end
            end
        end
    end

    % 每10步更新图像
    if mod(t, 10*dt) == 0
        % 更新活动粒子
        set(h_particles, 'XData', particles(~absorbed,1), ...
                         'YData', particles(~absorbed,2), ...
                         'ZData', particles(~absorbed,3));
        
        % 更新标题
        title_str = sprintf('区域污染物扩散 | 时间: %.1f 小时 | 活动粒子: %d', ...
                           t/3600, sum(~absorbed));
        title(title_str);

        % 添加沉积粒子
        hold on;
        absorbed_idx = find(absorbed);
        if ~isempty(absorbed_idx)
            scatter3(deposit_position(absorbed_idx,1), ...
                     deposit_position(absorbed_idx,2), ...
                     deposit_position(absorbed_idx,3), ...
                     20, 'r', 'filled');
        end

        % 风场可视化
        [X_quiv, Y_quiv] = meshgrid(linspace(0, 1000000, 20), linspace(0, 1000000, 20));
        [U_quiv, V_quiv, ~] = vortex_velocity(X_quiv(:), Y_quiv(:), vortex_center, vortex_strength, max_wind_speed, wind_direction);
        quiver3(X_quiv(:), Y_quiv(:), max_terrain+500*ones(size(X_quiv(:))), ...
                U_quiv, V_quiv, zeros(size(U_quiv)), 2, 'k', 'LineWidth', 1);
        hold off;
        
        drawnow;
    end
end

%% 旋转风场函数
function [u, v, w] = vortex_velocity(x, y, vortex_center, vortex_strength, max_wind_speed, wind_direction)
    dx = x - vortex_center(1);
    dy = y - vortex_center(2);
    r = sqrt(dx.^2 + dy.^2);
    r(r < 1000) = 1000;  % 避免中心处奇异值

    % 计算切向速度 (随距离衰减)
    tangential_velocity = vortex_strength ./ r;
    tangential_velocity = min(tangential_velocity, max_wind_speed);

    u = -wind_direction * tangential_velocity .* dy ./ r;
    v = wind_direction * tangential_velocity .* dx ./ r;
    
    % 垂直风速 - 在涡旋中心附近有上升气流
    w = zeros(size(x));
    center_distance = r;
    w(center_distance < 200000) = 0.5 * (1 - center_distance(center_distance < 200000)/200000);
end
