%% 三维山地地形模型 - 旋转风场
% 定义地形参数
A1 = 300; f1 = 0.0001;   % 第一个正弦波：振幅300m，频率0.0001 (周期10km)
A2 = 200; f2 = 0.0003;   % 第二个正弦波：振幅200m，频率0.0003 (周期3.33km)
terrain = @(x,y) A1 * sin(2*pi*f1*x) .* sin(2*pi*f1*y) + ...
                 A2 * sin(2*pi*f2*x) .* sin(2*pi*f2*y);

% 定义污染源参数 (直角坐标系)
sourceX = [25000, 75000];   % 污染源X坐标 (m)
sourceY = [25000, 75000];   % 污染源Y坐标 (m)
sourceZ = [500, 500];       % 污染源高度 (m)，初始高度高于地形
Q = [1.0, 0.5];             % 排放速率 (kg/s)

% 模拟参数
T_end = 24*3600;        % 总时长24小时
dt = 60;                % 时间步长60秒
numP = 10000;           % 粒子总数

% 初始化粒子
particles = zeros(numP, 3);
absorbed = false(numP, 1);  % 吸收状态标记

% 记录沉积时间和沉积位置
deposit_time = nan(numP, 1);
deposit_position = nan(numP, 3);

% 添加反弹粒子标志
rebounded = false(numP, 1);

% 分配粒子到源
for i = 1:length(sourceX)
    idx = (1 + (i-1)*numP/length(sourceX)) : (i*numP/length(sourceX));
    particles(idx, :) = repmat([sourceX(i), sourceY(i), sourceZ(i)], ...
                              length(idx), 1);
end

% 旋转风场参数
vortex_center = [50000, 50000]; % 涡旋中心位置 (x,y)
vortex_strength = 10000;        % 涡旋强度 (控制旋转速度)
max_wind_speed = 50.0;          % 最大风速 (m/s)
wind_direction = 1;             % 1=逆时针, -1=顺时针
Kx = 10; Ky = 10; Kz = 1;        % 扩散系数 (m²/s)

%% 创建地形可视化
% 生成地形网格
x_range = 0:200:100000;
y_range = 0:200:100000;
[X_terrain, Y_terrain] = meshgrid(x_range, y_range);
Z_terrain = terrain(X_terrain, Y_terrain);

% 创建3D图形
figure('Position', [100, 100, 1200, 800]);
h_terrain = surf(X_terrain, Y_terrain, Z_terrain, 'FaceAlpha', 0.6, 'EdgeColor', 'none');
colormap(jet);
hold on;

% 初始化粒子图
h_particles = scatter3(particles(:,1), particles(:,2), particles(:,3), 10, 'filled', 'MarkerFaceColor', 'b');

% 标记涡旋中心
plot3(vortex_center(1), vortex_center(2), max(Z_terrain(:))+100, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');

% 图形设置
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
title('山地地形污染物扩散模拟 (旋转风场)');
view(35, 25);  % 设置视角
axis tight;    % 自动调整轴范围
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

        % 更新粒子位置
        particles(active_idx, 1) = active_particles(:,1) + dx;
        particles(active_idx, 2) = active_particles(:,2) + dy;
        particles(active_idx, 3) = active_particles(:,3) + dz;

        % 检测地形碰撞
        for i = find(active_idx)'
            x = particles(i, 1);
            y = particles(i, 2);
            terrain_height = terrain(x, y);

            if particles(i, 3) < terrain_height
                % 新增机制：90% 概率沉积，10% 概率反弹
                if rand < 0.9
                    absorbed(i) = true;
                    deposit_time(i) = t;
                    deposit_position(i, :) = [x, y, terrain_height];
                else
                    rebounded(i) = true;
                    particles(i,3) = terrain_height + 0.1;  % 反弹微高度
                end
            end
        end
    end

    % 每10步更新图像
    if mod(t, 10*dt) == 0
        % 更新粒子位置数据
        set(h_particles, 'XData', particles(~absorbed,1), ...
                         'YData', particles(~absorbed,2), ...
                         'ZData', particles(~absorbed,3));

        % 更新标题
        title_str = sprintf('山地地形污染物扩散 | 时间: %.1f 小时 | 活动粒子: %d | 涡旋强度: %.1e', ...
                           t/3600, sum(~absorbed), vortex_strength);
        title(title_str);

        % 添加被吸收粒子的标记 (红色)
        hold on;
        absorbed_idx = find(absorbed);
        if ~isempty(absorbed_idx)
            terrain_heights = arrayfun(@(i) terrain(particles(i,1), particles(i,2)), absorbed_idx);
            scatter3(particles(absorbed_idx,1), ...
                     particles(absorbed_idx,2), ...
                     terrain_heights, ...
                     20, 'r', 'filled');
        end

        % 添加风场可视化 (箭头)
        [X_quiv, Y_quiv] = meshgrid(linspace(0, 100000, 15), linspace(0, 100000, 15));
        [U_quiv, V_quiv, ~] = vortex_velocity(X_quiv(:), Y_quiv(:), vortex_center, vortex_strength, max_wind_speed, wind_direction);
        quiver3(X_quiv(:), Y_quiv(:), max(Z_terrain(:))+50*ones(size(X_quiv(:))), ...
                U_quiv, V_quiv, zeros(size(U_quiv)), 2, 'k', 'LineWidth', 1);
        hold off;

        drawnow;
    end
end

%% 旋转风场函数
% 该函数计算给定位置的旋转风速分量
function [u, v, w] = vortex_velocity(x, y, vortex_center, vortex_strength, max_wind_speed, wind_direction)
    dx = x - vortex_center(1);
    dy = y - vortex_center(2);
    r = sqrt(dx.^2 + dy.^2);
    r(r < 1) = 1;

    tangential_velocity = vortex_strength ./ r;
    tangential_velocity = min(tangential_velocity, max_wind_speed);

    u = -wind_direction * tangential_velocity .* dy ./ r;
    v = wind_direction * tangential_velocity .* dx ./ r;

    % 垂直风速 - 在涡旋中心附近有上升气流
    w = zeros(size(x));
    center_distance = sqrt((x - vortex_center(1)).^2 + (y - vortex_center(2)).^2);
    w(center_distance < 10000) = 0.1 * (1 - center_distance(center_distance < 10000)/10000);
end
