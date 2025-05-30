%% 三维山地地形模型
% 定义地形参数
A1 = 300; f1 = 0.0001;   % 第一个正弦波：振幅500m，频率0.0001 (周期10km)
A2 = 200; f2 = 0.0003;   % 第二个正弦波：振幅200m，频率0.0003 (周期3.33km)
terrain = @(x,y) A1 * sin(2*pi*f1*x) .* sin(2*pi*f1*y) + ...
                 A2 * sin(2*pi*f2*x) .* sin(2*pi*f2*y);

% 定义污染源参数 (直角坐标系)
sourceX = [0, 10000];   % 污染源X坐标 (m)
sourceY = [0, 0];       % 污染源Y坐标 (m)
sourceZ = [500, 500];   % 污染源高度 (m)，初始高度高于地形
Q = [1.0, 0.5];         % 排放速率 (kg/s)

% 模拟参数
T_end = 24*3600;        % 总时长24小时
dt = 60;                % 时间步长60秒
numP = 10000;           % 粒子总数

% 初始化粒子
particles = zeros(numP, 3);
absorbed = false(numP, 1);  % 吸收状态标记

% 分配粒子到源
for i = 1:length(sourceX)
    idx = (1 + (i-1)*numP/length(sourceX)) : (i*numP/length(sourceX));
    particles(idx, :) = repmat([sourceX(i), sourceY(i), sourceZ(i)], ...
                              length(idx), 1);
end

% 风场参数 - 对数风廓线
u_ref = 0.3;            % 参考高度处的风速 (m/s)
z_ref = 100;            % 参考高度 (m)
z0 = 0.1;               % 地表粗糙度长度 (m)
k = 0.4;                % 卡曼常数
Kx = 10; Ky = 10; Kz = 1;  % 扩散系数 (m²/s)

%% 创建地形可视化
% 生成地形网格
x_range = -100000:200:100000;
y_range = -100000:200:100000;
[X_terrain, Y_terrain] = meshgrid(x_range, y_range);
Z_terrain = terrain(X_terrain, Y_terrain);

% 创建3D图形
figure;
h_terrain = surf(X_terrain, Y_terrain, Z_terrain, 'FaceAlpha', 0.6, 'EdgeColor', 'none');
colormap(jet);
hold on;

% 初始化粒子图
h_particles = scatter3(particles(:,1), particles(:,2), particles(:,3), 10, 'filled', 'MarkerFaceColor', 'b');
hold off;

% 图形设置
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
title('山地地形污染物扩散模拟 (对数风廓线)');
colorbar;
view(35, 25);  % 设置视角
axis tight;
grid on;

%% 主模拟循环
for t = 0:dt:T_end
    % 只更新未被吸收的粒子
    active_idx = ~absorbed;
    numActive = sum(active_idx);
    
    if numActive > 0
        % 获取活动粒子位置
        active_particles = particles(active_idx, :);
        
        % 计算活动粒子对应的地形高度
        active_terrain_heights = terrain(active_particles(:,1), active_particles(:,2));
        
        % 计算离地高度 (确保不小于0.1m)
        h_above_ground = max(active_particles(:,3) - active_terrain_heights, 0.1);
        
        % 对数风廓线计算风速
        u_local = (u_ref * log(max(h_above_ground, z0+0.01)/z0)) / log(z_ref/z0);
        v_local = zeros(numActive, 1);
        w_local = zeros(numActive, 1);
        
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
            
            % 计算当前位置地形高度
            terrain_height = terrain(x, y);
            
            % 检查是否被地形吸收
            if particles(i, 3) < terrain_height
                absorbed(i) = true;
            end
        end
    end
    
    % 更新图形 (每10步更新一次以提高性能)
    if mod(t, 10*dt) == 0
        % 更新粒子位置数据
        set(h_particles, 'XData', particles(~absorbed, 1), ...
                         'YData', particles(~absorbed, 2), ...
                         'ZData', particles(~absorbed, 3));
        
        % 更新标题
        title(['山地地形污染物扩散 | 时间: ', num2str(t/3600, '%.1f'), ' 小时 | 活动粒子: ', num2str(sum(~absorbed))]);
        
        % 添加被吸收粒子的标记 (红色)
        hold on;
        absorbed_idx = find(absorbed);
        if ~isempty(absorbed_idx)
            % 计算被吸收粒子位置的地形高度
            terrain_heights = arrayfun(@(i) terrain(particles(i,1), particles(i,2)), absorbed_idx);
            
            scatter3(particles(absorbed_idx,1), ...
                    particles(absorbed_idx,2), ...
                    terrain_heights, ...  % 使用计算出的地形高度
                    20, 'r', 'filled');
        end
        hold off;
        
        drawnow;
    end
end
