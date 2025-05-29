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

%% 高斯烟羽模型模拟
% 设置点源参数
Q0 = 100.0;    % 排放强度 (kg/s)
u0 = 5.0;    % 风速 (m/s)
H = 50;      % 烟羽有效释放高度 (m)

% 扩散参数经验公式 (D类稳定度)
a_y = 0.5; b_y = 0.8;   % σ_y 参数
a_z = 0.5; b_z = 0.8;   % σ_z 参数

% 构建固定网格 (单位m)
x = 0:50:100000;     % 下风向距离
y = -20000:50:20000; % 横向距离
[X, Y] = meshgrid(x, y);

% 预处理：避免除零错误
X(X == 0) = 1e-5;

% 初始化图形窗口
figure;
h_surf = surf(X/1000, Y/1000, zeros(size(X)), 'EdgeColor', 'none');
xlabel('下风向距离 (km)'); ylabel('横向距离 (km)'); 
title('高斯烟羽地面浓度分布');
colorbar;
view(2);  % 俯视图
axis tight;
clim([0, 1e-5]); % 设置合理的颜色范围

% 设置坐标轴范围保持固定
xlim([min(x)/1000, max(x)/1000]);
ylim([min(y)/1000, max(y)/1000]);

% 模拟参数
dt = 3600; % 时间步长(秒)
T_end = 24*3600; % 总时长

for t = 0:dt:T_end
    % 当前烟羽前缘位置
    x0 = u0 * t;
    
    % 计算扩散参数 (随距离变化)
    sig_y = a_y * X.^b_y;
    sig_z = a_z * X.^b_z;
    
    % 正确的高斯烟羽公式 (含地面反射)
    C = (Q0./(pi*u0*sig_y.*sig_z)) .* ...
        exp(-(Y.^2)./(2*sig_y.^2)) .* ...
        exp(-(H^2)./(2*sig_z.^2));
    
    % 屏蔽未到达区域的浓度
    C(X > x0) = 0;
    
    % 更新图形
    set(h_surf, 'ZData', C);
    set(h_surf, 'CData', C); % 同时更新颜色数据
    title(['高斯烟羽浓度 @ t = ', num2str(t/3600, '%.1f'), ' 小时']);
    drawnow;
end

%% 模拟数据输出
% 输出浓度场
[idx_i, idx_j] = meshgrid(1:size(C,1), 1:size(C,2));
lon_vals = X(:); lat_vals = Y(:);
conc_vals = C(:);
T2 = table(lon_vals, lat_vals, conc_vals, 'VariableNames', {'X_m','Y_m','Concentration'});
writetable(T2, 'gaussian_concentration.csv');
