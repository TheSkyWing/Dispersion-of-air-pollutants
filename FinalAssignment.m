% 定义污染源参数
sourceLon = [100, 120];   % 污染源经度
sourceLat = [30, 10];     % 污染源纬度
sourceAlt = [0, 0];       % 污染源高度 (m)
Q = [1.0, 0.5];           % 每个源的排放速率 (kg/s)

% 定义模拟时间和步长
T_end = 24*3600; dt = 3600;  % 模拟24小时，以秒为单位，步长10分钟

%% 拉格朗日粒子扩散模拟
% 粒子数目与初始化位置（所有粒子初始均位于源点处）
numP = 10000;
% 将粒子分配到不同源（均匀分配示例）
particles = zeros(numP, 3);  % 增加高度维度
for i = 1:length(sourceLon)
    idx = (1 + (i-1)*numP/length(sourceLon)) : (i*numP/length(sourceLon));
    particles(idx,:) = repmat([sourceLon(i), sourceLat(i), sourceAlt(i)], numP/length(sourceLon), 1);
end

% 简化假设: 使用常量风场 (可改为读取ERA5等数据)
u = 5; v = 0; w = 0.1;  % 水平风速和垂直风速 (m/s)
Kx = 10; Ky = 10; Kz = 1;  % 水平和垂直扩散系数 (m^2/s)

% 初始化图形窗口
figure; worldmap('World');
load coastlines
plotm(coastlat, coastlon, 'k');  % 绘制海岸线
hold on;  % 保持图形，以便动态更新
scatterm(particles(:,2), particles(:,1), 2, 'filled');  % 初始粒子位置 (lat, lon)
title('拉格朗日粒子模拟污染物分布');
hold off;

for t = 0:dt:T_end
    % 粒子随机行走：平均位移 + 湍流扩散
    particles(:,1) = particles(:,1) + u*dt/1000 + sqrt(2*Kx*dt)*randn(numP,1)/1000;  % 经度
    particles(:,2) = particles(:,2) + v*dt/1000 + sqrt(2*Ky*dt)*randn(numP,1)/1000;  % 纬度
    particles(:,3) = particles(:,3) + w*dt + sqrt(2*Kz*dt)*randn(numP,1);  % 高度

    % 更新图形
    clf;  % 清除当前图形
    worldmap('World');
    load coastlines
    plotm(coastlat, coastlon, 'k');  % 绘制海岸线
    hold on;
    scatterm(particles(:,2), particles(:,1), 2, 'filled');  % 更新粒子位置 (lat, lon)
    title(['拉格朗日粒子模拟污染物分布，时间: ', num2str(t/3600), ' 小时']);
    hold off;
    drawnow;  % 刷新图形窗口
end

%% 高斯烟羽模型模拟
% 设置点源参数
Q0 = 100.0;    % 排放强度 (kg/s)
u0 = 3.0;    % 风速 (m/s)
H = 50;      % 烟羽有效释放高度 (m)

% 扩散参数经验公式 (D类稳定度)
a_y = 0.5; b_y = 0.8;   % σ_y 参数
a_z = 0.5; b_z = 0.8;   % σ_z 参数

% 构建固定网格 (单位m)
x = 0:1000:50000;     % 下风向距离
y = -20000:1000:20000; % 横向距离
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
dt = 30; % 时间步长(秒)
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
