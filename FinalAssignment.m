% 定义污染源参数
sourceLon = [100, 120];   % 污染源经度（可设置多个源）
sourceLat = [30, 10];     % 污染源纬度
Q = [1.0, 0.5];           % 每个源的排放速率 (kg/s)
% 定义模拟时间和步长
T_end = 24*3600; dt = 600;  % 模拟24小时，以秒为单位，步长10分钟

%%拉格朗日粒子扩散模拟，
% 实际应用中，可替换风速u,v为读入的网格风场，并根据高度层或地形进行改进；也可使用更精细的随机过程模拟垂直扩散和沉降等效应。
% 粒子数目与初始化位置（所有粒子初始均位于源点处）
numP = 10000;
% 将粒子分配到不同源（均匀分配示例）
particles = zeros(numP,2); 
for i=1:length(sourceLon)
    idx = (1 + (i-1)*numP/length(sourceLon)) : (i*numP/length(sourceLon));
    particles(idx,:) = repmat([sourceLon(i), sourceLat(i)], numP/length(sourceLon), 1);
end
% 简化假设: 使用常量风场 (可改为读取ERA5等数据)
u = 5; v = 0;    % 水平风速 (m/s)
K = 10;          % 粒子水平扩散系数 (m^2/s)
for t = 0:dt:T_end
    % 粒子随机行走：平均位移 + 湍流扩散
    particles(:,1) = particles(:,1) + u*dt/1000 + sqrt(2*K*dt)*randn(numP,1)/1000;
    particles(:,2) = particles(:,2) + v*dt/1000 + sqrt(2*K*dt)*randn(numP,1)/1000;
end
% 粒子模拟结束，particles(:,1)为经度，(:,2)为纬度

%%高斯烟羽模型模拟
%在二维网格上计算了高斯烟羽模型的地面浓度分布。可以绘制C的二维图像或轮廓，观察随距离减弱的规律。多源情形下，可对每个源分别计算浓度并叠加。
% 设置单个点源参数
Q0 = 1.0;    % 排放强度 (kg/s)
u0 = 5.0;    % 风速 (m/s)
H = 50;      % 烟羽有效释放高度 (m)
% 计算扩散参数（简化常数或经验公式）
sig_y = 500; sig_z = 200;  
% 构建网格（单位m）
[x,y] = meshgrid(0:1000:50000, -20000:1000:20000);
% 计算浓度分布 (假设z=0地面浓度)
C = (Q0./(2*pi*u0*sig_y*sig_z)) ...
    .* exp(-y.^2/(2*sig_y^2)) .* (exp(-(H).^2/(2*sig_z^2))+exp(-(H).^2/(2*sig_z^2)));
% C是网格点浓度 (kg/m^3)，可转为µg/m^3等

%%模拟结果可视化
%可用worldmap、plotm等函数在地图上标出粒子位置或等浓度线；也可使用mesh、surf等绘制浓度场平面图。MATLAB支持将这些图像输出为动画（利用getframe和VideoWriter）展示污染物随时间传播过程。
% 示例：绘制拉格朗日粒子终点分布（散点图）
figure; worldmap('World');
load coastlines
plotm(coastlat, coastlon, 'k');  % 绘制海岸线
scatterm(particles(:,2), particles(:,1), 2, 'filled');  % 粒子位置 (lat, lon)
title('拉格朗日粒子模拟污染物分布');
% 示例：绘制高斯模型浓度场
figure; surf(x/1000, y/1000, C, 'EdgeColor','none');
xlabel('下风向距离 (km)'); ylabel('横向距离 (km)'); zlabel('浓度');
title('高斯烟羽地面浓度分布');
colorbar; view(2);  % 俯视图

%%模拟数据输出
%使用writetable将模拟结果保存为CSV格式表格：粒子模拟保存经纬度坐标，高斯模型保存网格坐标和浓度值。可根据需要输出不同字段的CSV文件，方便后续分析和可视化。
% 假设我们要输出拉格朗日粒子最终位置的表格
T = table(particles(:,1), particles(:,2), 'VariableNames', {'Longitude','Latitude'});
writetable(T, 'particle_endpoints.csv');
% 输出高斯模型浓度场为CSV
[idx_i, idx_j] = meshgrid(1:size(C,1), 1:size(C,2));
lon_vals = x(:); lat_vals = y(:);
conc_vals = C(:);
T2 = table(lon_vals, lat_vals, conc_vals, 'VariableNames', {'X_m','Y_m','Concentration'});
writetable(T2, 'gaussian_concentration.csv');
