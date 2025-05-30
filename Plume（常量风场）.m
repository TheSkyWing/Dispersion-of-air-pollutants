%% 三维高斯烟羽模型模拟（含烟气抬升）
% 设置点源参数
Q0 = 20.0;      % 排放强度 (kg/s)
u0 = 5.0;        % 风速 (m/s)
T_g = 423;       % 烟气温度 (K)
T_a = 293;       % 环境温度 (K)
d = 5.0;         % 烟囱出口直径 (m)
v_s = 5.0;      % 烟气出口速度 (m/s)
H_stack = 50;   % 烟囱物理高度 (m)

% 计算浮力通量 (Briggs公式)
g = 9.81;        % 重力加速度 (m/s²)
F_b = g * v_s * (d/2)^2 * (T_g - T_a)/T_g;

% 扩散参数经验公式 (D类稳定度)
a_y = 0.5; b_y = 0.8;   % σ_y 参数
a_z = 0.5; b_z = 0.8;   % σ_z 参数

% 构建三维网格 (单位m)
x = 0:10:10000;      % 下风向距离
y = -1000:10:1000;   % 横向距离
z = 0:10:1000;        % 垂直高度
[X, Y, Z] = meshgrid(x, y, z);

% 预处理：避免除零错误
X(X == 0) = 1e-5;

% 初始化图形窗口
figure;
set(gcf, 'Position', [100, 100, 1200, 800]);

% 子图1: 水平切片 (z=0 地面)
subplot(2,2,1);
h_surf_ground = surf(x/1000, y/1000, zeros(length(y), length(x)), 'EdgeColor', 'none');
xlabel('下风向距离 (km)'); ylabel('横向距离 (km)'); 
title('地面浓度 (z=0)');
colorbar;
view(2);
axis tight;
clim([0, 1e-5]);

% 子图2: 垂直剖面 (y=0 中心面)
subplot(2,2,2);
h_surf_profile = surf(x/1000, z, zeros(length(z), length(x)), 'EdgeColor', 'none');
xlabel('下风向距离 (km)'); ylabel('高度 (m)'); 
title('垂直剖面浓度 (y=0)');
colorbar;
view(2);
axis tight;
clim([0, 1e-5]);

% 子图3: 烟气抬升高度
subplot(2,2,3:4);
h_plume = plot(0, H_stack, 'r-', 'LineWidth', 2);
xlabel('下风向距离 (km)'); ylabel('有效源高 (m)');
title('烟气抬升轨迹');
xlim([min(x)/1000, max(x)/1000]);
ylim([H_stack-50, H_stack+450]);
grid on;

% 模拟参数
dt = 3600; % 时间步长(秒)
T_end = 6*3600; % 总时长缩短为6小时

% 预分配抬升高度数组
x_vals = x;
H_eff = zeros(size(x_vals));

for t = 0:dt:T_end
    % 当前烟羽前缘位置
    x0 = u0 * t;
    
    % 计算烟气抬升高度 (Briggs公式)
    for i = 1:length(x_vals)
        x_val = x_vals(i);
        
        % 抬升阶段 (x < 3.5x*)
        if x_val <= 0
            H_eff(i) = H_stack;
        else
            % 稳定度参数 (D类)
            s = 0.02; % 稳定度参数 (s⁻²)
            
            % 浮力主导抬升
            if F_b > 55
                deltaH = 38.7 * F_b^(3/5) / u0 * (s^(-1/5));
                x_f = 119 * F_b^(2/5); % 抬升完成距离
            else
                deltaH = 21.4 * F_b^(3/4) / u0;
                x_f = 49 * F_b^(5/8); % 抬升完成距离
            end
            
            % 抬升过程
            if x_val < x_f
                H_eff(i) = H_stack + deltaH * (x_val/x_f)^(1/3);
            else
                H_eff(i) = H_stack + deltaH;
            end
        end
    end
    
    % 创建全网格的有效高度数组
    H_eff_grid = interp1(x_vals, H_eff, X, 'linear', 'extrap');
    
    % 计算扩散参数 (随距离变化)
    sig_y = a_y * X.^b_y;
    sig_z = a_z * X.^b_z;
    
    % 三维高斯烟羽公式 (含地面反射和抬升高度)
    term1 = exp(-(Y.^2)./(2*sig_y.^2));
    term2 = exp(-(Z - H_eff_grid).^2./(2*sig_z.^2)) + ...
            exp(-(Z + H_eff_grid).^2./(2*sig_z.^2));
    C = Q0./(2*pi*u0.*sig_y.*sig_z) .* term1 .* term2;
    
    % 屏蔽未到达区域的浓度
    C(X > x0) = 0;
    
    % 提取两个特征剖面
    % 1. 地面浓度 (z=0)
    C_ground = squeeze(C(:, :, 1));  % z的第一层对应高度0
    
    % 2. 垂直剖面 (y=0)
    [~, y0_idx] = min(abs(y));  % 找到y=0的索引
    C_profile = squeeze(C(y0_idx, :, :))';
    
    % 更新图形
    subplot(2,2,1);
    set(h_surf_ground, 'ZData', C_ground, 'CData', C_ground);
    title(['地面浓度 @ t = ', num2str(t/3600, '%.1f'), ' 小时']);
    
    subplot(2,2,2);
    set(h_surf_profile, 'ZData', C_profile, 'CData', C_profile);
    title(['垂直剖面 @ y=0 @ t = ', num2str(t/3600, '%.1f'), ' 小时']);
    
    subplot(2,2,3:4);
    set(h_plume, 'XData', x_vals/1000, 'YData', H_eff);
    title(['烟气抬升轨迹 @ t = ', num2str(t/3600, '%.1f'), ' 小时 | ΔH = ', ...
          num2str(max(H_eff)-H_stack, '%.1f'), ' m']);
    drawnow;
end
