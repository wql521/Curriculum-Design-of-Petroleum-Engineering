clear;clc;

% 基础参数
%含油面积(m^2)
A = 4*1e6;
%含油厚度(m)
h=8;
%平均含油饱和度
So_average = 0.68;
%平均含水饱和度
Sw_average= 0.32;
%储层渗透率(mD)
K=19*5;
%有效渗透率(mD)
k =19*5;
%平均孔隙度
Phi = 0.26;
%地下原油粘度(mPa.s)
mu_oil = 30;
%地下水粘度(mPa.s)
mu_water = 0.6;
%注采压差(MPa)
Delta_p = 2;
%输入原油密度 (kg/m^3)
rho_o = 850;
%基本建设总投资和投资回收期内操作费用总和(元)
Money_m = 800*(rho_o/1000);
%输入油体积系数
B_o = 1.2;
%输入水体积系数
B_w = 1;

%原油价格(元/t)
money_oil = 3000;
%原油价格(元/m^3)(单位换算)
money_oil = money_oil*(rho_o/1000);

%破裂压力梯度(MPa/km)
gradient = 0.02;
%平均油层中部深度(km)
H = 1.2;


%加载数据
load data.mat
%含水饱和度
Sw = data(:,1);
%油的相对渗透率
kro = data(:,2);
%水的相对渗透率
krw = data(:,3);

fw = F_W(krw,kro,mu_oil,mu_water);
% 数值求导 (差分法)
fw_prime = diff(fw) ./ diff(Sw);
fw_prime(end+1) = 0;


%束缚水饱和度
Swc = 0.32;
%残余油饱和度对应的相对渗透率
kro_Swc = 0.65;
%油藏长度(m)
L= 1000;

%见水后油井的含水饱和度
S_we = 0.44:0.014:0.67;

%水突破的饱和度
Swf = 0.44;

% 见水相对渗透率
Swm = 0.5;
krw_Swm=interp1(Sw, krw, Swm, 'linear', 'extrap');
% 半径 (ft)
d = 1000;
% 油井半径 (ft)
rw = 0.5; 
% 外半径 (ft)
rf = 2000; 
% 皮氏系数
m = 0.5; 
% 见水后的相对渗透率
kro_Swe = zeros(size(S_we));
for i = 1:length(S_we)
    kro_Swe(i) =  interp1(Sw, kro, S_we(i), 'linear', 'extrap');
end

% 定义拟合的多项式系数
p1 = 1.0169e+03;
p2 = -1.5370e+03;
p3 = 802.3758;
p4 = -272.8625;
p5 = 83.7612;
% 定义多项式函数 f(x)
f_K_kro = @(x) p1 * x.^4 + p2 * x.^3 + p3 * x.^2 + p4 * x + p5;

K_kro_Swe = f_K_kro(kro_Swe);




% 地质储量计算（kg）
N =dzcljs(A,h,Phi,So_average,rho_o, B_o);
disp('地质储量：(kg)')
disp(N)


% 井网密度计算
%最终采收率
%单位面积储量(kg/m^2)
V = N/A;

% 合理的井网密度(口/km^2)
n_max = 0;
while n_max == 0
    n_0 = 100;
    n_1 = 0;
    for n=n_1:0.001:n_0
        % 计算 E_R 最终采收率
        ER = hsnlf(k, mu_oil, n);
        %计算新增可采储量
        
        % 计算
        delta_ER_V = kccljs(V, k, mu_oil, n);
        %disp(delta_ER_V)

        %新增可采储量的价值
        Money_new = delta_ER_V * money_oil;

        Money = Money_m*delta_ER_V +1200*(h+900);
       

        if abs(Money - Money_new) < 100
            n_max = ceil(n);
        end
    end
    n_1 = n_0 +1;
    n_0 = n_0 + 100;

end


disp('井网密度：')
disp(n_max)

%

%地层破裂压力计算(MPa)
P_H = dcplyl(H,gradient);

A1 = 1414.2;
L1 = 353.6;
%油井初始产量计算(立方米/天)
qi = yjcscljs(k, kro_Swc,A1, Delta_p, mu_oil, L1);


% 排状注水开发指标预测方法
%油井见水时间
% 计算油井见水时间
%无因次注水量
t = 1:1:100;


res_jsq  = jsq(kro, krw, mu_oil, mu_water,kro_Swc,t,qi,Phi,A1,L1,Sw);

%
%见水后开发指标

fw_prime_Swe_matrix = fw_prime_Swe(Sw, krw, kro, mu_water, mu_oil, S_we);

rec_jsh = jsh(fw_prime_Swe_matrix,Sw,kro,S_we,Swc,fw,Swf,mu_oil,mu_water,krw,kro_Swc);

%

% 面积注水开发指标预测

%见水前开发指标
% 调用函数
res_mjjsq = mjjsq(K, h, kro_Swc, krw_Swm, mu_oil, mu_water, Delta_p, d, rw, rf, m,kro,krw);


% fw_Swe 函数
fw_Swe = @(Swe_val) 0.5; 

% 调用函数
res_mjjsh=mjjsh(K_kro_Swe,k, h, krw_Swm, kro_Swe, mu_water, mu_oil, Delta_p, d, rw, m,S_we);


%% 公式函数

% 地质储量计算
function N_d = dzcljs(A, h, Phi,S_o,rho_o, B_o)
    % A - 含油面积 (km^2)
    % h - 含油厚度 (m)
    % Phi - 平均孔隙度
    % S_o - 含油饱和度（矩阵）
    % rho_o - 原油密度 (g/cm^3)
    % B_o - 油体积系数 B_o (无量纲)

    % 地质储量计算
    N_d = 100 * A * h * Phi * S_o * rho_o / B_o;
end


%根据胡斯努林法改进和完善的公式
function ER = hsnlf(K, mu_oil, n)
    % K - 储层渗透率 (mD)
    % mu_oil - 地下原油粘度 (mPa.s)
    % n - 井网密度

    k_mu_ratio = K ./ mu_oil;
    ER = (0.698 + 0.16625 .* log10(k_mu_ratio)) .* exp(-0.792 ./ n .* k_mu_ratio.^(-0.253));
end

%计算每平方公里新增可采储量
function delta_ER_V = kccljs(V, K, mu_oil, n)

    % V - 单位面积储量 (x10^4t/km^2)
    % K - 储层渗透率 (mD)
    % mu_oil - 地下原油粘度 (mPa.s)
    % n - 井网密度

    k_mu_ratio = K ./ mu_oil;
    
    % 计算指数部分
    exp_term1 = exp(-0.792 ./ (n + 1) .* k_mu_ratio.^(-0.253));
    exp_term2 = exp(-0.792 ./ n .* k_mu_ratio.^(-0.253));
    
    % 计算 ΔE_R
    delta_ER = (0.698 + 0.16625 .* log10(k_mu_ratio)) .* (exp_term1 - exp_term2);
    
    % 计算 ΔE_R * V
    delta_ER_V = V .* delta_ER;
end


%地层破裂压力计算
function PH = dcplyl(H,gradient)    
    % 计算破裂压力 P_H
    PH = gradient *(1200);
end


%油井初始产量计算
function qi = yjcscljs(k, kro_Swc, A, Delta_p, mu_oil, L)
    % 计算油井的初始产量 q_i
    qi = (k * kro_Swc * A * Delta_p) / (mu_oil * L);
end


%见水前开发指标
function res = jsq(kro, krw, mu_oil, mu_water,kro_Swc,t,qi,Phi,A,L,Sw)

    V_Df = 1:1:100;
    % 计算 f_{w}^{\prime}(S_{w_f})
    fw_prime_Swf = 1./V_Df;

    % 初始化结果数组
    I_Swf = zeros(size(fw_prime_Swf));

    % 每个 fw_prime_Swf 进行积分
    for i = 1:length(fw_prime_Swf)
        % 定义积分函数，积分从 0 到 fw_prime_Swf(i)
        I_Swf(i) = integral(@(fw) 1 ./ (interp1(Sw, kro, fw, 'linear', 'extrap') + ...
                          (mu_oil / mu_water) * interp1(Sw, krw, fw, 'linear', 'extrap')), ...
                          0, fw_prime_Swf(i));
    end

    E = kro_Swc .* I_Swf - fw_prime_Swf;

    %无因次时间
    t_D = (qi./t)./ (Phi*A*L);

    %油井见水时间
    t_Df = (1./fw_prime_Swf).^2 .* (fw_prime_Swf+E./2);

    V_t = qi .* t;

    q_t = qi;

    %一维水驱油前缘饱和度推进方程
    x = ((fw_prime_Swf ./ (Phi*A))) .* V_t;

    %无因次累计注水量
    V_D = V_t ./ (Phi*A*L);

    %无因次瞬时产液量
    q_D = q_t ./ qi;

    % 将所有结果存入结构体
    res = struct();
    res.I = I_Swf;
    res.t_Df = t_Df;           % 见水时间
    res.E = E;                 % E 参数
    res.t_D = t_D;             % 无因次时间
    res.V_t = V_t;             % 累积注水体积
    res.q_t = q_t;             % 瞬时产液量
    res.x = x;                 % 水驱油前缘推进距离
    res.V_D = V_D;             % 无因次累计注水量
    res.q_D = q_D;             % 无因次瞬时产液量
 
end

%见水后开发指标
function result = jsh(fw_prime_Swe_matrix,Sw,kro,S_we_range,Swc,fw,Swf,mu_oil,mu_water,krw,kro_Swc)
    %油井井排见水后无因此注水量为
    V_D = 1./fw_prime_Swe_matrix;

    % 初始化结果数组
    q_D_array = zeros(size(S_we_range));
    V_oD_array = zeros(size(S_we_range));
    F_Swe = zeros(size(S_we_range));
    t_De = zeros(size(S_we_range));
    
    for i = 1:length(S_we_range)
        S_we = S_we_range(i);
        
        % 计算 f_{w}^{\prime}(S_{w e})
        fw_prime_Swe = interp1(Sw, fw_prime_Swe_matrix, S_we, 'linear', 'extrap');

        % 计算 f_w(S_we)
        fw_Swe = interp1(Sw, fw, S_we, 'linear', 'extrap');
        
        % 计算 k_{ro}(S_{w e})
        k_ro_Swe = interp1(Sw, kro, S_we, 'linear', 'extrap');
        
        % 计算 I(S_{w e})
        I_Swe = integral(@(Sw_val) 1 ./ (interp1(Sw, kro, Sw_val, 'linear', 'extrap') + ...
                                          (mu_oil / mu_water) * interp1(Sw, krw, Sw_val, 'linear', 'extrap')), ...
                          0, S_we);
        
        % 计算 q_D
        q_D_array(i) = fw_prime_Swe / (k_ro_Swe * I_Swe);
        % 计算 V_oD
        V_oD_array(i) = S_we + (1 - fw_Swe) / fw_prime_Swe - Swc;

        % 计算 F(S_we) 的数值积分
        fw_prime_Swf = interp1(Sw, fw_prime_Swe_matrix, Swf, 'linear', 'extrap');
        F_Swe(i) = integral(@(fw_prime_Sw_val) ...
            (I_Swe ./ (fw_prime_Sw_val.^3)), fw_prime_Swf, fw_prime_Swe);

        % 计算无因次时间 t_De
        t_De(i) = -kro_Swc * F_Swe(i);
    end

    % 将结果存入结构体
    result = struct();
    result.q_D_array = q_D_array;     % 无因次瞬时产液量
    result.V_oD_array = V_oD_array;   % 无因次累计采出量
    result.F_Swe = F_Swe;             % F(S_we) 参数
    result.t_De = t_De;               % 无因次时间


end


function fw_prime_Swe_matrix = fw_prime_Swe(Sw, krw, kro, mu_water, mu_oil, SWE)
    % 计算 f_w 对于每个 Sw 的值
    f_w = (krw ./ mu_water) ./ (krw ./ mu_water + kro ./ mu_oil);
    
    % 使用差分方法计算 f_w 的导数
    fw_prime = diff(f_w) ./ diff(Sw);
    
    % 差分结果是 n-1 个值，需要补全为 n 个值，通常在边界点使用插值或重复
    fw_prime = [fw_prime; fw_prime(end)]; % 保持与 Sw 长度一致
    disp(fw_prime)
    
    % 初始化结果矩阵，尺寸与 Swe 一致
    fw_prime_Swe_matrix = zeros(17,1);
    
    % 遍历 S_we 的每个元素
    for i = 1:17
        % 获取当前 S_we 值
        S = SWE(i);
    
        % 使用插值法找到给定 Swe 处的导数值
        fw_prime_Swe_matrix(i) = interp1(Sw, fw_prime, S, 'linear', 'extrap');
    end

end




%见水前开发指标(面积注水)
function results = mjjsq(K, h, kro_Swc, krw_Swm, mu_oil, mu_water, Delta_p, d, rw, rf, m,kro,krw)


    V_Df = 1:1:17;
    % 初始化结果数组
    q_L_1_array = zeros(length(K), 1);
    R1_array = zeros(length(K), 1);
    R2_array = zeros(length(K), 1);
    R3_array = zeros(length(K), 1);
    q_L_2_array = zeros(length(K), 1);
    
    % 循环计算每个 K 值的 q_L_1, R1, R2, R3, q_L_2
    for i = 1:length(K)
        % 当前 K 值
        current_K = K(i);
        
        % 计算 q_L_1
        q_L_1 = (pi * current_K * h * kro_Swc * Delta_p) / (mu_oil * (log(d / rw) - 0.619));
        q_L_1_array(i) = q_L_1;
        
        % 计算 R1
        R1 = (mu_water / (2 * pi * current_K * h * krw_Swm)) * log(rf / rw);
        R1_array(i) = R1;
        
        % 计算 R2
        R2 = (mu_oil / (2 * pi * current_K * h * kro_Swc)) * log(d / rf);
        R2_array(i) = R2;
        
        % 计算 R3
        R3 = (mu_oil / (2 * pi * current_K * h * kro_Swc)) * (1 / m) * log(d / (2 * (m + 1) * rw));
        R3_array(i) = R3;
        
        % 计算分母
        term1 = (mu_water / mu_oil) * (kro_Swc / krw_Swm) * log(rf / rw);
        term2 = log(d / rf);
        term3 = (1 / m) * log(d / (2 * (m + 1) * rw));    
        denominator = mu_oil * (term1 + term2 + term3);  
        
        % 计算 q_L_2
        q_L_2 = (2 * pi * current_K * h * kro_Swc * Delta_p) / denominator;
        q_L_2_array(i) = q_L_2;
    end


    % 计算 q_o
    q_o = q_L_2_array ./ m;


      % 计算 f_{w}^{\prime}(S_{w_f})
    fw_prime_Swf = 1./V_Df;

    % 初始化结果数组
    I_Swf = zeros(size(fw_prime_Swf));

    % 每个 fw_prime_Swf 进行积分
    for i = 1:length(fw_prime_Swf)
        % 定义积分函数，积分从 0 到 fw_prime_Swf(i)
        I_Swf(i) = integral(@(fw) 1 ./ (interp1(fw_prime_Swf, kro, fw, 'linear', 'extrap') + ...
                          (mu_oil / mu_water) * interp1(fw_prime_Swf, krw, fw, 'linear', 'extrap')), ...
                          0, fw_prime_Swf(i));
    end

    E = kro_Swc .* I_Swf -fw_prime_Swf;


    A = (mu_water / mu_oil) * (kro_Swc / krw_Swm);

    B = 1;

    C = (1 / m) * log(d / (2 * (m + 1) * rw));

    D = C + B * log(d) - A * log(rw);

    t = E .* ((rf^2 / 2) * ((A - B) * (log(rf) - 0.5) + D));


    % 将结果放到结构体中
    results.q_L_1 = q_L_1_array;
    results.R1 = R1_array;
    results.R2 = R2_array;
    results.R3 = R3_array;
    results.q_L_2 = q_L_2_array;
    results.q_o = q_o;
    results.t = t;

end


%见水后开发指标(面积注水)
function results = mjjsh(K_kro_Swe,k, h, krw_Swm, kro_Swe, mu_water, mu_oil, Delta_p, d, rw, m,S_we)
    % 初始化结果数组
    num_K = length(K_kro_Swe);
    R1_array = zeros(num_K, 1);
    R2_array = zeros(num_K, 1);
    qL_array = zeros(num_K, 1);
    qpL_array = zeros(num_K, 1);
    qL_high_water_array = zeros(num_K, 1);

    % 循环计算每个 K 值的结果
    for i = 1:num_K
        % 当前 K 值
        current_K = K_kro_Swe(i);

        current_Swe = S_we(i);

        % 计算 R1
        R1 = (mu_water / (2 * pi * current_K * h * krw_Swm)) * log(d / rw);
        R1_array(i) = R1;

        % 计算 R2
        R2 = (mu_oil / (2 * pi * current_K * h * kro_Swe(i))) * (1 / m) * log(d / (2 * (m + 1) * rw));
        R2_array(i) = R2;

        % 计算 qL
        term1 = (mu_water / mu_oil) * (krw_Swm / kro_Swe(i));
        term2 = log(d / rw);
        term3 = (1 / m) * log(d / (2 * (m + 1) * rw));
        denominator = mu_oil * (term1 + term2 + term3);
        qL = (2 * pi * current_K * h * krw_Swm * Delta_p*h) / denominator;
        qL_array(i) = qL;

        % 计算 qpL
        qpL = qL / m;
        qpL_array(i) = qpL;

        % 计算 qL_high_water
        qL_high_water = (2 * pi * current_K * h * krw_Swm * Delta_p*h) / ...
            (mu_oil * (log(d / rw) + (krw_Swm / kro_Swe(i)) * (1 / m) * log(d / (2 * (m + 1) * rw))));
        qL_high_water_array(i) = qL_high_water;


    end



    % 定义多项式系数
    p1 = -225.8213;
    p2 = 615.9447;
    p3 = -635.4226;
    p4 = 298.4858;
    p5 = -57.3667;
    p6 = 2.8995;
    
    % 定义计算函数
    F_Swe = @(S_we) p1 * S_we.^5 + p2 * S_we.^4 + p3 * S_we.^3 + p4 * S_we.^2 + p5 * S_we + p6;
    
    % 计算 S_we 范围内每个值的 f_Swe
    f_Swe = F_Swe(S_we);

    q_o = (qL_high_water_array ./m ) .* (1-f_Swe);


    % 将结果存储到结构体中
    results.R1 = R1_array;
    results.R2 = R2_array;
    results.qL = qL_array;
    results.qpL = qpL_array;
    results.qL_high_water = qL_high_water_array;
    results.qo = q_o;
end



%含水率曲线
function [fw] = F_W(krw,kro,mu_oil,mu_water)
    f_w = (krw./mu_water)./((krw./mu_water)+(kro./mu_oil));
    %返回值
    fw = f_w;
end


