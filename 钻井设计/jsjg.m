clc;clear;
%%设计井身机构
%%地层三压力数据导入
load sanyali.mat
S_f=0.03; %地层破裂压力当量密度安全允许值
anquanpolieyali=polieyali-S_f;
%绘制地层孔隙压力当量密度曲线、地层破裂压力当量密度曲线、安全地层破裂压力当量密度值曲线。
%绘制
plot(anquanpolieyali,js,'-o');
hold on;
plot(kongxiyali,js);
hold on;
plot(polieyali,js);
% 设置轴标签
xlabel('当量密度g/cm^3');
ylabel('井深');
% 显示网格线
grid on;
set(gca, 'YDir', 'reverse');

%最大孔隙压力当量密度
kongxiyali_max = max(kongxiyali);
kongxiyali_max_js = 3558;

%钻井液密度
deltp = 0.1; %钻井液密度附加值
zjmd=max(kongxiyali) + deltp;
tantayali_max = max(tantayali);
%最大钻井液密度
zjmd_max = zjmd;
if tantayali_max > zjmd
    zjmd_max = tantayali_max;
end
%最大井内压力当量密度（正常作业时）
S_g= 0.035;%激动压力当量密度，
jnylmd_max = zjmd_max+ S_g;
%最大井内压力当量密度（放生溢流时）
S_k= 0.08;%溢流允许值
jnylmd_max_yl = zjmd_max + (kongxiyali_max_js/15)*S_k;

%%
% 确定中间套管下入深度初选点D21
S_b = 0.04; %抽吸压力系数
desired_depth = 0;
% 逐步计算等效密度并验证条件
for J = 1:length(js)
    % 计算
    rec = kongxiyali_max + S_b + S_f + (kongxiyali_max_js / js(J)) * S_k;
    % 验证条件
    if polieyali(J) >= rec
        zhongjian = [];
        for p = 1:length(polieyali)
            if polieyali(p) == polieyali(J)
                zhongjian(end+1) = js(p);
            end
        end
        % disp(zhongjian)
        % disp(polieyali(J))
        desired_depth = max(zhongjian);
        break; % 找到第一个符合条件的井深即可退出循环
    end
end
D21 = desired_depth;
%校核中间套管
kongxiyali_D21 = kongxiyali(find(js == D21));
kongxiyali_min = min(kongxiyali);
D_n= 3870;
delt_p = 0.00981*D_n*(kongxiyali_D21+S_b-kongxiyali_min);
if delt_p > 12
    disp("重算所允许的最大地层压力当量密度");
end

%%
%确定尾管下入深度D3
polieyali_D21 = polieyali(find(js == D21));
% 逐步计算等效密度并验证条件
desired_depth_3  =0;
o = find(js == D21);
for j = o:length(js)
    % 计算
    rec2 = polieyali_D21 - S_b - S_f - (js(j) / D21) * S_k;
    % disp(rec2)
    % 验证条件
    if polieyali(j) < rec2
        desired_depth_3 = js(j);
        break; % 找到第一个符合条件的井深即可退出循环
    end
end


%%
%计算表层套管
%中间套管以上的最大孔隙压力当量密度
kongxiyali_max_D21 = max(kongxiyali(1:o));

desired_depth1 = 0;
% 逐步计算等效密度并验证条件
for i = 1:o
    % 计算
    re = kongxiyali_max_D21 + S_b + S_f + (kongxiyali_max_js / js(i)) * S_k;
    % 验证条件
    if polieyali(i) >= re
        desired_depth1 = js(i);
        break; % 找到第一个符合条件的井深即可退出循环
    end
end
D11 = desired_depth1;

disp(['表层套管的深度为：', num2str(D11), ' 米']);
disp(['中间套管的深度为：', num2str(D21), ' 米']);

%%
%钻头尺寸和套管尺寸选择（毫米）
%一开
zuantou_1 = 444.5;
taoguan_1 = 339.7;
%二开
zuantou_2 = 311.2;
taoguan_2 = 244.5;
%三开
zuantou_3 = 200.0;
taoguan_3 = 139.7;

%%
%套管柱的设计

%生产套管
%确定第一段套管的钢级和壁厚
%最大外挤载荷
P_oc = 0.00981*zjmd*(max(js));
%允许抗外挤强度
S_d = 1.125; %抗外挤安全系数
P_c = P_oc*S_d;
%钢级N-80 
%屈服强度
p_n_80 = 552;
%壁厚
%t_1 = P_c * taoguan_3 / (2 * p_n_80 - P_c);
t_1 = 9.17;
%内径
r_1 = 121.4;
%抗挤强度
p_c_1 = 60.9;

%第二段套管的钢级和厚度
%钢级N-80 
%壁厚
t_2 = 7.72;
%内径
r_2 = 124.3;
%抗挤强度
p_c_2 = 43.3;
%确定第二段套管下入深度
d2 = p_c_2/(zjmd*0.00981*S_d);
%第一段套管使用长度
L1 = max(js)-d2;
%第一段套管的根数
n1 = ceil(L1/9.1);
L1_new = n1*9.1;
%第二段套管实际下的井深
d2_new = max(js)-L1_new;
%套管所受拉力
F_m_1 = 29.763*9.8*(1-(zjmd/7.85))*L1_new*0.001;

%校核第二段套管及其使用长度
%套管根数
n2 = ceil(d2_new/9.1);
%第二段套管的实际使用长度
L2_new = n2*9.1;
%第二段套管所受的最大拉应力为
F_m_2 = 25.298*9.8*(1-(zjmd/7.85))*L2_new*10^(-3);

%技术套管
%第一段套管的钢级和壁厚
%钢级N-80 
%壁厚
t_3 = 11.05;
%内径
r_3 = 222.4;
%抗挤强度
p_c_3 = 26.269;

%第二段套管的钢级和壁厚
%钢级N-80 
%壁厚
t_4 = 10.03;
%内径
r_4 = 224.4;
%抗挤强度
p_c_4 = 21.304;
%第二段套管下入深度
d4 = p_c_4/(zjmd*0.00981*S_d);
%第一段使用长度
L3=1500-d4;
%第一段套管根数
n3 =ceil(L3/9.1);
%第一段套管使用长度
L3_new = n3*9.1;
%第二段套管实际下的井深
d4_new = 1500-L3_new;

%第二段套管根数
n4 = ceil(d4/9.1);
%第二段套管的实际使用长度
L4_new = n4*9.1;
%第二段套管所受的最大拉应力为
F_m_4 = 59.52*9.8*(1-(zjmd/7.85))*L4_new*10^(-3);

%表层套管柱的设计
%钢级N-80 
%壁厚
t_5 = 13.06;
%内径
r_5 = 313.6;
%抗挤强度
p_c_5 = 18.3619;
%表层套管的根数
n5 = floor(D11/9.1);
%实际使用的长度
L5_new = n5*9.1;
%表层套管受到的最大拉应力
F_m_5 = 107.14*9.8*(1-(zjmd/7.85))*L5_new*10^(-3);

%%
%钻具设计

% 浮力系数
K_B = 0.8462;
% 安全系数
S_N = 1.25;
% 抗拉强度安全系数
S_t = 1.3;
% 抗挤强度安全系数
S_c = 1.125;
% 屈服强度与拉伸应力的比值
ql = 1.52;
% 最大操作压力 (MPa)
MOP = 400;

% 钻铤和钻杆参数
% 表层套管段
r1_in = 71.4; % 内径 (mm)
r1_on = 203.2; % 外径 (mm)
qc_1 = 223.5; % 每米钻铤在空气中的重力 (kg/m)
W_max_1 = 40; % 最大钻压 (吨)

R1_in = 151.54; % 钻杆内径 (mm)
R1_on = 168.8;  % 钻杆外径 (mm)
qp_1 = 37.54;   % 单位长度钻杆在空气中的重力 (kg/m)
F_y_1 = 2177;   % 最小屈服强度下的抗拉力 (kN)
p_1 = 33;       % 抗挤强度 (MPa)
h1 = 8.38;      % 壁厚 (mm)

% 中间套管段
r2_in = 71.4;
r2_on = 203.2;
qc_2 = 223.5;
W_max_2 = 60;

R2_in = 151.54;
R2_on = 168.8;
qp_2 = 37.54;
F_y_2 = 2177;
h2 = 8.38;

% 生产套管段
r3_in = 71.4;
r3_on = 203.2;
qc_3 = 223.5;
W_max_3 = 80;

R3_in = 108.62;
R3_on = 127;
qp_3 = 29.05;
F_y_3 = 1760;
h3 = 9.19;

% 套管段深度
D11 = 714; % 表层套管段深度（米）
D21 = 1500; % 中间套管段深度（米）
D31 = 4025; % 生产套管段深度（米）

% 第一开表层套管段
Lc_1 = ceil(((S_N * W_max_1 * 1000) / (qc_1 * 9.81 * K_B))/9.1)*9.1; % 计算钻铤长度
Fa_1 = min([(0.9 * F_y_1) / S_t, (0.9 * F_y_1) / ql, 0.9 * F_y_1 - MOP]); % 计算抗拉强度校核
Lp_1 = D11 - Lc_1; % 计算钻杆长度，确保总长度为段深度

if Lp_1 < 0 % 检查计算结果合理性
    Lc_1 = D11;
    Lp_1 = 0;
end

fprintf('第一开表层套管段：\n');
fprintf('钻铤长度：%.2f 米\n', Lc_1);
fprintf('钻杆长度：%.2f 米\n', Lp_1);

% 第二开中间套管段
Lc_2 = ceil(((S_N * W_max_2 * 1000) / (qc_2 * 9.81 * K_B))/9.1)*9.1; % 计算钻铤长度
Fa_2 = min([(0.9 * F_y_2) / S_t, (0.9 * F_y_2) / ql, 0.9 * F_y_2 - MOP]); % 计算抗拉强度校核
Lp_2 = (D21 - D11) - Lc_2; % 计算钻杆长度，确保总长度为段深度

if Lp_2 < 0 % 检查计算结果合理性
    Lc_2 = D21 - D11;
    Lp_2 = 0;
end

fprintf('第二开中间套管段：\n');
fprintf('钻铤长度：%.2f 米\n', Lc_2);
fprintf('钻杆长度：%.2f 米\n', Lp_2);

% 第三开生产套管段
Lc_3 = ceil(((S_N * W_max_3 * 1000) / (qc_3 * 9.81 * K_B))/9.1)*9.1; % 计算钻铤长度
Fa_3 = min([(0.9 * F_y_3) / S_t, (0.9 * F_y_3) / ql, 0.9 * F_y_3 - MOP]); % 计算抗拉强度校核
Lp_3 = (D31 - D21) - Lc_3; % 计算钻杆长度，确保总长度为段深度

if Lp_3 < 0 % 检查计算结果合理性
    Lc_3 = D31 - D21;
    Lp_3 = 0;
end

fprintf('第三开生产套管段：\n');
fprintf('钻铤长度：%.2f 米\n', Lc_3);
fprintf('钻杆长度：%.2f 米\n', Lp_3);

%%
%钻井液设计
%井筒内钻井液体积
V = (pi/4)*((0.4445)^2)*710+(pi/4)*((0.3112)^2)*(1500-710)+(pi/4)*((0.2)^2)*(4025-1500);

