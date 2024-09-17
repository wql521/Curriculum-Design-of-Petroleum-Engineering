clc;
clear;
load sanyali.mat;
%绘制地层孔隙压力当量密度曲线、地层破裂压力当量密度曲线、安全地层破裂压力当量密度值曲线。
%绘制
plot(tantayali,js);
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
%设置图例
legend('地层坍塌压力当量密度曲线', '地层孔隙压力当量密度曲线','地层破裂压力当量密度曲线')