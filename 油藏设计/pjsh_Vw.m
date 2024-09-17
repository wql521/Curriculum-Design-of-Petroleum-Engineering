function f_x = pjsh_Vw(x)
    % 定义多项式拟合的系数
    p1 = 0.9333;
    p2 = -44.5925;

    % 根据一次多项式公式计算 f(x)
    f_x = p1 * x + p2;
end
