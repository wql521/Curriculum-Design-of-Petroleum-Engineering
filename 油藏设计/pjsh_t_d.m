function f_x = pjsh_t_d(x)
    % 定义多项式拟合的系数
    p1 = 2.4486e+04;
    p2 = 3.2097e+03;

    % 根据一次多项式公式计算 f(x)
    f_x = p1 * x + p2;
end
