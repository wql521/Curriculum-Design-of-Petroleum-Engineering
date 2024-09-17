function y = p_ljcy2(x)
    % 定义拟合的系数
    p1 = 8.3966e-05;
    p2 = -0.7940;
    p3 = 2.2630e+03;
    
    % 计算拟合函数 f(x)
    y = p1 * x.^2 + p2 * x + p3;
end
