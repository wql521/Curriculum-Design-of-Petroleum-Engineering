function y = p_ljcy1(x)
    % 定义拟合的系数
    a0 = -5.3489e+09;
    a1 = 5.3489e+09;
    b1 = 4.9565e+06;
    w = 1.7956e-07;
    
    % 计算拟合函数 f(x)
    y = a0 + a1 * cos(x * w) + b1 * sin(x * w);
end
