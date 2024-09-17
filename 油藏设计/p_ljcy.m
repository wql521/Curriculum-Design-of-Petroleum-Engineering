function y = p_ljcy(x)
    % 定义拟合的系数
    a = 638.6957;
    b = 1.4804e-04;
    c = -1.1640e+03;
    d = -1.8817e-04;
    
    % 计算拟合函数 f(x)
    y = a * exp(b * x) + c * exp(d * x);
end
