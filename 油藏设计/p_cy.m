function y = p_cy(x)
    % 定义拟合的系数
    a = 1.9600e+03;
    b = 0.0564;
    c = 4.6576;
    d = 0.6732;
    
    % 计算拟合函数 f(x)
    y = a * exp(b * x) + c * exp(d * x);
end
