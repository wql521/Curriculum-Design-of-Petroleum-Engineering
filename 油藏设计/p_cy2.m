function y = p_cy2(x)
    % 定义拟合的系数
    a = 3.1182e+03;
    b = -0.2146;
    c = 2.8039e+03;
    d = -0.0502;
    
    % 计算拟合函数 f(x)
    y = a * exp(b * x) + c * exp(d * x);
end
