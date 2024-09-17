function y = p_hsl(x)
    % 定义拟合的系数
    p1 = 0.0002;
    p2 = 0.2142;
    
    % 计算拟合函数 f(x)
    y = p1 * x + p2;
end
