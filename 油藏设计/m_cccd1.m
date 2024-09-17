function y = m_cccd1(x)
    % 定义拟合的系数
    p1 = 0.0005;
    p2 = 0.0189;
    p3 = -0.5670;
    
    % 计算拟合函数 f(x)
    y = p1 * x.^2 + p2 * x + p3;
end
