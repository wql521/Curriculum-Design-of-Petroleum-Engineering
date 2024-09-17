function y = m_Np(x)
    % 定义拟合的系数
    p1 = 0.0006;
    p2 = 0.0258;
    p3 = -0.7693;
    
    % 计算拟合函数 f(x)
    y = p1 * x.^2 + p2 * x + p3;
end
