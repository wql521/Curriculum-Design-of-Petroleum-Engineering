function y = m_hslll(x)
    % 函数 m_hslll 计算 log2(x) 曲线拟合的值
    % 输入:
    %   x - 自变量
    % 输出:
    %   y - 因变量值
    
    % 拟合参数
    a = 36.5530;
    b = -137.6095;
    
    % 计算函数值
    y = a * log2(x) + b;
end
