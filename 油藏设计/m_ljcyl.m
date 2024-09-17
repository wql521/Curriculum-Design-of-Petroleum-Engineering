function y = m_ljcyl(x)
    % 函数 m_ljcyl 计算双指数曲线拟合的值
    % 输入:
    %   x - 自变量
    % 输出:
    %   y - 因变量值
    
    % 拟合参数
    a = 2.5361e+03;
    b = 0.0025;
    c = -2.6014e+03;
    d = -0.0373;
    
    % 计算函数值
    y = a * exp(b * x) + c * exp(d * x);
end
