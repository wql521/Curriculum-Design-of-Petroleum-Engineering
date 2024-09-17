function y = m_ql(x)
    % 定义拟合的系数
    a = 7.5512;
    b = 28.8378;
    
    % 计算拟合函数 f(x)
    y = a * log(x) + b;
end
