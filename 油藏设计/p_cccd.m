function y = p_cccd(x)
    % 定义拟合的系数
    a = 24.6699;
    b = -187.8388;
    
    % 计算拟合函数 f(x)
    % 确保 x > 0，因为对数函数在 x <= 0 时无定义
    if any(x <= 0)
        error('输入 x 必须大于 0');
    end
    
    y = a * log(x) + b;
end
