function y = m_ljcyll(x)
    % 函数 m_hsllll 计算傅里叶级数的一次拟合值
    % 输入:
    %   x - 自变量
    % 输出:
    %   y - 因变量值

    % 拟合参数
    a0 = -1.4390e+08;
    a1 = 1.4390e+08;
    b1 = 2.1169e+05;
    w = 1.6604e-05;

    % 计算函数值
    y = a0 + a1 * cos(x * w) + b1 * sin(x * w);
end
