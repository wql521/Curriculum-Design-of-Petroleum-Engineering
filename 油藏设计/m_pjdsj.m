function y = m_pjdsj(x)
    % 函数 m_pjdsj 计算二次多项式拟合值
    % 输入:
    %   x - 自变量
    % 输出:
    %   y - 因变量值

    % 拟合参数
    p1 = -0.0003;
    p2 = 0.2140;
    p3 = -0.6403;

    % 计算函数值
    y = p1 * x.^2 + p2 * x + p3;
end
