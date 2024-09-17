function y = m_bzd(x)
    % 定义参数
    a = 4.0079e+03;
    b = -3.9876;

    % 计算指数曲线拟合函数 f(x) = a * exp(b * x)
    % 输入:
    % x - 自变量

    % 输出:
    % y - 计算得到的结果

    % 计算结果
    y = a .* exp(b .* x);
end
