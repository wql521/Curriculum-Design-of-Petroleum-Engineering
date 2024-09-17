function f_x = pjsh_qo(x)
    % 定义指数拟合的系数
    a = -0.0000;
    b = 0.0162;
    c = 171.9764;
    d = -0.0028;

    % 根据双指数公式计算 f(x)
    f_x = a * exp(b * x) + c * exp(d * x);
end
