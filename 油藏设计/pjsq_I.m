function f_x = pjsq_I(x)
    %输入f'w
    % 定义傅里叶模型的系数
    a0 = 7.1806;
    a1 = -7.1872;
    b1 = 1.0090;
    w = 0.1588;

    % 根据傅里叶公式计算 f(x)
    f_x = a0 + a1 * cos(w * x) + b1 * sin(w * x);
end
