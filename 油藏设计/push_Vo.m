function f_x = push_Vo(x)
    % 定义对数拟合的系数
    a = 15.9128;
    b = -20.0000;

    % 根据对数公式计算 f(x)
    f_x = a * log(x) + b;
end
