function y = pq_1(x)
    % 定义拟合的系数
    a1 = 0;
    b1 = 38.0773;
    c1 = 0.1968;
    
    a2 = 2.2976e+04;
    b2 = 42.8437;
    c2 = 26.8706;
    
    a3 = 0;
    b3 = 38.9156;
    c3 = 0.3889;
    
    a4 = 2.9950e+03;
    b4 = 18.3052;
    c4 = 8.9054;
    
    % 计算拟合函数 f(x)
    y = a1 .* exp(-((x - b1) / c1).^2) + ...
        a2 .* exp(-((x - b2) / c2).^2) + ...
        a3 .* exp(-((x - b3) / c3).^2) + ...
        a4 .* exp(-((x - b4) / c4).^2);
end
