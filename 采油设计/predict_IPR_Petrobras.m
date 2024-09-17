function p_wf = predict_IPR_Petrobras(p_r, p_wftest, p_b, q_test, f_w,q_t)
    %计算井底流压函数
    %输入参数：
    %P_r:油藏压力(Pa)
    %p_wftest:测试点的井底流压(Pa)
    %p_b:油藏条件下的泡点压力(Pa)
    %q_test:井底流压p_wftest下的总产液量(m^3/d)
    %f_w:体积含水率
    %q_t:给定的产液量


    if p_wftest >= p_b
        J1 = q_test / (p_r - p_wftest);
    elseif p_wftest < p_b
        % 计算 A
        A = 1 - 0.2 * (p_wftest / p_b) - 0.8 * (p_wftest / p_b)^2;      
        % 计算 J1
        J1 = q_test / ((1 - f_w) * (p_r - p_b + p_b / 1.8 * A) + f_w * (p_r - p_wftest));
    end


    % 计算 qb
    qb = J1 * (p_r - p_b);
    
    % 计算 q_omax
    q_omax = qb + (J1 * p_b) / 1.8;

    disp(A)
    disp(J1)
    disp(qb)
    disp(q_omax)
 

    % 初始化 p_wf 为与 q_t 大小相同的向量
    p_wf = zeros(size(q_t));
    
    % 遍历 q_t 中的每一个元素
    for idx = 1:length(q_t)
        current_q_t = q_t(idx);
        
        if current_q_t > 0 && current_q_t <= qb
            % 当 0 < q_t <= q_b 时
            p_wf(idx) = p_r - current_q_t / J1;
        elseif qb < current_q_t && current_q_t <= q_omax
            % 当 q_b < q_t <= q_omax 时
            p_wf(idx) = f_w * (p_r - current_q_t / J1) + 0.125 * (1 - f_w) * p_b * (-1 + sqrt(81 - 80 * ((current_q_t - qb) / (q_omax - qb))));
        elseif current_q_t > q_omax
            % 当 q_omax < q_t 时
            p_wf(idx) = f_w * (p_r - q_omax / J1) + ((current_q_t - q_omax) * (8 * f_w - 9)) / J1;
        end
    end
end
