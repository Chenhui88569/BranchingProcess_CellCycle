function  output = deal_with_conv_output(conv_output, h, dx, idx, g_est2)
    output = ( conv_output(1:idx) .* dx )   - h(1:idx);
    output = sum(output.^2);
    lambda = 1e-2;
    objective_smoothness = @(g_est2) sum((g_est2(3:end) - 2 * g_est2(2:end-1) + g_est2(1:end-2)).^2);
    output  = output +     lambda*  objective_smoothness(g_est2);
end