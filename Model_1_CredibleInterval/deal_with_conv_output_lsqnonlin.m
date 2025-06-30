function  output = deal_with_conv_output_lsqnonlin(conv_output, h, dx, idx)
    output = ( conv_output(1:idx) .* dx )   - h(1:idx);
end