function adjusted_vector = adjust_vector_length(input_vector, obj)
    % Adjust the length of 'input_vector' to match 'obj.num_t'
    % If 'input_vector' is shorter, pad with zeros.
    % If 'input_vector' is longer, truncate it.

    len = length(input_vector);
    
    if len < obj.num_t
        % Pad with zeros if the input vector is shorter
        adjusted_vector = [input_vector, zeros(1, obj.num_t - len)];
    elseif len > obj.num_t
        % Truncate if the input vector is longer
        adjusted_vector = input_vector(1:obj.num_t);
    else
        % If lengths are the same, no adjustment needed
        adjusted_vector = input_vector;
    end
end