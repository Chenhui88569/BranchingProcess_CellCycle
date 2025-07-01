function random_numbers= generateRandomNumbersAndHistogram(obj,cdf,num_samples)
    % Calculate the cumulative distribution function (CDF)
    x_values = obj.t_span;

    % Generate random numbers using inverse transform sampling
    random_numbers = zeros(1, num_samples);
    for i = 1:num_samples
        % Generate a random number between 0 and 1
        u = rand();

        % Find the index in the CDF where u falls
        index = find(cdf >= u, 1);

        % Map the random number to the corresponding value from the PDF
        random_numbers(i) = x_values(index);
    end


end

