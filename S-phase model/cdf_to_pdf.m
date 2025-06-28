function pdf_values = cdf_to_pdf(x, cdf_values)
    % Function to convert CDF to PDF for continuous data points

    % Check if 'x' and 'cdf_values' have the same length
    if length(x) ~= length(cdf_values)
        error('x and cdf_values must have the same length.');
    end

    % Calculate the difference in x values
    delta_x = diff(x);

    % Calculate the difference in CDF values
    delta_cdf = diff(cdf_values);

    % Calculate the PDF by taking the difference in CDF values and dividing by the difference in x values
    pdf_values = delta_cdf ./ delta_x;
end
