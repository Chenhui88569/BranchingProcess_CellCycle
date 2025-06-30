function cdf = pdf_to_cdf(pdf, dx)
    % Function to convert PDF to CDF
    % pdf - input PDF values (array)
    % dx  - spacing between PDF values (scalar)
    
    cdf = cumtrapz(pdf) * dx;
    
    % If you want to normalize the CDF to make sure it goes up to 1
    cdf = cdf / cdf(end);
end