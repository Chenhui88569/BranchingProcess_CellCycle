function outputArg = myfun_ImproperIntegral_Cal(t_span,Integrand_Y)
num_t = length(t_span);
outputArg = zeros(1, num_t);
for i = 0:num_t-2
    outputArg(i+1) = trapz( t_span(1:end-i),  Integrand_Y(i+1:end)  );
end
outputArg(end) = 0;
end

