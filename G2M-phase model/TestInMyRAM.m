function  prior_val_prod = TestInMyRAM(next_para,lb,ub)
theta  = next_para;
prior_vec = zeros(length(theta), 1);
for i = 1: length(theta)
        
    prior_vec(i) = unifpdf( theta(i), lb(i), ub(i)  ) ;      
end
prior_val_prod = prod(  prior_vec);

end