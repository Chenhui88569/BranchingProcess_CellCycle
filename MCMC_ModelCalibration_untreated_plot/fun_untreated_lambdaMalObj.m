function p_Mal  = fun_untreated_lambdaMalObj(theta)

a = theta(1:3);
b = theta(4);
a_all = sum(a);

%skewness_coeff =   2/sqrt(a_all)     ;
lambda_fun  = @(m, mu, CV) log(m)/mu  + (CV^2)* ( log(m)^2)/(2*mu);  %+ CV^2*(3*CV^2 -  skewness_coeff *CV )*log(m)^2/(6*mu)  ;
variance = a_all/b^2;
m = 2;
mu = a_all/b;
CV = sqrt(variance)/mu;
lambda_Mal =  lambda_fun(m,mu, CV);

p1_Mal  = 2*(1- (b^a(1)/(b+lambda_Mal)^a(1) ) );
p2_Mal = 2*( ( b/(b+lambda_Mal) )^a(1) - ( ( b/(b+lambda_Mal) )^a(1) )  * ( b/(b+lambda_Mal) )^a(2)   ); % needs modification 
p3_Mal = 2*( ( ( b/(b+lambda_Mal) )^a(1) )  * ( b/(b+lambda_Mal) )^a(2)  -( b/(b+lambda_Mal) )^a(1) * ( b/(b+lambda_Mal) )^a(2)  * ( b/(b+lambda_Mal) )^a(3)); % needs modification
p_Mal = [ p1_Mal  ; p2_Mal;p3_Mal   ];

end