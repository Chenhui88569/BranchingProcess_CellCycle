function [mat_Gm_tilde_mixture,I_mat_G_tilde_mixture]  = myfun_Gm_tilde_mixture(num_gen, ChoiceOfDist_untreated, ...
    ChoiceOfDist_treated_Dd,ChoiceOfDist_treated_D1, ChoiceOfDist_treated_U, a, b_scale,para_unknown, t_span, StepSize)

%% \tildeG matrix stores of waiting time prolonged by the checkpoint activation
num_t = length(t_span);
%d_G1 = para_unknown(end-2) ;  d_S = para_unknown(end-1) ;  d_G2M = para_unknown(end);
m12 = 1;
m1d = 0;
m23 = 1;
m2d = 0;
m31 = 2;
m3d = 0;

% %% G1
% [m1, G1_pdf, G1_cdf]  =  PdfcdfCal_treated(ChoiceOfDist_untreated,  t_span , a(1) , b_scale );
% %% G2
% [m2, S_pdf, S_cd]  =   PdfcdfCal_treated(ChoiceOfDist_untreated,  t_span , a(2) , b_scale ) ;
% %% G3
% [m3, G2M_pdf, G2M_cdf]  =   PdfcdfCal_treated(ChoiceOfDist_untreated,  t_span , a(3) ,b_scale );

nl = 2;
para_unknown = para_unknown(:);
mixture_para = num2cell( [para_unknown(20:28);  para_unknown(9)] );
[factor_G1,  factor_S,  factor_G2M, p1_G1, p2_G1, p1_S, p2_S, p1_G2M, p2_G2M, U] = mixture_para{:} ;
%% G1
alpha = [a(1)   a(1)*factor_G1];
beta = [b_scale   b_scale  ];
pesos = [p1_G1./(p1_G1 + p2_G1 ), p2_G1./(p1_G1 + p2_G1 )];
%pesos = [4/5, 1/10, 1/10];
tilde_G1_pdf =  0;
for conta = 1: length(alpha)
    %plot(X,pesos(conta)*gampdf(X,alpha(conta),beta(conta)),'m');
    tilde_G1_pdf =  tilde_G1_pdf + pesos(conta)*gampdf(t_span,alpha(conta),beta(conta));
end
tilde_G1_cdf = zeros(1,num_t);
for t = 2:num_t
        tilde_G1_cdf(t)   = trapz(0:StepSize: t_span(t) ,  tilde_G1_pdf(1:t));
end
% cov_mat_G1 = ones(1,1,2);
% cov_mat_G1(1,1,1) = 10;
% cov_mat_G1(1,1,2) = 10;
% gm_G1 = gmdistribution([ m1/1.5; m1*3   ], cov_mat_G1);
% tilde_G1_pdf= pdf(gm_G1, t_span');
% tilde_G1_cdf= cdf(gm_G1, t_span');

%% S
alpha = [a(2)   a(2)*factor_S];
beta = [b_scale   b_scale  ];
pesos = [p1_S./(p1_S + p2_S ), p2_S./(p1_S + p2_S )];
%pesos = [4/5, 1/10, 1/10];
tilde_S_pdf =  0;
for conta = 1:nl
    %plot(X,pesos(conta)*gampdf(X,alpha(conta),beta(conta)),'m');
    tilde_S_pdf =  tilde_S_pdf + pesos(conta)*gampdf(t_span,alpha(conta),beta(conta));
end
tilde_S_cdf = zeros(1,num_t);
for t = 2:num_t
        tilde_S_cdf(t)   = trapz(0:StepSize: t_span(t) ,  tilde_S_pdf(1:t));
end

% cov_mat_S = ones(1,1,2);
% cov_mat_S(1,1,1) = 10;
% cov_mat_S(1,1,2) = 10;
% gm_S = gmdistribution([ m2/1.5; m2*3   ], cov_mat_S);
% tilde_S_pdf= pdf(gm_S, t_span');
% tilde_S_cdf= cdf(gm_S, t_span');

%% G2M
alpha = [a(3)  a(3)*factor_G2M  a(3)*0.5];
beta = [b_scale   b_scale  b_scale  ];
p3_G2M = 0.3;
pesos = [p1_G2M./(p1_G2M + p2_G2M + p3_G2M ), p2_G2M./(p1_G2M + p2_G2M + p3_G2M ), p3_G2M./(p1_G2M + p2_G2M + p3_G2M )];
tilde_G2M_pdf =  0;
for conta = 1:nl
    %plot(X,pesos(conta)*gampdf(X,alpha(conta),beta(conta)),'m');
    tilde_G2M_pdf =  tilde_G2M_pdf + pesos(conta)*gampdf(t_span,alpha(conta),beta(conta));
end
tilde_G2M_cdf = zeros(1,num_t);
for t = 2:num_t
        tilde_G2M_cdf(t)   = trapz(0:StepSize: t_span(t) ,  tilde_G2M_pdf(1:t));
end

% cov_mat_G2M = ones(1,1,2);
% cov_mat_G2M(1,1,1) = 10;
% cov_mat_G2M(1,1,2) = 10;
% gm_G2M = gmdistribution([ m3/1.5; m3*3   ], cov_mat_G2M);
% tilde_G2M_pdf= pdf(gm_G2M, t_span');
% tilde_G2M_cdf= cdf(gm_G2M, t_span');

%% 
mat_G_tilde_mixture = zeros(3*num_gen+1, 3*num_gen+1 , num_t);
%single dose
for k = 1:num_gen
    mat_G_tilde_mixture( 1+ 3*(k-1),  1+ 3*(k-1) ,:   ) = tilde_G1_cdf;
    mat_G_tilde_mixture( 2+ 3*(k-1),  2+ 3*(k-1),:   )  = tilde_S_cdf; 
    mat_G_tilde_mixture( 3+ 3*(k-1),  3+ 3*(k-1),:   )  = tilde_G2M_cdf;      
end     
[~, D_pdf, D_cdf]  =   PdfcdfCal_treated('Exponential',  t_span , U ,  [] );
mat_G_tilde_mixture(3*num_gen+1 , 3*num_gen+1,:   )  =  D_cdf; 
% If the imediate new phase is G1 at time of the treatment, then the PDF of
% that phase will be that of normal G1 

%% matrix M
temp = diag( repmat([m12, m23, m31],1,num_gen )  ,1);
mat_m_tilde = [temp(1: 3*num_gen, 1: 3*num_gen),  repmat( [m1d; m2d; m3d], num_gen,1 )    ; 
    zeros(1, 3*num_gen+1 )  ];
I_mat_G_tilde_mixture = zeros(3*num_gen+1,3*num_gen+1 , num_t);
for k = 1:num_gen
        I_mat_G_tilde_mixture( 1+ 3*(k-1),  1+ 3*(k-1) ,:   ) = 1- tilde_G1_cdf;
        I_mat_G_tilde_mixture( 2+ 3*(k-1),  2+ 3*(k-1),:   )  = 1- tilde_S_cdf; 
        I_mat_G_tilde_mixture( 3+ 3*(k-1),  3+ 3*(k-1),:   )  = 1- tilde_G2M_cdf; 
end
I_mat_G_tilde_mixture( 3*num_gen+1 , 3*num_gen+1,:   )  =  1 - D_cdf;
% I_mat_3D( 3*num_gen+1 , 3*num_gen+1,: )  = ones(1, num_t);

%% matrix Gm tilde
mat_Gm_tilde_mixture = zeros(3*num_gen+1,3*num_gen+1 , num_t);
for i = 1: size(mat_m_tilde,1)
        if  i < 3*num_gen-1
            mat_Gm_tilde_mixture(i,i+1 ,: ) = reshape( mat_G_tilde_mixture( i,  i,:   ) , 1, num_t  ) * mat_m_tilde(i,i+1); %mat_m(i,j) scalar value
        end
        if  i < 3*num_gen
            mat_Gm_tilde_mixture(i, end ,: ) = reshape( mat_G_tilde_mixture( i,  i,:   ) , 1, num_t  ) * mat_m_tilde(i,end); %mat_m(i,j) scalar value
        end
end

end