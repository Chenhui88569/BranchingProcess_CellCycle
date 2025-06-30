function [mat_Gm_tilde,I_mat_G_tilde]  = myfun_Gm_tilde(num_gen, ChoiceOfDist_untreated, ...
    ChoiceOfDist_treated_Dd,ChoiceOfDist_treated_D1, ChoiceOfDist_treated_U, a, b_scale, para_unknown, t_span, StepSize)
%% \tildeG matrix stores of waiting time prolonged by the checkpoint activation
num_t = length(t_span);
d_G1 = para_unknown(end-2) ;  d_S = para_unknown(end-1) ;  d_G2M = para_unknown(end);
m12 = 1-d_G1;
m1d = 0;
m23 = 1-d_S;
m2d = 0;
m31 = 2*(1-d_G2M);
m3d = 0;
[G1_cdf, tilde_G1_cdf, tilde_S_cdf, tilde_G2M_cdf, G1_pdf, f_1to2, f_2to3, f_3to1,f_1tod,f_2tod, f_3tod] = myfun_treated_PDF_TransitionTime(ChoiceOfDist_untreated, ...
    ChoiceOfDist_treated_Dd,ChoiceOfDist_treated_D1, ChoiceOfDist_treated_U, [a;b_scale], para_unknown, t_span, StepSize);
mat_G_tilde = zeros(3*num_gen+1, 3*num_gen+1 , num_t);
%single dose
for k = 1:num_gen
    mat_G_tilde( 1+ 3*(k-1),  1+ 3*(k-1) ,:   ) = tilde_G1_cdf;
    mat_G_tilde( 2+ 3*(k-1),  2+ 3*(k-1),:   )  = tilde_S_cdf; 
    mat_G_tilde( 3+ 3*(k-1),  3+ 3*(k-1),:   )  = tilde_G2M_cdf;      
end     
[~, D_pdf, D_cdf]  =   PdfcdfCal_treated('Exponential',  t_span , 4,  [] );
mat_G_tilde(3*num_gen+1 , 3*num_gen+1,:   )  =  D_cdf; 
% If the imediate new phase is G1 at time of the treatment, then the PDF of
% that phase will be that of normal G1 

%% matrix M
temp = diag( repmat([m12, m23, m31],1,num_gen )  ,1);
mat_m_tilde = [temp(1: 3*num_gen, 1: 3*num_gen),  repmat( [m1d; m2d; m3d], num_gen,1 )    ; 
    zeros(1, 3*num_gen+1 )  ];
I_mat_G_tilde = zeros(3*num_gen+1,3*num_gen+1 , num_t);
for k = 1:num_gen
        I_mat_G_tilde( 1+ 3*(k-1),  1+ 3*(k-1) ,:   ) = 1- tilde_G1_cdf;
        I_mat_G_tilde( 2+ 3*(k-1),  2+ 3*(k-1),:   )  = 1- tilde_S_cdf; 
        I_mat_G_tilde( 3+ 3*(k-1),  3+ 3*(k-1),:   )  = 1- tilde_G2M_cdf; 
end
I_mat_G_tilde( 3*num_gen+1 , 3*num_gen+1,:   )  =  1 - D_cdf;
% I_mat_3D( 3*num_gen+1 , 3*num_gen+1,: )  = ones(1, num_t);

%% matrix Gm tilde
mat_Gm_tilde = zeros(3*num_gen+1,3*num_gen+1 , num_t);
for i = 1: size(mat_m_tilde,1)
        if  i < 3*num_gen-1
            mat_Gm_tilde(i,i+1 ,: ) = reshape( mat_G_tilde( i,  i,:   ) , 1, num_t  ) * mat_m_tilde(i,i+1); %mat_m(i,j) scalar value
        end
        if  i < 3*num_gen
            mat_Gm_tilde(i, end ,: ) = reshape( mat_G_tilde( i,  i,:   ) , 1, num_t  ) * mat_m_tilde(i,end); %mat_m(i,j) scalar value
        end
end

end