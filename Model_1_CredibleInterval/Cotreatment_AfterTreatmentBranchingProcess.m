function mat_M = Cotreatment_AfterTreatmentBranchingProcess( obj, mat_m ,  mat_G  ,   g_est2_cdf, phase)
num_t_AfterTreatment =  obj.t_span_AfterTreatment;

%% restructure mat G

if  phase == 0 % G2M phase
    mat_m  = mat_m(6: end, 6 , end);
    mat_G  =  mat_G(6: end, 6 , end,: );
    mat_G(1,1,:) =   g_est2_cdf(1:num_t_AfterTreatment ); 
    TotalState = size(mat_m,1);
elseif  phase == 1 % G1 phase
    mat_G(1,1,:) =   g_est2_cdf(1:num_t_AfterTreatment ); 
    TotalState = size(mat_m,1);
elseif  phase == 2 % S 
    mat_m  = mat_m(3: end, 3 , end);
    mat_G  =  mat_G(3: end, 3 , end,: );
    mat_G(1,1,:) =   g_est2_cdf(1:num_t_AfterTreatment ); 
    TotalState = size(mat_m,1);
end
%% calculate Gm, I-G
mat_I = reshape( repmat( eye(TotalState ), 1, num_t_AfterTreatment ) ,  TotalState ,TotalState  , num_t_AfterTreatment  );
I_mat_G =  mat_I - mat_G ;

mat_Gm = zeros(TotalState  ,TotalState  , num_t_AfterTreatment );
for i = 1: TotalState 
        for j =  1: TotalState 
            mat_Gm(i,j ,: ) = reshape( mat_G( i,  i,:   ) , 1, num_t_AfterTreatment  ) * mat_m(i,j); %mat_m(i,j) scalar value
            %mat_Gm(i, end ,: ) = reshape( mat_G( i,  i,:   ) , 1, num_t  ) * mat_m(i,end); %mat_m(i,j) scalar value
        end
    
end

%%
ctr = 0;
FlagTreatment = 1; 
while ctr<=  TotalState 
    if ctr  == 0
        % k=0
        next_mat_M = I_mat_G;
    elseif ctr == 1    % k = 1
        curr_mat_M = next_mat_M;
        next_T_k_iterate = mat_Gm;
        
        [ next_mat_M_term_k, ~ ]  =  fun_NeumannSeriesSolvingForM(mat_Gm, TotalState , [],  I_mat_G, ctr, StepSize,FlagTreatment);
        next_mat_M = curr_mat_M+next_mat_M_term_k;
    else  %k >= 2
        curr_mat_M = next_mat_M;
        curr_T_k_iterate = next_T_k_iterate;
        [ next_mat_M_term_k, iterate_T_k ]  =  fun_NeumannSeriesSolvingForM(mat_Gm,TotalState , curr_T_k_iterate,  I_mat_G, ctr, StepSize,FlagTreatment);
        next_T_k_iterate  = iterate_T_k ; %update iterate
        next_mat_M  = curr_mat_M + next_mat_M_term_k;
    end
    ctr = ctr +1;
end
%%
mat_M =  next_mat_M;

end