function mat_M = G2Mphasedrug_AfterTreatmentBranchingProcess_G2M_2( mat_m, G2_arrest_UR,  G2_arrest_FR,  G2_arrest_MS ,  tilde_G1,   tilde_S, G2M_2 , G1_cdf, S_cdf, G2M_cdf,  ...
num_GenOfHealthycells ,StepSize,  t_span_AfterTreatment , num_RroundsOfTreatedcells ) 
%M = readtable('Branching process after treatment-G2M phase drug.xlsx','Sheet','Transition matrixG1_2','Range','B2:L12') ;

%%

True_num_GenOfHealthycells  = num_GenOfHealthycells -1;
% Initialize a matrix filled with symbolic zeros
% Create symbolic variable
%syms q_1 q_2 q_3 md_UR md_FR md_MS md_UR_G1
% mat_G1 =  mat_m(3:end, 3:end);
% [~,size1] = size(mat_G1);
% mat_m =  zeros(size1 +8,  size1 +8 );
% mat_m(1:8, 1:8)  = mat_G1(1:8, 1:8);
% mat_m(9:end, 9:end) =  mat_G1;
mat_m = mat_m(3:end, 3:end);
TotalState_G2M = size(mat_m,1);
%% calculate Gm, I-G
num_t_BeforeTreatment = length( G2M_2 );
mat_G = zeros(TotalState_G2M ,TotalState_G2M , num_t_BeforeTreatment   );
num_t_AfterTreatment =  length(  t_span_AfterTreatment );


mat_G(1,1,:) =  G2M_2;
mat_G(2,2, : ) = G1_cdf;
mat_G(3,3,: ) = G2_arrest_UR;
mat_G(4,4,: ) =  G2_arrest_FR;
mat_G(5,5,: ) = G2_arrest_MS;
mat_G(6,6,: ) = tilde_G1;
mat_G(7,7,: ) = tilde_S;
mat_G(8,8, :) = S_cdf;
mat_G(9,9, :) = G2M_cdf;

for i = 1: num_RroundsOfTreatedcells-1
    startnum = 10 +( i-1 )*8;
    endnum = 10 + i*8-1;
    mat_G(startnum   : endnum, startnum :  endnum,:    ) = mat_G( 2:9, 2:9,:);
end


newG2MWithoutdrug =  8*num_RroundsOfTreatedcells  +1;
%mat_G(10:17,10:17,: ) = mat_G(2:9, 2:9,:);
if True_num_GenOfHealthycells > 0
    for i = 1:True_num_GenOfHealthycells
        startnum = newG2MWithoutdrug+1 + (i-1)*3;
        mat_G(startnum , startnum,: ) = G1_cdf;
        mat_G(startnum+1, startnum+1,: ) = S_cdf;
        mat_G(startnum+2, startnum+2,: ) = G2M_cdf;
    end
end
mat_G = mat_G(:,:,1:num_t_AfterTreatment );
[~, ~, G_apop_cdf]  =   PdfcdfCal_treated('Exponential',  t_span_AfterTreatment  ,5 ,  [] );
mat_G(end,  end,:   )  =    G_apop_cdf; % time for eliminating the apoptotic cells

mat_I = reshape( repmat( eye(TotalState_G2M), 1,num_t_AfterTreatment) ,  TotalState_G2M, TotalState_G2M , num_t_AfterTreatment );
I_mat_G =  mat_I - mat_G ;

mat_Gm = zeros(TotalState_G2M , TotalState_G2M , num_t_AfterTreatment  );
for i = 1: TotalState_G2M
        for j =  1: TotalState_G2M
            mat_Gm(i,j ,: ) = reshape( mat_G( i,  i,:   ) , 1, num_t_AfterTreatment  ) * mat_m(i,j); %mat_m(i,j) scalar value
            %mat_Gm(i, end ,: ) = reshape( mat_G( i,  i,:   ) , 1, num_t  ) * mat_m(i,end); %mat_m(i,j) scalar value
        end    
end

%%
ctr = 0;
FlagTreatment = 1; 
while ctr<=  TotalState_G2M
    if ctr  == 0
        % k=0
        next_mat_M = I_mat_G;
    elseif ctr == 1    % k = 1
        curr_mat_M = next_mat_M;
        next_T_k_iterate = mat_Gm;
        
        [ next_mat_M_term_k, ~ ]  =  fun_NeumannSeriesSolvingForM(mat_Gm, TotalState_G2M, [],  I_mat_G, ctr, StepSize,FlagTreatment);
        next_mat_M = curr_mat_M+next_mat_M_term_k;
    else  %k >= 2
        curr_mat_M = next_mat_M;
        curr_T_k_iterate = next_T_k_iterate;
        [ next_mat_M_term_k, iterate_T_k ]  =  fun_NeumannSeriesSolvingForM(mat_Gm,TotalState_G2M, curr_T_k_iterate,  I_mat_G, ctr, StepSize,FlagTreatment);
        next_T_k_iterate  = iterate_T_k ; %update iterate
        next_mat_M  = curr_mat_M + next_mat_M_term_k;
    end
    ctr = ctr +1;
end
%%
next_mat_M(next_mat_M<0) = 0;
mat_M =  next_mat_M;

end