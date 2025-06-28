function mat_M = Sphasedrug_AfterTreatmentBranchingProcess_G1_2( G1_block, mat_m, S_arrest_UR,  S_arrest_FR,  tilde_G2M, G1_2 , G1_cdf, S_cdf, G2M_cdf,...
    num_GenOfHealthycells ,StepSize,  t_span_AfterTreatment , num_RroundsOfTreatedcells) 
%M = readtable('Branching process after treatment-G2M phase drug.xlsx','Sheet','Transition matrixG1_2','Range','B2:L12') ;

True_num_GenOfHealthycells  = num_GenOfHealthycells -1;
% Initialize a matrix filled with symbolic zeros
% Create symbolic variable
%syms q_1 q_2 q_3 md_UR md_FR md_MS md_UR_G1
TotalState  = size( mat_m, 1);
 
%%
% calculate Gm, I-G
% create the tilde_G1, G2_arrest
num_t_BeforeTreatment  = length(  G1_2);
mat_G = zeros(TotalState , TotalState , num_t_BeforeTreatment );
num_t_AfterTreatment =  length(  t_span_AfterTreatment );

mat_G(1,1,:) = G1_2;
mat_G(2,2,: ) = G1_block ;
mat_G(3,3,: ) = S_cdf;
mat_G(4,4,: ) = G2M_cdf;
mat_G(5,5,: ) = S_arrest_UR;
mat_G(6,6,: ) =  S_arrest_FR;
mat_G(7,7,: ) = tilde_G2M;
%% Repeated block
for i = 1: num_RroundsOfTreatedcells-1
    startnum = 8 +( i-1 )*7;
    endnum =  8 + i*7-1;
    mat_G(startnum , startnum, :     ) =  G1_cdf ;
    mat_G(startnum+1   : endnum, startnum+1 :  endnum,:    ) = mat_G( 2 :7, 2:7 ,:);
end


% mat_G(8,8,: ) = G1_cdf;
% mat_G(9,9,: ) = G1_block ;
% mat_G(10,10,: ) = S_cdf;
% mat_G(11,11,: ) = G2M_cdf;
% mat_G(12,12,: ) = S_arrest_UR;
% mat_G(13,13,: ) =  S_arrest_FR;
% mat_G(14,14,: ) = tilde_G2M;
%%
if True_num_GenOfHealthycells > 0
    for i = 1:True_num_GenOfHealthycells
        startnum = 7*num_RroundsOfTreatedcells + 1  + (i-1)*3;
        mat_G(startnum , startnum,: ) = G1_cdf;
        mat_G(startnum+1, startnum+1,: ) = S_cdf;
        mat_G(startnum+2, startnum+2,: ) = G2M_cdf;
    end
end
mat_G = mat_G(:,:,1:num_t_AfterTreatment );

[~, ~, G_apop_cdf]  =   PdfcdfCal_treated('Exponential',  t_span_AfterTreatment ,5 ,  [] );
mat_G(end,  end,:   )  =    G_apop_cdf; % time for eliminating the apoptotic cells



mat_I = reshape( repmat( eye(TotalState), 1,num_t_AfterTreatment) ,  TotalState, TotalState, num_t_AfterTreatment );
I_mat_G =  mat_I - mat_G ;

mat_Gm = zeros(TotalState , TotalState ,num_t_AfterTreatment) ;
for i = 1: TotalState
        for j =  1: TotalState
            mat_Gm(i,j ,: ) = reshape( mat_G( i,  i,:   ) , 1, num_t_AfterTreatment   ) * mat_m(i,j); %mat_m(i,j) scalar value
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
        
        [ next_mat_M_term_k, ~ ]  =  fun_NeumannSeriesSolvingForM(mat_Gm, TotalState, [],  I_mat_G, ctr, StepSize,FlagTreatment);
        next_mat_M = curr_mat_M+next_mat_M_term_k;
    else  %k >= 2
        curr_mat_M = next_mat_M;
        curr_T_k_iterate = next_T_k_iterate;
        [ next_mat_M_term_k, iterate_T_k ]  =  fun_NeumannSeriesSolvingForM(mat_Gm,TotalState, curr_T_k_iterate,  I_mat_G, ctr, StepSize,FlagTreatment);
        next_T_k_iterate  = iterate_T_k ; %update iterate
        next_mat_M  = curr_mat_M + next_mat_M_term_k;
    end
    ctr = ctr +1;
end
%%
mat_M =  next_mat_M;

end