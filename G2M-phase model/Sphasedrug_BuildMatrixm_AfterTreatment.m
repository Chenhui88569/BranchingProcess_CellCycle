function  mat_m = Sphasedrug_BuildMatrixm_AfterTreatment(Dose, obj, num_GenOfHealthycells, num_RroundsOfTreatedcells)

q_1 = obj.para_unknown(1);
q_2 = obj.para_unknown(2);
q_3 = obj.para_unknown(3);
%q_1 =  q_1./(q_1+   q_2+ q_3  )   ;    q_2 =  q_1./(q_1+   q_2+ q_3  )   ;  q_3 =  q_1./(q_1+   q_2+ q_3  )   ;
m_d_FR_max =  obj.para_unknown(4);
m_d_UR_max =  obj.para_unknown(5);
m_d_UR_G2M_max =  obj.para_unknown(6);
m_d_G1block_max =  obj.para_unknown(7);
n_S =  obj.para_unknown(15);
EC50_d =  obj.para_unknown(16);

m_d_fun = @(m_max) m_max * ((Dose / EC50_d)^n_S) / (1 + (Dose / EC50_d)^n_S);

md_UR = m_d_fun(m_d_UR_max );
md_FR = m_d_fun(m_d_FR_max );
md_UR_G2M =  m_d_fun(m_d_UR_G2M_max);
m_d_G1block =  m_d_fun(m_d_G1block_max);

True_num_GenOfHealthycells  = num_GenOfHealthycells -1;

TotalState =  7*num_RroundsOfTreatedcells + True_num_GenOfHealthycells*3+1;
mat_m = zeros(TotalState , TotalState );

% Fill the matrix based on your table
mat_m(1, 2) = q_1;
mat_m(1, 3) = 1-q_1;
mat_m(2, 3) = 1- m_d_G1block;
mat_m(2, end) =  m_d_G1block;

mat_m(3, 4) = (1 - q_2 - q_3 );
mat_m(3, 5) = q_2;
mat_m(3, 6) = q_3;

mat_m(4, 8) = 2;

mat_m(5, 7) = 2 * (1 - md_UR);
mat_m(5, end) = md_UR;

mat_m(6, 4) =  (1 - md_FR);
mat_m(6, end) = md_FR;
mat_m(7, 8) = 2*(1 -  md_UR_G2M);
mat_m(7, end) =  md_UR_G2M;

%% repeated block

for i = 1: num_RroundsOfTreatedcells-1
    startnum = 8 +( i-1 )*7;
    endnum = 8 + i*7;
    mat_m(startnum   : endnum, startnum :  endnum   ) = mat_m( 1:8, 1:8);
    mat_m(startnum   : endnum, end  ) =  mat_m( 1:8, end);
end

% mat_m(8, 9) = q_1;
% mat_m(8, 10) = 1-q_1;
% mat_m(9, 10) = 1- m_d_G1block;
% mat_m(9, end) =  m_d_G1block;
% 
% mat_m(10, 11) = (1 - q_2 - q_3 );
% mat_m(10, 12) = q_2;
% mat_m(10, 13) = q_3;
% 
% mat_m(12, 14) = 1 - md_UR;
% mat_m(12, end) = md_UR;
% 
% mat_m(13, 11) =  1 - md_FR;
% mat_m(13, end) = md_FR;
% 
% mat_m(14, end) = md_UR_G2M ;

%mat_m(:,  7*num_RroundsOfTreatedcells + True_num_GenOfHealthycells*3) = 0;
    

if True_num_GenOfHealthycells > 0
    %mat_m(11, 15) = 2;
    %mat_m(14, 15) = 2*(1- md_UR_G2M);
    for i = 1:True_num_GenOfHealthycells
        startnum = 7*num_RroundsOfTreatedcells + 1 + (i-1)*3;
        mat_m(startnum , startnum+1) = 1;
        mat_m(startnum+1 , startnum+2) = 1;
        if startnum+2 == TotalState-1
            break
        end
        mat_m(startnum+2 , startnum+3) = 2;
    end
end