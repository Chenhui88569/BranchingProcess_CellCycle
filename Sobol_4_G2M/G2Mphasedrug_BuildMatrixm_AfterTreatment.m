function  mat_m = G2Mphasedrug_BuildMatrixm_AfterTreatment(Dose, obj, num_GenOfHealthycells, num_RroundsOfTreatedcells)

q_1 = obj.para_unknown(1);
q_2 = obj.para_unknown(2);
q_3 = obj.para_unknown(3);
q_4 = obj.para_unknown(4);
%q_1 =  q_1./(q_1+   q_2+ q_3  )   ;    q_2 =  q_1./(q_1+   q_2+ q_3  )   ;  q_3 =  q_1./(q_1+   q_2+ q_3  )   ;
m_d_FR_max =  obj.para_unknown(5);
m_d_UR_max =  obj.para_unknown(6);
m_d_MS_max =  obj.para_unknown(7);
m_d_UR_G1_max=  obj.para_unknown(8);
m_d_UR_S_max =  obj.para_unknown(9);
n_G2M =  obj.para_unknown(19);
EC50_d =  obj.para_unknown(20);

m_d_fun = @(m_max) m_max * ((Dose / EC50_d)^n_G2M) / (1 + (Dose / EC50_d)^n_G2M);

md_UR = m_d_fun(m_d_UR_max );
md_FR = m_d_fun(m_d_FR_max );
md_MS = m_d_fun(m_d_MS_max);
md_UR_G1 =  m_d_fun(m_d_UR_G1_max);
md_UR_S =  m_d_fun(m_d_UR_S_max);

True_num_GenOfHealthycells  = num_GenOfHealthycells -1;

TotalState  = 2 + 8*num_RroundsOfTreatedcells  +1 + True_num_GenOfHealthycells*3 +1 ;
mat_m = zeros(TotalState , TotalState );

% Fill the matrix based on your table
mat_m(1, 2) = 1;
mat_m(2, 3) = 1;

mat_m(3, 4) = 2 * (1 - q_1 - q_2 - q_3);
mat_m(3, 5) = q_1;
mat_m(3, 6) = q_2;
mat_m(3, 7) = q_3;

mat_m(4, 10) = 1;

mat_m(5, 8) = 2 * (1 - md_UR);
mat_m(5, end) = md_UR;

mat_m(6, 4) = 2 * (1 - md_FR);
mat_m(6, end) = md_FR;

mat_m(7, 8) = 1 - md_MS;
mat_m(7, end) = md_MS;

mat_m(8, 9) =1- md_UR_G1-q_4;
mat_m(8, 10) = q_4;
mat_m(8, end) = md_UR_G1;

mat_m(9, 11) = 1-md_UR_S ;
mat_m(9, end) = md_UR_S ;

mat_m(10, 11) = 1 ; %11 : G2M phase 



for i = 1: num_RroundsOfTreatedcells-1
    startnum = 11 +( i-1 )*8;
    endnum = 11 + i*8;
    mat_m(startnum   : endnum, startnum :  endnum   ) = mat_m( 3:11, 3:11);
    mat_m(startnum:endnum, end  ) = mat_m( 3:11, end);
    
end


newG2MWithoutdrug = 2 + 8*num_RroundsOfTreatedcells  +1;
if True_num_GenOfHealthycells > 0
    mat_m(newG2MWithoutdrug  , newG2MWithoutdrug +1  ) = 2;
    for i = 1:True_num_GenOfHealthycells
        startnum = newG2MWithoutdrug +1  + (i-1)*3;
        mat_m(startnum , startnum+1) = 1;
        mat_m(startnum+1 , startnum+2) = 1;
        if startnum+2 == TotalState-1
            break
        end
        mat_m(startnum+2 , startnum+3) = 2;
    end
end