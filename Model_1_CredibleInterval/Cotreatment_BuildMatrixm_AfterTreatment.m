function   mat_m = Cotreatment_BuildMatrixm_AfterTreatment(obj)

T = readcell('transition_matrix.xlsx', 'Sheet' ,'co-treatment',  'Range', 'B3:U22');%ReadVariableNames=false
%T = standardizeMissing(T, {});
% Define values for each variable (for the sake of this example, let's use random numbers)
%values = rand(2, 3); % create a 3x3 matrix of random values

[rows, cols] = size(T);
numTotalStatesPerRound =  rows -1 ;
True_num_GenOfHealthycells  = obj.num_GenOfHealthycells -1;
TotalState =  numTotalStatesPerRound * obj.num_RroundsOfTreatedcells + True_num_GenOfHealthycells*3+1;
mat_m = zeros(TotalState , TotalState );
para_label = obj.para_label ;
 
num_para = length(para_label );
values = rand(num_para ,1);
EC50_d_G2M =  values(21)  ;
n_G2M = values(19)  ;
n_S = values(6)  ;
EC50_d_S =  values(37) ;

Dose_G2M = obj.Dose_1thDrug;
Dose_S = obj.Dose_2thDrug;

m_d_fun_G2M = @(m_max) m_max * ((Dose_G2M / EC50_d_G2M)^n_G2M) / (1 + (Dose_G2M/ EC50_d_G2M)^n_G2M);
m_d_fun_S = @(m_max) m_max * ((Dose_S / EC50_d_S)^n_S) / (1 + (Dose_S  / EC50_d_S)^n_S);

% Iterate through the table
for i = 1:rows
    for j = 1:cols
        
        expr = T{i,j}; % Extract variable name from the table
        if ismissing( expr )
            mat_m(i,j) = 0;
        else 
            expr_substituted  =  expr;
            for label = 1: length(para_label)
                newLabel = strrep( para_label(label), '^{max}', '');
                newLabel = strrep(  newLabel, '$', '');
                if contains( expr_substituted, newLabel )
                    if contains( para_label(label),   '^{max}') && label > 21 %S
                        truevalue = m_d_fun_S(values(label));
                    elseif  contains( para_label(label),   '^{max}') && label <= 21 %G2M
                         truevalue = m_d_fun_G2M(values(label));
                    else 
                         truevalue  = values(label);
                    end
                    expr_substituted = strrep(expr_substituted, newLabel, num2str(truevalue   )  );
                end
            end
            expr_substituted = strrep(expr_substituted, '$', '');
            if  contains( expr_substituted, '2('  )
               expr_substituted = strrep(expr_substituted, '2(', '2*(');
            end
            % Evaluate the substituted expression
            result = eval(expr_substituted);
        
            mat_m(i,j) =  result ;
        end
        
    end
    
end


%% repeated block

for i = 1: num_RroundsOfTreatedcells-1
    startnum = numTotalStatesPerRound + 1 +( i-1 )*numTotalStatesPerRound;
    endnum = numTotalStatesPerRound +1 + i*numTotalStatesPerRound;
    mat_m(startnum   : endnum, startnum :  endnum   ) = mat_m( 1: numTotalStatesPerRound + 1 , 1: numTotalStatesPerRound + 1);
    mat_m(startnum   : endnum, end  ) =  mat_m( 1: numTotalStatesPerRound + 1, end);
end    

if True_num_GenOfHealthycells > 0
    %mat_m(11, 15) = 2;
    %mat_m(14, 15) = 2*(1- md_UR_G2M);
    for i = 1:True_num_GenOfHealthycells
        startnum = numTotalStatesPerRound *num_RroundsOfTreatedcells + 1 + (i-1)*3;
        mat_m(startnum , startnum+1) = 1;
        mat_m(startnum+1 , startnum+2) = 1;
        if startnum+2 == TotalState-1
            break
        end
        mat_m(startnum+2 , startnum+3) = 2;
    end
end

end