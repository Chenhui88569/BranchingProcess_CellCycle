function  mat_m = SequentialTreatment_BuildMatrixm_AfterTreatment(obj)
Dose_G2M =  obj.Dose_1thDrug;
Dose_S = obj.Dose_2thDrug;

EC50_d_G2M =  values(21)  ;
n_G2M = values(19)  ;
n_S = values(6)  ;
EC50_d_S =  values(37) ;

m_d_fun_G2M = @(m_max) m_max * ((Dose_G2M / EC50_d_G2M)^n_G2M) / (1 + (Dose_G2M/ EC50_d_G2M)^n_G2M);
m_d_fun_S = @(m_max) m_max * ((Dose_S / EC50_d_S)^n_S) / (1 + (Dose_S  / EC50_d_S)^n_S);

filename = 'transition_matrix.xlsx';
if obj.TreatmentEffectOnPhases_1thDrug == 3 && obj.TreatmentEffectOnPhases_2thDrug == 2 % G2M phase first
    mat_m_template  = readcell(filename, 'Sheet' ,'G2M phase first',  'Range', 'B2:AD30');
    TotalState  = size(mat_m_template,1);
    mat_m = zeros(TotalState , TotalState );
  
    for i = 1:rows
        for j = 1:cols
            
            expr = mat_m_template{i,j}; % Extract variable name from the table
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
    
    
    
 
elseif obj.TreatmentEffectOnPhases_1thDrug == 2 && obj.TreatmentEffectOnPhases_2thDrug == 3
    mat_m_template  = readcell(filename, 'Sheet' ,'S phase first',  'Range', 'B2:AB28');
    TotalState  = size(mat_m_template,1);
    mat_m = zeros(TotalState , TotalState );
  
    for i = 1:rows
        for j = 1:cols
            
            expr = mat_m_template{i,j}; % Extract variable name from the table
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
    
    
    

elseif obj.TreatmentEffectOnPhases_1thDrug == 3 && obj.TreatmentEffectOnPhases_2thDrug == 0 % normal monotherapy
    mat_m = G2Mphasedrug_BuildMatrixm_AfterTreatment(Dose, obj, obj.num_GenOfHealthycells, obj.num_RroundsOfTreatedcells);

elseif obj.TreatmentEffectOnPhases_1thDrug == 2 && obj.TreatmentEffectOnPhases_2thDrug == 0 % normal monotherapy
    mat_m = Sphasedrug_BuildMatrixm_AfterTreatment(Dose, obj, obj.num_GenOfHealthycells, obj.num_RroundsOfTreatedcells);

end