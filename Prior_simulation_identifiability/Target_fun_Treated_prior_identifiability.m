function  [exception, Target_value,  diff_norm_sum ]= Target_fun_Treated_prior_identifiability(theta, Model_Index ,outfile, model_obj,lb,ub, ...
    synthetic_data, para_change_exclude_idx,  para_fixed)
% SimData is a matrix with dimension of Time * 3


theta_full = assemble_theta(theta, para_fixed, para_change_exclude_idx); 
model_obj.para_unknown = theta_full ;
[exception, SimData, ~, ~,~, ~,   ~ , ~,~]  =   model_obj.model_simulation(  );
if isnan(SimData)
    Target_value = NaN;
    diff_norm_sum = [];
else
    num_dataset = size(SimData,2);
    ThetaNoise = theta(end-2:end);
    
    %fprintf(outfile, '\t Sim by new para: G1= %.3f , S = %.3f , G2M = %.3f \n', SteadyState(1), SteadyState(2), SteadyState(3) );
    %% specify the prior distribution for each parameter
    % The variances of normal random variable used to 2.
    if Model_Index == 1
        Current_Data  = [39.78    35.02    25.20
            28.1102   26.6131   45.2766
            18.1365   20.6562   61.2073
            17.4534   16.9032   65.6434
            14.9088   13.8842   71.2070
            ]*0.01  ;  % Time * 3
        num_Data  = size(Current_Data,1 );
        
        % b: shape
        % lambda : scale parameter
        prior_vec = zeros(length(theta), 1);
        
        for i = 1: length(theta)
            
            prior_vec(i) = unifpdf( theta(i), lb(i), ub(i)  ) ;
        end
        
    elseif  Model_Index == 3
        prior_vec = zeros(length(theta), 1);
        
        for i = 1: length(theta)
            
            prior_vec(i) = unifpdf( theta(i), lb(i), ub(i)  ) ;
        end
    elseif Model_Index == 41
        Current_Data  = [   40.5645161290322    31.209677419354865    28.70967741935491
            32.5	18.79032258	47.82258065
            16.37096774	12.82258065	70.24193548
            25.32258065	30.88709677	44.43548387
            17.5	51.69354839	29.91935484 ] *0.01;
        Current_Data  = Current_Data - (sum(Current_Data ,2)-1)/3;
        num_Data  = size(Current_Data,1 );
        
        % b: shape
        % lambda : scale parameter
        prior_vec = zeros(length(theta), 1);
        
        for i = 1: length(theta)
            
            prior_vec(i) = unifpdf( theta(i), lb(i), ub(i)  ) ;
        end
        %%
    elseif Model_Index == 42
        Current_Data  = [  40.40322581	31.0483871	28.5483871
            59.43548387	38.9516129	2.5
            63.30645161	29.91935484	4.758064516
            72.5	21.69354839	4.435483871
            86.69354839	10.16129032	3.306451613 ] *0.01;
        Current_Data  = Current_Data - (sum(Current_Data ,2)-1)/3;
        num_Data  = size(Current_Data,1 );
        % b: shape
        % lambda : scale parameter
        prior_vec = zeros(length(theta), 1);
        
        for i = 1: length(theta)
            
            prior_vec(i) = unifpdf( theta(i), lb(i), ub(i)  ) ;
        end
        
    end

    Current_Data_full  =  synthetic_data; 
    Current_Data  = interp1( model_obj.t_span_AfterTreatment , Current_Data_full',  model_obj.data_time ,'pchip') ; 
   
    num_Data =  size(Current_Data,1 ); 

    prior_val  = sum( log(prior_vec));
    
    Target_LE_temp_col =  zeros(num_dataset,1);
    temp_diff = SimData - Current_Data;
    diff_norm_sum = 0;
    %% calculate likehood for each dataset
    for i = 1:num_dataset
        diff_norm = sumsqr(temp_diff(:,i));
        diff_norm_sum  = diff_norm_sum  + diff_norm;
        fprintf(outfile, '\t diff_norm (%d) : %.3f  \n', i,  vecnorm(temp_diff(:,i),1)  );
        
        Noise_Term = ThetaNoise(i); %variance
%         Sigma =  sqrt(Noise_Term )* eye( num_Data );
%         Target_LE_temp_1 = mvnpdf(SimData(:,i), Sigma, Current_Data(:,i) )
%         %% 
        C_vec = (1./ (2*pi*Noise_Term) ).^(num_Data./2); %vector % constant factor 2pi doesn't influence the result.
        Target_LE_temp = C_vec*...
            exp(- 1/(2*Noise_Term) *diff_norm)    ;  % eliminate the multiplicative constant
        Target_LE_temp_col(i) = Target_LE_temp;
    end
    
    %% multiply likehood with prior%
    Target_value   =  sum( log(Target_LE_temp_col ) ) + prior_val;
    %Target_value   =  sum( log(Target_LE_temp ) );
end
end