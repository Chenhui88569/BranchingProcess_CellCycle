function [ next_mat_M_term_k, iterate_T_k ]  =  fun_NeumannSeriesSolvingForM(mat_Gm, TotalState, curr_T_k_iterate,  I_mat_G, ctr, StepSize, FlagTreatment)


[num_dim,~,num_t ] =  size(mat_Gm);
if ctr  == 1
    if  FlagTreatment ==0
        next_mat_M_term_k  =   zeros(TotalState,TotalState, num_t);
        for i = 1:num_dim
            %for  j  = 1: size(mat_m,1)
            if i+ctr <  num_dim % curr_T_k_iterate(i,i+ctr-1,:),1,num_t) is non-zero, computed in the pevious interation.
                temp_mat_M2_conv_full =  conv( reshape(I_mat_G(i+ctr,i+ctr,:),1,num_t),   TakeDerivative(reshape(  mat_Gm(i,i+ctr,:),1,num_t),StepSize ),'full')*StepSize;
                % temp_mat_M2_conv_full =  conv( reshape(  iterate_T_k(i,j,:),1,num_t),  reshape(I_mat_G(j,j,:),1,num_t) ,'full')*StepSize;
                temp_mat_M2_conv_same =   temp_mat_M2_conv_full(1:num_t);
                next_mat_M_term_k(i,i+ctr,:) =  temp_mat_M2_conv_same ;
            end
            %end
        end
        iterate_T_k = [];
    end
    %% With treatment
    if  FlagTreatment == 1
        next_mat_M_term_k  =   zeros(TotalState,TotalState, num_t);
        for i = 1:num_dim
            for  j  = 1:num_dim
                temp_mat_M2_conv_full =  conv( reshape(I_mat_G(j,j,:),1,num_t),   TakeDerivative(reshape(  mat_Gm(i,j,:),1,num_t),StepSize ),'full')*StepSize;
                % temp_mat_M2_conv_full =  conv( reshape(  iterate_T_k(i,j,:),1,num_t),  reshape(I_mat_G(j,j,:),1,num_t) ,'full')*StepSize;
                temp_mat_M2_conv_same =   temp_mat_M2_conv_full(1:num_t);
                next_mat_M_term_k(i,j,:) =  temp_mat_M2_conv_same ;
            end
        end
        iterate_T_k = [];
    end
    
end


if ctr >= 2
    if FlagTreatment ==0
        % used in the inner loop
        iterate_T_k =  zeros(TotalState,TotalState , num_t);
        for i = 1: num_dim
            if i+ctr <  num_dim % curr_T_k_iterate(i,i+ctr-1,:),1,num_t) is non-zero, computed in the pevious interation.
                temp_mat_M1_conv_full = conv( reshape( curr_T_k_iterate(i,i+ctr-1,:),1,num_t),   TakeDerivative(reshape(  mat_Gm(i+ctr-1,i+ctr,:),1,num_t),StepSize),'full')*StepSize;
                temp_mat_M1_conv_same =  temp_mat_M1_conv_full(1:num_t);
                iterate_T_k(i,i+ctr,:) = temp_mat_M1_conv_same;
            end
        end
        
        next_mat_M_term_k  =   zeros(TotalState,TotalState, num_t);
        for i = 1:num_dim
            %for  j  = 1: size(mat_m,1)
            if i+ctr <  num_dim % curr_T_k_iterate(i,i+ctr-1,:),1,num_t) is non-zero, computed in the pevious interation.
                temp_mat_M2_conv_full =  conv( reshape(I_mat_G(i+ctr,i+ctr,:),1,num_t),   TakeDerivative(reshape(  iterate_T_k(i,i+ctr,:),1,num_t), StepSize),'full')*StepSize;
                % temp_mat_M2_conv_full =  conv( reshape(  iterate_T_k(i,j,:),1,num_t),  reshape(I_mat_G(j,j,:),1,num_t) ,'full')*StepSize;
                temp_mat_M2_conv_same =   temp_mat_M2_conv_full(1:num_t);
                next_mat_M_term_k(i,i+ctr,:) =  temp_mat_M2_conv_same ;
            end
            %end
        end
    end
     %% with treatment
     
    use_gpu = 1;  % or 0


    if FlagTreatment == 1
        if use_gpu == 1
            % Allocate on GPU
            iterate_T_k = gpuArray.zeros(TotalState, TotalState, num_t);
            
            for i = 1:num_dim
                for j = 1:num_dim
                    temp_mat_M1_conv_full = gpuArray.zeros(1, 2*num_t - 1);
                    for k = 1:num_dim
                        a = reshape(curr_T_k_iterate(i,k,:), 1, num_t);
                        b = TakeDerivative(reshape(mat_Gm(k,j,:), 1, num_t), StepSize);
                        temp_mat_M1_conv_full = temp_mat_M1_conv_full + convn(a, b, 'full') * StepSize;
                    end
                    temp_mat_M1_conv_same = temp_mat_M1_conv_full(1:num_t);
                    temp_mat_M1_conv_same(temp_mat_M1_conv_same < 0) = 0;
                    iterate_T_k(i,j,:) = temp_mat_M1_conv_same;
                end
            end
        
            next_mat_M_term_k_gpu = gpuArray.zeros(TotalState, TotalState, num_t);
            
            for i = 1:num_dim
                for j = 1:num_dim
                    a = reshape(I_mat_G(j,j,:), 1, num_t);
                    b = TakeDerivative(reshape(iterate_T_k(i,j,:), 1, num_t), StepSize);
                    temp_mat_M2_conv_full = convn(a, b, 'full') * StepSize;
                    temp_mat_M2_conv_same = temp_mat_M2_conv_full(1:num_t);
                    temp_mat_M2_conv_same(temp_mat_M2_conv_same < 0) = 0;
                    next_mat_M_term_k_gpu(i,j,:) = temp_mat_M2_conv_same;
                end
            end
        
            % Convert final result to CPU
            next_mat_M_term_k = gather(next_mat_M_term_k_gpu);
        
        else
            % CPU-only version
            iterate_T_k = zeros(TotalState, TotalState, num_t);
            
            for i = 1:num_dim
                for j = 1:num_dim
                    temp_mat_M1_conv_full = zeros(1, 2*num_t - 1);
                    for k = 1:num_dim
                        a = reshape(curr_T_k_iterate(i,k,:), 1, num_t);
                        b = TakeDerivative(reshape(mat_Gm(k,j,:), 1, num_t), StepSize);
                        temp_mat_M1_conv_full = temp_mat_M1_conv_full + conv(a, b, 'full') * StepSize;
                    end
                    temp_mat_M1_conv_same = temp_mat_M1_conv_full(1:num_t);
                    temp_mat_M1_conv_same(temp_mat_M1_conv_same < 0) = 0;
                    iterate_T_k(i,j,:) = temp_mat_M1_conv_same;
                end
            end
        
            next_mat_M_term_k = zeros(TotalState, TotalState, num_t);
            
            for i = 1:num_dim
                for j = 1:num_dim
                    a = reshape(I_mat_G(j,j,:), 1, num_t);
                    b = TakeDerivative(reshape(iterate_T_k(i,j,:), 1, num_t), StepSize);
                    temp_mat_M2_conv_full = conv(a, b, 'full') * StepSize;
                    temp_mat_M2_conv_same = temp_mat_M2_conv_full(1:num_t);
                    temp_mat_M2_conv_same(temp_mat_M2_conv_same < 0) = 0;
                    next_mat_M_term_k(i,j,:) = temp_mat_M2_conv_same;
                end
            end
        end


    end
    
   
end
    
end