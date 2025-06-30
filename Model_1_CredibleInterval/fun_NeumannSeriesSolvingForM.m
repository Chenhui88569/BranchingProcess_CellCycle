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
    if FlagTreatment == 1
        iterate_T_k =  zeros(TotalState,TotalState , num_t);
        for i = 1: num_dim
                for j = 1:num_dim
                    temp_mat_M1_conv_full  = 0;
                    for k = 1:num_dim
                        temp_mat_M1_conv_full =  temp_mat_M1_conv_full + conv( reshape( curr_T_k_iterate(i,k,:),1,num_t),   TakeDerivative(reshape(  mat_Gm(k,j,:),1,num_t),StepSize),'full')*StepSize;
                    end
                    temp_mat_M1_conv_same =  temp_mat_M1_conv_full(1:num_t);
                    temp_mat_M1_conv_same(temp_mat_M1_conv_same <0) = 0;
                    iterate_T_k(i,j,:) = temp_mat_M1_conv_same;
                    
                    
%                     temp_mat_M1_conv_full_death = conv( reshape( curr_T_k_iterate(i,i+ctr-1,:),1,num_t),   TakeDerivative(reshape(mat_Gm(i+ctr-1,end,:),1,num_t),StepSize),'full' )*StepSize;
%                     temp_mat_M1_conv_death_same =   temp_mat_M1_conv_full_death(1:num_t);
%                     temp_mat_M1_conv_death_same( temp_mat_M1_conv_death_same< 0)= 0;
%                     iterate_T_k(i,end,:) =  temp_mat_M1_conv_death_same;
                end
        end
        
       next_mat_M_term_k  =   zeros(TotalState,TotalState , num_t);
       for i = 1:num_dim
            for  j  = 1: num_dim
                temp_mat_M2_conv_full =  conv( reshape(I_mat_G(j,j,:),1,num_t),   TakeDerivative(reshape(  iterate_T_k(i,j,:),1,num_t) ,StepSize),'full' )*StepSize;
                % temp_mat_M2_conv_full =  conv( reshape(  iterate_T_k(i,j,:),1,num_t),  reshape(I_mat_G(j,j,:),1,num_t) ,'full', StepSize, )*StepSize;
                 temp_mat_M2_conv_same =   temp_mat_M2_conv_full(1:num_t);
                 temp_mat_M2_conv_same(temp_mat_M2_conv_same<0)  = 0;
%                 InterVar = 0 ;
%                 for m = 1: size(mat_m,1)
%                     temp_mat_M2_conv_full = conv( reshape(  iterate_T_k(i,m,:),1,num_t),  reshape(I_mat_G(m,j,:),1,num_t) ,'full', StepSize, )*StepSize;
%                     temp_mat_M1_conv_same =   temp_mat_M2_conv_full(1:num_t);
%                     InterVar   =   InterVar +     temp_mat_M1_conv_same;
%                 end
                next_mat_M_term_k(i,j,:) =  temp_mat_M2_conv_same ;
            end
        end
        
    end
    
   
end
    
end