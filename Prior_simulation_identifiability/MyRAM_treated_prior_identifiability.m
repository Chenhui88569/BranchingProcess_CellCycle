function  [Para_col_accepted,  acc_rate_col,    S_curr, Para_col_all,Diff_norm_col]  = MyRAM_treated_prior_identifiability(  T, num_para, Para_col_accepted , S_curr,myfun,Var_name, outputfile_dir, Model_Index, num_HistPoints,lb,ub, ...
    synthetic_data, para_change_exclude_idx,  para_fixed, outfile, sample, curr_T, Diff_norm_col, Para_col_all )



%% 
num_figure = ceil(num_para/8);
for i = 1:  num_figure
    figure(i);
end

%IncreaseGen_flag = 0;
% f_FalsePlot = figure(num_figure+1);
theta0 = Para_col_accepted(curr_T ,:);
model_obj = create_model(Model_Index, theta0, [] );
[~, Target_curr,  diff_norm_sum]   = myfun( theta0 , Model_Index ,outfile, model_obj, lb,ub,  synthetic_data, para_change_exclude_idx,  para_fixed  )  ;

for i = 1: length(Para_col_accepted(1,:))
    fprintf( outfile, 'lb(%d)= %f  , ub(%d)=%f \n ', i,  lb(i), i, ub(i));
end
% %Var_name  =  {'V_{max}' , 'Q_{21}' ,    'V_1',  'V_2',  'K_m','\sigma'};
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
[~,idx] = sort( [FigList.Number]);
FigList  = FigList(idx);
%FigList(end) = [];



alpha_bar = 0.234; %

acc_rate_col = zeros(1, T);
gamma= 0.5;

%Para_col_all = zeros(size(Para_col_accepted ));
%Diff_norm_col = zeros(T,1);
%outputfile_dir = 'MCMC_cellular/';
%% 
while  curr_T- num_HistPoints+1 <= T 
    local_model_obj = model_obj;
    rng('shuffle')
    fprintf(outfile,'curr_T = %d \n ', curr_T );
    curr_para = Para_col_accepted(curr_T,:);
    curr_para = curr_para';
    %FlagSteadyState = false;
    outercount_Prod = 0;
    outercount_NaN = 0;
    %intercount  = 0;
    ValidFlag = 0;
    
    while ~ValidFlag
         U = randn(num_para,1);
         next_para = curr_para +  S_curr*U;
         prior_val_prod = TestInMyRAM(next_para, lb,ub);
         while   any(next_para < 0)  ||  prior_val_prod  == 0 ||    ( sum(next_para(1:3))>1 &&   local_model_obj.TreatmentEffectOnPhases == 3 )  || ( sum(next_para([4,8]))>1 &&   local_model_obj.TreatmentEffectOnPhases == 3 )  ||   ( sum(next_para(2:3))>1 &&   local_model_obj.TreatmentEffectOnPhases == 2)   
                U = randn(num_para,1);
                next_para = curr_para +  S_curr*U;
                prior_val_prod = TestInMyRAM(next_para , lb,ub);
                outercount_Prod  =  outercount_Prod +1;
          end
          try 
                   [exception, Target_next,  diff_norm_sum]  = myfun( next_para, Model_Index, outfile ,local_model_obj , lb,ub, synthetic_data, para_change_exclude_idx,  para_fixed)  ;
                   if isempty(exception)
                        ValidFlag = 1;
                   end
          catch  e
                   outercount_NaN =  outercount_NaN+1;
                   %msgText = getReport(exception,'hyperlinks','off');
                   fprintf(outfile, e.message);
                   ValidFlag = 0;
                   continue

          end
      
    end
    Diff_norm_col( curr_T) = diff_norm_sum;
%     while  ~FlagSteadyState
%         IncreaseGen_flag = 0;
%         outercount  =   outercount  +1;
%         intercount =  intercount + 1;
%         if     intercount   == 1 
%             [Target_next, FlagSteadyState ]  = myfun( next_para, Model_Index, outfile, IncreaseGen_flag, lb,ub  )  ;
%         elseif  intercount  == 2
%             IncreaseGen_flag = 1;
%             [Target_next, FlagSteadyState ]  = myfun( next_para, Model_Index, outfile, IncreaseGen_flag, lb,ub  )  ;
%         elseif  intercount  > 2
%             U = randn(num_para,1);
%             next_para = curr_para +  S_curr*U;
%             while  ~prod( next_para > 0)
%                 U = randn(num_para,1);
%                 next_para = curr_para +  S_curr*abs(U);
%             end
%             IncreaseGen_flag  = 0;
%             [Target_next, FlagSteadyState ]  = myfun( next_para, Model_Index, outfile,IncreaseGen_flag    )  ;
%         end
%     end
    fprintf(outfile,' \t Find next_para in %d attempts(prod: %d, nan: %d) \n ',  outercount_Prod + outercount_NaN ,outercount_Prod ,outercount_NaN );

    Para_col_all(curr_T+1,:)  =  next_para';
    alpha_i =  min([1,   exp(Target_next  - Target_curr)       ]) ;
    u = rand;
    if  u <  alpha_i   % Target_next/Target_curr
        Para_col_accepted(curr_T+1,:)  =  next_para'; 
        acc_rate_col(curr_T - num_HistPoints+1) = 1;
        fprintf(outfile, '\t accept , %d \n', curr_T );
        Target_curr =     Target_next ;
    else
        Para_col_accepted(curr_T+1,:)  =   curr_para; 
        Target_curr = Target_curr;
        fprintf(outfile, '\t reject , %d \n', curr_T );
    end
    
    
    factor  =    U*U'./vecnorm( U,2)^2;

    temp_mat = S_curr*(  eye(num_para) + min( [1,curr_T^(-gamma)] )*(    alpha_i  - alpha_bar  )*factor  )*S_curr';
   
    S_curr = chol(temp_mat, 'lower');
    
    if rem(curr_T,3) == 0 %100*floor(curr_T/100) = curr_T
        save( fullfile( outputfile_dir ,  strcat('Para_treated_', num2str(Model_Index), 'sample', num2str(sample),   '.mat' ) ) , 'Para_col_accepted'  , 'acc_rate_col',  'S_curr' , 'Para_col_all','Diff_norm_col' ,'-v7.3');
        fprintf(outfile,'Saved .mat file %d \n ', curr_T);
    end
    
    %% plot trajectory for the parameters
    flag = 0 ;
    for iFig = 1:length(FigList)
        FigHandle = FigList(iFig);
        set(0, 'CurrentFigure', FigHandle);
        for j = 1:8
            subplot(2,4, j)
            idx = (iFig-1)*8+j;
            plot(Para_col_accepted(1: curr_T, idx),'k-','LineWidth',1.5)
            xlabel('iteration')
            %axis([1+T*fix(curr_T/T),T+T*fix(curr_T/T),0,15]);
            grid off
            title(Var_name(idx),  'Interpreter','Latex')
            drawnow
            if idx == num_para
                flag  = 1;
                break
            end
        end
        if flag == 1
            break
        end 
    end
    curr_T  = curr_T  +1 ;
end
%fclose(outfile);
for iFig = 1:length(FigList)
        FigHandle = FigList(iFig);
        exportgraphics(  FigHandle ,  fullfile(outputfile_dir, ['Traj_', num2str( iFig), '_sample', num2str(sample),'.png'] ) ,'Resolution' , 300   );
        savefig(  FigHandle ,  fullfile(outputfile_dir, ['Traj_', num2str( iFig), '_sample', num2str(sample), '.fig'] ));
end
  