classdef create_model_CombinationTreatment
    properties
        CellLine_Model_Index
        outputfile_dir
        clinical_time
        MaxT
        num_gen
        StepSize
        t_span
        t_span_AfterTreatment
        Para_untreated
        para_unknown
        num_t
        initial_States_convergence_RemainingTimeDist
        Dose_1thDrug
        Dose_2thDrug
        TreatmentEffectOnPhases_1thDrug
        TreatmentEffectOnPhases_2thDrug
        num_GenOfHealthycells 
        num_RroundsOfTreatedcells
        CotreatmentFlag 
        para_label
        RoundsFor1thDrug
        RoundsFor2thDrug
        
    end
    
    
    methods
        
        function obj = create_model_CombinationTreatment(CellLine_Model_Index, para_unknown, output_dir, Dose_1thDrug,Dose_2thDrug )
            %  output_dir is used to store the M matrix plot
            addpath('ReadYAML-master/ReadYAML-master/')
            
            ymlflie = ['DatasetsInfo'  num2str( CellLine_Model_Index )  '.yml'];
            datastr = ReadYaml(ymlflie);
            obj.CellLine_Model_Index =  CellLine_Model_Index;
            obj.num_gen = datastr.num_gen;
            %obj.Model_Index =  Model_Index;
            if  ~isempty(output_dir)
                obj.outputfile_dir =  output_dir;
            else
                obj.outputfile_dir =   ['Figs_', MCMC_sample_treated_DIR];
            end
            
            obj.clinical_time=    datastr.clinical_time;
            obj.MaxT = obj.clinical_time + 4 ;
            %obj.TreatmentEffectOnPhases = datastr.TreatmentEffectOnPhases ;
            obj.StepSize = datastr.StepSize;
            t_span_full = 0: obj.StepSize: obj.MaxT;
            if t_span_full(end) ~=  obj.MaxT
                t_span_full  = [t_span_full obj.MaxT];
            end
            obj.t_span = t_span_full ; % t_span before treatment
            
            obj.Para_untreated = datastr.Para_untreated;
            obj.num_t = length(t_span_full );
            obj.para_unknown =   para_unknown;
            data_treated_dir  = datastr.States_convergence_RemainingTimeDist;
            load(   data_treated_dir , 'States_convergence_RemainingTimeDist');
            obj.initial_States_convergence_RemainingTimeDist =   States_convergence_RemainingTimeDist;
                 
            
            ymlflie = 'Cotreatment.yml';
            datastr = ReadYaml(ymlflie);
            obj.t_span_AfterTreatment = 0: obj.StepSize: datastr.t_span_AfterTreatment_MatT;
            
            obj.TreatmentEffectOnPhases_1thDrug = datastr.TreatmentEffectOnPhases_1thDrug;
            obj.TreatmentEffectOnPhases_2thDrug = datastr.TreatmentEffectOnPhases_2thDrug;
            obj.num_RroundsOfTreatedcells = datastr.num_RroundsOfTreatedcells;
            obj.num_GenOfHealthycells =  datastr.num_GenOfHealthycells;
            obj.CotreatmentFlag  = datastr.CotreatmentFlag ;
            if  obj.CotreatmentFlag  == 1 
                 obj.RoundsFor1thDrug = 0;
                 obj.RoundsFor2thDrug = 0;
            else
                 obj.RoundsFor1thDrug = 2;
                 obj.RoundsFor2thDrug = 2;
            end
            obj.para_label=  ["q_{1,G2M}"
                "q_{2,G2M}"
                "q_{3,G2M}"
                "q_{4,G2M}"
                "m_{d,FR,G2M}^{max}"
                "m_{d,UR,G2M}^{max}"
                "m_{d,MS,G2M}^{max}"
                "m_{d,UR,G1}^{max}"
                "m_{d,UR,S}^{max}"
                "\lambda_{d,max,G1}"
                "\lambda_{G2arrest,UR, max}"
                "\lambda_{G2arrest,FR, max}"
                "\lambda_{G2arrest,MS, max}"
                "\lambda_{G1arrest, max}"
                "\lambda_{d,max,S}"
                "\lambda_{Sarrest, max}"
                "{b,G2M}"
                "{b}_{d,G2M}"
                "n_{G2M}"
                "EC_{50,d,G2M}"
                "EC_{50,G2M}"
                "$a_{S}$"
                "$b_{S}$"
                "q_{1,S}"
                "q_{2,S}"
                "q_{3,S}"
                "m_{d,FR,S}^{max}"
                "m_{d,UR,S}^{max}"
                "m_{d,UR,G2M}^{max}"
                "m_{d,G1block}^{max}"
                "\lambda_{d,max, G2M}"
                "\lambda_{Sarrest,UR,max}"
                "\lambda_{Sarrest,FR,max}"
                "\lambda_{G2Marrest,max}"
                "\lambda_{G1block, max}"
                "{b,S}"
                "{b}_{d,S}"
                "n_{S}"
                "EC_{50,d,S}"
                "EC_{50,S}"
                "$a_{G2M}$"
                "$b_{G2M}$"
                ]; 
                obj.Dose_1thDrug =  Dose_1thDrug ;
                obj.Dose_2thDrug = Dose_2thDrug ;
        end
    
        function  [a, b_scale, mat_G  ,  mat_Gm,  mat_m  , I_mat_G,   G1_pdf, S_pdf, G2M_pdf, G1_cdf, S_cdf, G2M_cdf ] =  mat_init(obj)
        
        %%  ChoiceOfDist
             ChoiceOfDist_untreated = 'Gamma';
        
            %% Time steps
            obj.StepSize= 0.1;
            a = obj.Para_untreated(1:3);
            b =  obj.Para_untreated(4);
            b_scale = 1/b;
            
            %% G1
            [G1_pdf, G1_cdf]  =  PdfCdfCal(ChoiceOfDist_untreated,  obj.t_span , a(1) ,b_scale);
            %% S
            [S_pdf, S_cdf]  =  PdfCdfCal(ChoiceOfDist_untreated,  obj.t_span , a(2) ,b_scale);
            %% G2
            [G2M_pdf, G2M_cdf]  =  PdfCdfCal(ChoiceOfDist_untreated,  obj.t_span , a(3) ,b_scale);
            %%
            
            m12 = 1 ;
            m23 = 1;
            m31 = 2;
            
            TotalState = 3*obj.num_gen;
            %% matrix m
            temp = diag( repmat([m12, m23, m31],1,obj.num_gen )  ,1);
            %mat_m= temp;
            mat_m= [temp(1: TotalState, 1: TotalState),   zeros(TotalState,1 )    ;
                zeros(1, TotalState+1 )  ];
            
            %% matrix G
            mat_G = zeros(TotalState+1, TotalState+1 ,obj.num_t );
            mat_G(1,1,:) = G1_cdf;
            mat_G(2,2,:) = S_cdf;
            mat_G(3,3,:) = G2M_cdf;
            %single dose
            for k = 2:obj.num_gen
                mat_G( 1+ 3*(k-1),  1+ 3*(k-1) ,:   ) = G1_cdf;
                mat_G( 2+ 3*(k-1),  2+ 3*(k-1),:   )  =  S_cdf;
                mat_G( 3+ 3*(k-1),  3+ 3*(k-1),:   )  =  G2M_cdf;
            end
            %mat_G(end,  end,:   )  =   G1_cdf;
            
            %% matrix I-G
            mat_I = reshape( repmat( eye(TotalState+1), 1,  obj.num_t) ,  TotalState+1, TotalState+1 , obj.num_t );
            I_mat_G =  mat_I - mat_G ;
            %% matrix Gm
            mat_Gm = zeros(TotalState+1,TotalState+1, obj.num_t);
            for i = 1: size(mat_m,1)
                if  i < TotalState-1
                    mat_Gm(i,i+1 ,: ) = reshape( mat_G( i,  i,:   ) , 1, obj.num_t  ) * mat_m(i,i+1); %mat_m(i,j) scalar value
                    %mat_Gm(i, end ,: ) = reshape( mat_G( i,  i,:   ) , 1, num_t  ) * mat_m(i,end); %mat_m(i,j) scalar value
                end
                
            end
        end
        %%
        
        function   full_mat_M_original =  M_beforetreatement(obj)
            %% Calulate the M matrix before the treatment
            ctr = 0; % the number of iterates
            [ ~,~, ~,  mat_Gm,  ~  , I_mat_G,   ~, ~, ~, ~, ~, ~ ] =  obj.mat_init();
            TotalState = 3*obj.num_gen+1;
            while ctr<=  TotalState
                if ctr  == 0
                    % k=0
                    next_mat_M = I_mat_G;
                elseif ctr == 1    % k = 1
                    curr_mat_M = next_mat_M;
                    next_T_k_iterate = mat_Gm;
                    
                    [ next_mat_M_term_k, ~ ]  =  fun_NeumannSeriesSolvingForM(mat_Gm, TotalState, [],  I_mat_G, ctr, obj.StepSize,0);
                    next_mat_M = curr_mat_M+next_mat_M_term_k;
                else  %k >= 2
                    curr_mat_M = next_mat_M;
                    curr_T_k_iterate = next_T_k_iterate;
                    [ next_mat_M_term_k, iterate_T_k ]  =  fun_NeumannSeriesSolvingForM(mat_Gm, TotalState, curr_T_k_iterate,  I_mat_G, ctr, obj.StepSize,0);
                    next_T_k_iterate  = iterate_T_k ; %update iterate
                    next_mat_M  = curr_mat_M + next_mat_M_term_k;
                end
                ctr = ctr +1;
            end
            %%
            mat_M =  next_mat_M;
            
            %% scale idx to time zero
            full_mat_M_original = mat_M;
            
            
        end
        %%
        function  [mat_M_convergence_BeforeTreatment_complete,   mat_M_G1_BeforeTreatment, States_convergence,   M_aftertreatment_col] = M_aftertreatment(obj)
            [ a, b_scale, mat_G  ,  mat_Gm,  mat_m  , I_mat_G,   G1_pdf, S_pdf, G2M_pdf, G1_cdf, S_cdf, G2M_cdf ] =  obj.mat_init();
            idx = find(obj.t_span == obj.clinical_time);  % Find where the trajectory converges
            cutoff =  idx;
            
            
            TotalState = 3*obj.num_gen;
            Gen_G1_Nonzero_col = cell(TotalState, 2);
            
            
            %         OriginalGm = mat_Gm;
            %         OriginalG = mat_G;
            %         OriginalM =  mat_m;
            %         Original_I_mat_G = I_mat_G ;
            %
            %         ChoiceOfDist_treated_Dd = 'Weibull' ;
            %         ChoiceOfDist_treated_D = 'Weibull';% 'Loglogistic';
            %         ChoiceOfDist_treated_U = 'Exponential';
            
            
            iterate_T_k_col = cell(1, 3*obj.num_gen,1);
            mat_M_col = cell(1, 3*obj.num_gen,1  );
            flag_G1  = 0;
            Threshold_DeterminingGenAtTreatment = 2;
            ctr = 0;
            while ctr<=  TotalState
                if ctr  == 0
                    % k=0
                    next_mat_M = I_mat_G;
                elseif ctr == 1    % k = 1
                    curr_mat_M = next_mat_M;
                    next_T_k_iterate = mat_Gm;
                    
                    [ next_mat_M_term_k, ~ ]  =  fun_NeumannSeriesSolvingForM(mat_Gm, TotalState+1, [],  I_mat_G, ctr, obj.StepSize,0);
                    next_mat_M = curr_mat_M+next_mat_M_term_k;
                    Gen_G1_Nonzero_col{ctr,1} = find( next_mat_M(1,:, cutoff ) > Threshold_DeterminingGenAtTreatment );
                    Gen_G1_Nonzero_col{ctr,2} = next_mat_M(1,:, cutoff );
                    
                    iterate_T_k_col{1,ctr+1} = mat_Gm;
                    mat_M_col{1, ctr+1} = next_mat_M;
                else  %k >= 2
                    if flag_G1 == 0 %|| flag_S ==0 || flag_G2M ==0
                        %% The treatment has not hitted.
                        % This should come before flag_G1, flag_S, flag_G2M turn 1
                        curr_mat_M = next_mat_M;
                        curr_T_k_iterate = next_T_k_iterate;
                        [ next_mat_M_term_k, iterate_T_k ]  =  fun_NeumannSeriesSolvingForM(mat_Gm,  TotalState+1, curr_T_k_iterate,  I_mat_G, ctr, obj.StepSize,0);
                        next_T_k_iterate  = iterate_T_k ; %update iterate
                        next_mat_M  = curr_mat_M + next_mat_M_term_k;
                        
                        % if flag_G1 == 0
                        iterate_T_k_col{1,ctr+1} =  next_T_k_iterate;
                        mat_M_col{1, ctr+1} = next_mat_M;
                        
                    end
                    %% find which states are progressing at the time of convergence for G1, S, G2M phase
                    if flag_G1 == 0% || flag_S ==0 || flag_G2M ==0
                        
                        Gen_G1_Nonzero_col{ctr,1} = find( next_mat_M(1,:, cutoff ) > Threshold_DeterminingGenAtTreatment );
                        Gen_G1_Nonzero_col{ctr,2} = next_mat_M(1, : , cutoff );
                        
                        %% % statisfy the condition of convergence
                        if  flag_G1  == 0 &&   isequal( Gen_G1_Nonzero_col{ctr,1} , Gen_G1_Nonzero_col{ctr-1,1}) && ~isempty( Gen_G1_Nonzero_col{ctr,1} )  && ~isempty( Gen_G1_Nonzero_col{ctr-1,1} )
                            flag_G1 = 1;
                            temp_G1_convergence =   Gen_G1_Nonzero_col{ctr,1};
                            %                         locs_col = zeros(1, length(temp_G1_convergence ));
                            %                         for i = 1: length(temp_G1_convergence )
                            %                             temp_vec = reshape( next_mat_M(1, temp_G1_convergence(i) , : ),1, obj.num_t   );
                            %                             [pks,locs]  = findpeaks(  temp_vec , 'NPeaks',1) ;
                            %                             locs_col(i) = locs;
                            %                         end
                            %  calculate only once
                            % StartState_G1 =  temp_G1_convergence(end)+1;
                            %StateAtTreatment =  temp_G1_convergence(end) ;
                            States_convergence  = temp_G1_convergence;
                            break
                        end
                    end
                end
                ctr = ctr +1;
            end
            
            mat_M_convergence_BeforeTreatment_complete =   reshape( next_mat_M(1, States_convergence , 1:end), length( States_convergence ) , obj.num_t);
            mat_M_G1_BeforeTreatment =  next_mat_M; % change me!
            
            %%
            num_CyclingState = length( States_convergence );
            %         initial_States_convergence_RemainingTimeDist = zeros( num_CyclingState ,  obj.num_t  );
            M_aftertreatment_col =  cell(num_CyclingState,1 ); 
            
            if obj.CotreatmentFlag  == 1 % Cotreatment
                mat_m =  Cotreatment_BuildMatrixm_AfterTreatment(obj);
                TotalState = size(mat_m ,1);
                mat_G = Cotreatment_BuildMatrixG_AfterTreatment(obj,G1_pdf, S_pdf, G2M_pdf, G1_cdf,G2M_cdf,S_cdf,TotalState);

                
                for i = 1:num_CyclingState
                    CyclingSate =  States_convergence(i);
                    phase = mod( CyclingSate ,3 );
                    %%
                    g_est2_cdf =  obj.initial_States_convergence_RemainingTimeDist(i,:);
                    g_est2_cdf = adjust_vector_length(g_est2_cdf, obj);
                    g_est2_cdf(g_est2_cdf>0.95) = 1; 
                    mat_M =  Cotreatment_AfterTreatmentBranchingProcess(obj, mat_m ,  mat_G  ,   g_est2_cdf, phase);
                    M_aftertreatment_col{i} = mat_M(1,:,:);
                end
                
            end
          %  if obj.CotreatmentFlag  == 0 % Cotreatment
                
                
%                 
%                 if obj.TreatmentEffectOnPhases_1thDrug == 3 % G2M
%                      
%                      
%                 elseif obj.TreatmentEffectOnPhases_1thDrug == 2 % G2M 
%                     
%                     
%                 end
%                 
%                 
%                 if obj.TreatmentEffectOnPhases_1thDrug == 2 % G2M 
%                     
%                      
%                 elseif obj.TreatmentEffectOnPhases_1thDrug == 3 % G2M 
%                     
%                 end
                    
                    
                %num_GenOfHealthycells = 2 ;% one round of healty cells is mandatory
%                 lambda_d_max_G1_G2M = obj.para_unknown(8);
%                 lambda_Sarrest_UR_max = obj.para_unknown(9);
%                 lambda_Sarrest_FR_max = obj.para_unknown(10);
%                 lambda_G2Marrest_max = obj.para_unknown(11);
%                 lambda_G1block_max = obj.para_unknown(12);
%                 tilde_b =  obj.para_unknown(13);
%                 tilde_b_d =  obj.para_unknown(14);
%                 n_S =  obj.para_unknown(15);
%                 EC50_d =  obj.para_unknown(16);
%                 EC50 =  obj.para_unknown(17);
%                 %Dose = 20*EC50_d;
%                 mat_m = Sphasedrug_BuildMatrixm_AfterTreatment( obj.Dose, obj, obj.num_GenOfHealthycells, obj.num_RroundsOfTreatedcells);
%                 
%                 %% create S arrest
%                 lambda_fun = @(lambda_max) lambda_max * ((obj.Dose / EC50)^n_S) / (1 + (obj.Dose / EC50)^n_S);
%                 
%                 %% create G1 block
%                 lambda_G1block  = lambda_fun(lambda_G1block_max );
%                 [~, ~, G1block ]  =   PdfcdfCal_treated('Weibull',  obj.t_span ,  lambda_G1block  , tilde_b  ); %A: scale, B , shape
%                           
%                  %% create S arrest
%                 lambda_Sarrest_UR  = lambda_fun(lambda_Sarrest_UR_max);
%                 [~, ~, S_arrest_UR ]  =   PdfcdfCal_treated('Weibull',  obj.t_span , lambda_Sarrest_UR  , tilde_b  ); %A: scale, B , shape
%                 lambda_Sarrest_FR  = lambda_fun(lambda_Sarrest_FR_max);
%                 [~, ~, S_arrest_FR ]  =   PdfcdfCal_treated('Weibull',  obj.t_span , lambda_Sarrest_FR  , tilde_b  ); %A: scale, B , shape
%                 %% Create tilde_G2M
%                 lambda_G2Marrest = lambda_fun(lambda_G2Marrest_max);
%                 lambda_d_max_G1_G2M =  lambda_fun(lambda_d_max_G1_G2M);
%                 [ m_D_G2Marrest ,D_G2Marrest_pdf , ~ ]  =   PdfcdfCal_treated('Weibull',  obj.t_span ,   lambda_G2Marrest   , tilde_b  );
%                 D_G2Marrest_pdf(D_G2Marrest_pdf  == inf) = 0;
%                 [ m_Dd , Dd_pdf , ~ ]  =   PdfcdfCal_treated('Weibull',  obj.t_span ,  lambda_d_max_G1_G2M  , tilde_b_d  );
%                 Dd_pdf( Dd_pdf == inf) = 0;
%                 ChoiceOfDist_treated_U = 'Exponential';
%                 [m_U, U_pdf, ~]  =   PdfcdfCal_treated(ChoiceOfDist_treated_U, obj.t_span , 1 ,[] );
%                 U_pdf( U_pdf == inf) = 0;
%                 tilde_G2M  = fun_treated_PDFWatingTime( m_U, m_Dd , m_D_G2Marrest, G1_pdf , U_pdf, Dd_pdf, D_G2Marrest_pdf ,  obj.num_t, obj.t_span , obj.StepSize   )  ;
%                 for i = 1:num_CyclingState
%                     CyclingSate =  States_convergence(i);
%                     phase = mod( CyclingSate ,3 );
%                     %%
%                     g_est2_cdf =  obj.initial_States_convergence_RemainingTimeDist(i,:);
%                     g_est2_cdf = adjust_vector_length(g_est2_cdf, obj);
%                     if  phase == 0
%                         mat_M =  Sphasedrug_AfterTreatmentBranchingProcess_G2M_2(G1block,   mat_m , S_arrest_UR,  S_arrest_FR,  tilde_G2M, g_est2_cdf, G1_cdf, S_cdf, G2M_cdf,  obj.num_GenOfHealthycells ,obj.StepSize,obj.t_span_AfterTreatment , obj.num_RroundsOfTreatedcells);
%                     elseif phase == 1
%                         mat_M =  Sphasedrug_AfterTreatmentBranchingProcess_G1_2( G1block, mat_m , S_arrest_UR,  S_arrest_FR,   tilde_G2M, g_est2_cdf, G1_cdf, S_cdf, G2M_cdf,  obj.num_GenOfHealthycells ,obj.StepSize,obj.t_span_AfterTreatment, obj.num_RroundsOfTreatedcells);
%                     elseif phase == 2
%                         mat_M =  Sphasedrug_AfterTreatmentBranchingProcess_S_2( G1block ,  mat_m , S_arrest_UR,  S_arrest_FR,  tilde_G2M, g_est2_cdf, G1_cdf, S_cdf, G2M_cdf,  obj.num_GenOfHealthycells ,obj.StepSize,obj.t_span_AfterTreatment, obj.num_RroundsOfTreatedcells);
%                     end
%                     M_aftertreatment_col{i} = mat_M(1,:,:);
%                 end
%                 
%             end
            
        end
     %%
     
     
     function [CellFractions, N_d, t_span_new] =  M_matrix_plot(obj)
         %accross generation
        % full_mat_M_original =  obj.M_beforetreatement();
         [ CellFractions, N_d, t_span_new, mat_M_convergence_BeforeTreatment_complete,   mat_M_G1_BeforeTreatment,   M_aftertreatment_col, States_convergence   ]  = model_simulation(obj );
         % cells in different generations after the treatment 
         fun_reshape = @(x,n) reshape(x, n, length(obj.t_span_AfterTreatment) );
         if  obj.CotreatmentFlag  == 1
             TotalGens = obj.num_RroundsOfTreatedcells+obj.num_GenOfHealthycells; 
             NumberofCellsInDifferentGensAfterTreatment = zeros( 3*TotalGens, length(obj.t_span_AfterTreatment)  );
             
             for  k  = 1: length(States_convergence)
                 CurrentState =  States_convergence(k);
                 AbsolutePhase = mod(  CurrentState ,3);
                 M_aftertreatment_CurrentState  =    M_aftertreatment_col{k};
                 TotalState = size(M_aftertreatment_CurrentState, 2 );
                 [a_G1  , a_S , a_G2M ,~ ,~, ~  ] = Cotreatment_CollectIdx_CrossGens(AbsolutePhase, obj);
                 
                 for i = 1:  TotalGens
                     NumberofCellsInDifferentGensAfterTreatment(3*(i-1) +1,:) =  NumberofCellsInDifferentGensAfterTreatment(3*(i-1) +1,:) + sum(fun_reshape(M_aftertreatment_CurrentState(1, a_G1{i},: ),  length(a_G1{i})      ),1) ;
                     NumberofCellsInDifferentGensAfterTreatment(3*(i-1) +2,:) =  NumberofCellsInDifferentGensAfterTreatment(3*(i-1) +2 ,:) + sum(fun_reshape(M_aftertreatment_CurrentState(1, a_S{i},: ),   length(a_S{i} )     ),1) ;
                     NumberofCellsInDifferentGensAfterTreatment(3*(i-1) +3,:) =  NumberofCellsInDifferentGensAfterTreatment(3*(i-1) +3 ,:) + sum(fun_reshape(M_aftertreatment_CurrentState(1, a_G2M{i},: ),  length( a_G2M{i}  )    ),1) ;
                     
                 end
             end
             f1 = figure;
             numGenerations =  TotalGens;
             % Generate a colormap for different generations
             colors = lines(numGenerations);
             lineStyles = {'-', ':', '--'};
             % Initialize the legend entries
             legendEntries = {};
             for gen = 1: numGenerations
                 for cellType = 1:3  % 1 for G1, 2 for S, 3 for G2M
                     
                     % Sample data (replace this with your actual data)
                     xData = t_span_new ;
                     yData = NumberofCellsInDifferentGensAfterTreatment((gen-1)*3+ cellType,:  ) ;
                     
                     % Plot the line with the appropriate style and color
                     plot(xData, yData, 'LineStyle', lineStyles{cellType}, 'Color', colors(gen,:), 'LineWidth',1.5);
                     hold on;
                     % Add entry to legend
                     cellTypeStr = {'G1', 'S', 'G2M'};
                     legendEntries{end + 1} = [cellTypeStr{cellType} ' - Gen ' num2str(gen)];
                     
                 end
                 
             end
             legend(legendEntries);
             xlabel('Time(hr)');
             ylabel('Y-axis label');
         end
         
         %% sequential treatment
         if  obj.CotreatmentFlag ==  0
             
             TotalGens =   4;
             NumberofCellsInDifferentGensAfterTreatment = zeros( 3* TotalGens, length(obj.t_span_AfterTreatment)    ); % The extra one state is for the additional G2M phase
              for  k  = 1: length(States_convergence)
                 CurrentState =  States_convergence(k);
                 AbsolutePhase = mod(  CurrentState ,3);
                 M_aftertreatment_CurrentState  =    M_aftertreatment_col{k};
                  
                 [a_G1 , a_S, a_G2M, ~ , ~,  ~  ] = SequentialTreatment_CollectIdx_CrossGens(AbsolutePhase, obj);
               
                 for i = 1:  TotalGens
                     if ~isempty( a_G1{i})
                         NumberofCellsInDifferentGensAfterTreatment(3*(i-1) +1,:) =  NumberofCellsInDifferentGensAfterTreatment(3*(i-1) +1,:) + sum(fun_reshape(M_aftertreatment_CurrentState(1, a_G1{i},: ),  length(a_G1{i})      ),1) ;
                     end
                     if ~isempty( a_S{i})
                         NumberofCellsInDifferentGensAfterTreatment(3*(i-1) +2,:) =  NumberofCellsInDifferentGensAfterTreatment(3*(i-1) +2 ,:) + sum(fun_reshape(M_aftertreatment_CurrentState(1, a_S{i},: ),   length(a_S{i} )     ),1) ;
                     end
                     
                     if ~isempty( a_G2M{i})
                         NumberofCellsInDifferentGensAfterTreatment(3*(i-1) +3,:) =  NumberofCellsInDifferentGensAfterTreatment(3*(i-1) +3 ,:) + sum(fun_reshape(M_aftertreatment_CurrentState(1, a_G2M{i},: ),  length( a_G2M{i}  )    ),1) ;
                     end
                 end
                 
              end
              f2 = figure;
              numGenerations = TotalGens ;
              % Generate a colormap for different generations
              colors = lines(numGenerations);
              lineStyles = {'-', ':', '--'};
              % Initialize the legend entries
              legendEntries = {};
              xData = t_span_new ;
              yData = NumberofCellsInDifferentGensAfterTreatment(1,:   ) ;
              
              % Plot the line with the appropriate style and color
              plot(xData, yData, 'LineStyle', lineStyles{3}, 'Color', colors(1,:), 'LineWidth',1.5);
              hold on;
              legendEntries{end + 1} = ['G2M- Gen ' num2str(1)];
              for gen = 2: numGenerations
                  for cellType = 1:3  % 1 for G1, 2 for S, 3 for G2M
                      
                      % Sample data (replace this with your actual data)
                      xData = t_span_new ;
                      yData = NumberofCellsInDifferentGensAfterTreatment((gen-1)*3+ cellType,:  ) ;
                      
                      % Plot the line with the appropriate style and color
                      plot(xData, yData, 'LineStyle', lineStyles{cellType}, 'Color', colors(gen,:), 'LineWidth',1.5);
                      hold on;
                      % Add entry to legend
                      cellTypeStr = {'G1', 'S', 'G2M'};
                      legendEntries{end + 1} = [cellTypeStr{cellType} ' - Gen ' num2str(gen+1)];       
                  end
                  
              end
              legend(legendEntries, 'FontSize', 15);
              xlabel('Time(hr)');
              ylabel('Y-axis label');
         end
         
         
         %%
         colorMatrix = [[0 0.4470 0.7410];[0.8500 0.3250 0.0980];[0.4940 0.1840 0.5560]];
         f4 = figure;
         hold on;
         for i = 1:3
             plot(t_span_new, CellFractions(i,:), 'Color', colorMatrix(i,:) ,'LineWidth',2);
         end
         legend({'G1', 'S', 'G2M'}, 'Location', 'northwest', 'FontSize', 15)
         xlabel('Hours')
         exportgraphics( f4 ,  fullfile(obj.outputfile_dir,'Sim.png' ) ,'Resolution' , 300   );
         
         
%          f5 = figure;
%          hold on;
%          for i = 1:3
%              plot(t_span_new(1:400), CellFractions(i,1:400), 'Color', colorMatrix(i,:) ,'LineWidth',2);
%          end
%          hold on;
%          legend({'G1', 'S', 'G2M'}, 'Location', 'best', 'FontSize', 15)
%          for i = 1:3
%              s = scatter(obj.data_time, obj.current_data(:,i), 80, 'filled','d');
%              s.MarkerEdgeColor = colorMatrix(i,:);
%              s.MarkerFaceColor = colorMatrix(i,:);
%          end
%          legend({'G1', 'S', 'G2M', 'data G1', 'data S', 'data G2M'}, 'Location', 'best', 'FontSize', 15)
%          xlabel('Hours')
%          exportgraphics( f5,  fullfile(obj.outputfile_dir,'Sim_short.png' ) ,'Resolution' , 300   );
        
         
         f6 = figure;
         plot(t_span_new, N_d)
         xlabel('Hours')
         ylabel('Cell numbers')
         exportgraphics( f6 ,  fullfile(obj.outputfile_dir,'Sim_d.png' ) ,'Resolution' , 300   );
         
         

     end
     
     %%
     function   [ CellFractions, N_d, t_span_new, mat_M_convergence_BeforeTreatment_complete,   mat_M_G1_BeforeTreatment,   M_aftertreatment_col, States_convergence   ]  = model_simulation(obj )
         %[M_afterTreatment_complete,   mat_M_G1,  next_mat_M_G1, States_convergence ] = obj.M_aftertreatment( );
         [mat_M_convergence_BeforeTreatment_complete,   mat_M_G1_BeforeTreatment, States_convergence,   M_aftertreatment_col] = obj.M_aftertreatment( );
         
         
         idx = find(obj.t_span == obj.clinical_time);  % Find where the trajectory converges
         
         N_col_G1_AcrossGen = 0;
         N_col_S_AcrossGen = 0;
         N_col_G2M_AcrossGen = 0;
         N_d = 0;
         t_span_new = obj.t_span_AfterTreatment;
         num_t_new = length(t_span_new);
         
         N_0_1 =   1.0845e+12 *0.3847;  %1.0845e+12 *   0.4126; % initial number of cells
         N_0_2 =     0;                         %1.0845e+12 * 0.3847 ;
         N_0_3 =     0;                      %1.0845e+12*     0.2028;
         N_0  =  [N_0_1; N_0_2 ; N_0_3];
         
         
         fun_reshape = @(x,n) reshape(x, n, num_t_new);
         if obj.CotreatmentFlag == 1 %G2M phase drug
             for k = 1: length(States_convergence)
                 num_cell_AtTreatment = mat_M_convergence_BeforeTreatment_complete(k,  idx);
                 CurrentState =  States_convergence(k);
                 AbsolutePhase = mod(  CurrentState ,3);
                 M_aftertreatment_CurrentState  =    M_aftertreatment_col{k};
                 %TotalState = size(M_aftertreatment_CurrentState, 2 );
                 
                 [~,~,~, G1_idx_col , S_idx_col,  G2M_idx_col  ]   = Cotreatment_CollectIdx_CrossGens(AbsolutePhase, obj);
                 
                 M_aftertreatment_S = fun_reshape(  M_aftertreatment_CurrentState(1,  S_idx_col , : ), length( S_idx_col  )  );
                 M_aftertreatment_G2M = fun_reshape(  M_aftertreatment_CurrentState(1,  G2M_idx_col , : ), length( G2M_idx_col  )  );
                 M_aftertreatment_G1 = fun_reshape(  M_aftertreatment_CurrentState(1,  G1_idx_col , : ), length( G1_idx_col  )  );
                 N_col_S_AcrossGen = N_col_S_AcrossGen +   num_cell_AtTreatment* sum(M_aftertreatment_S,1 ) ;
                 N_col_G2M_AcrossGen = N_col_G2M_AcrossGen +   num_cell_AtTreatment* sum(M_aftertreatment_G2M,1) ;
                 N_col_G1_AcrossGen = N_col_G1_AcrossGen +   num_cell_AtTreatment * sum( M_aftertreatment_G1,1) ;
                 N_d =  N_d + num_cell_AtTreatment* reshape(M_aftertreatment_CurrentState(1,end,:), 1, num_t_new);  %+  N_0(2)*reshape(mat_M_S(2,end,:), 1,num_t_new) + N_0(3)*reshape(mat_M_G2M(3,end,:), 1,num_t_new);
                
             end
         end
         
         if  obj.CotreatmentFlag == 0 %S phase drug
             for k = 1: length(States_convergence)
                 num_cell_AtTreatment = mat_M_convergence_BeforeTreatment_complete(k,  idx);
                 CurrentState =  States_convergence(k);
                 AbsolutePhase = mod(  CurrentState ,3);
                 M_aftertreatment_CurrentState  =    M_aftertreatment_col{k};
                 TotalState = size(M_aftertreatment_CurrentState, 2 );
                 
                 [~ ,~, ~, G1_idx_col , S_idx_col,  G2M_idx_col  ] = SequentialTreatment_CollectIdx_CrossGens(AbsolutePhase, obj);
                 
                 M_aftertreatment_S = fun_reshape(  M_aftertreatment_CurrentState(1,  S_idx_col , : ), length( S_idx_col  )  );
                 M_aftertreatment_G2M = fun_reshape(  M_aftertreatment_CurrentState(1,  G2M_idx_col , : ), length( G2M_idx_col  )  );
                 M_aftertreatment_G1 = fun_reshape(  M_aftertreatment_CurrentState(1,  G1_idx_col , : ), length( G1_idx_col  )  );
                 N_col_S_AcrossGen = N_col_S_AcrossGen +   num_cell_AtTreatment*  sum(M_aftertreatment_S,1 ) ;
                 N_col_G2M_AcrossGen = N_col_G2M_AcrossGen +   num_cell_AtTreatment* sum(M_aftertreatment_G2M,1) ;
                 N_col_G1_AcrossGen = N_col_G1_AcrossGen +   num_cell_AtTreatment * sum( M_aftertreatment_G1,1) ;
                 N_d =  N_d + reshape(M_aftertreatment_CurrentState(1,end,:), 1, num_t_new);  %+  N_0(2)*reshape(mat_M_S(2,end,:), 1,num_t_new) + N_0(3)*reshape(mat_M_G2M(3,end,:), 1,num_t_new);
                 
             end
         end
         
         N_col_DiffPhases =  [ N_col_G1_AcrossGen  ;
             N_col_S_AcrossGen;
             N_col_G2M_AcrossGen];
         
         
         N_total = sum( N_col_DiffPhases,1  ) ;
         CellFractions = N_col_DiffPhases./repmat(N_total,3,1);
         
         
  
     end


  end
    
    
end
