classdef create_model
    properties
        Model_Index
        outputfile_dir
        current_data
        clinical_time
        data_time
        MaxT
        TreatmentEffectOnPhases
        num_gen
        StepSize
        t_span
        t_span_AfterTreatment
        Para_untreated
        para_unknown
        num_t
        States_convergence_RemainingTimeDist
        Dose
        num_GenOfHealthycells 
        num_RroundsOfTreatedcells
    end
    
    
    methods
        
        function obj = create_model(Model_Index, para_unknown, output_dir)
            %  output_dir is used to store the M matrix plot
            addpath('ReadYAML-master/ReadYAML-master/')
            ymlflie = ['DatasetsInfo'  num2str( Model_Index )  '.yml'];
            datastr = ReadYaml(ymlflie);
            data_treated_dir  = datastr.data_treated_dir;
            load(data_treated_dir, 'DataTime', 'Current_Data');
            MCMC_sample_treated_DIR = datastr.MCMC_sample_treated_DIR;
            
            obj.Model_Index =  Model_Index;
            if  ~isempty(output_dir)
                obj.outputfile_dir =  output_dir;
            else
                obj.outputfile_dir =   ['Figs_', MCMC_sample_treated_DIR];
            end
            
            obj.clinical_time=    datastr.clinical_time;
            obj.current_data=  Current_Data;
            obj.data_time = DataTime;
            obj.MaxT = obj.clinical_time +   DataTime(end) + 4 ;
            obj.TreatmentEffectOnPhases = datastr.TreatmentEffectOnPhases ;
            obj.num_gen = datastr.num_gen;
            obj.StepSize = datastr.StepSize;
            t_span_full = 0: obj.StepSize: obj.MaxT;
            if t_span_full(end) ~=  obj.MaxT
                t_span_full  = [t_span_full obj.MaxT];
            end
            obj.t_span = t_span_full ;
            obj.t_span_AfterTreatment = 0: obj.StepSize: DataTime(end) + 4 ;
            obj.Para_untreated = datastr.Para_untreated;
            obj.num_t = length(t_span_full );
            obj.para_unknown =   para_unknown;
            data_treated_dir  = datastr.States_convergence_RemainingTimeDist;
            load(   data_treated_dir , 'States_convergence_RemainingTimeDist');
            obj.States_convergence_RemainingTimeDist =   States_convergence_RemainingTimeDist;
            obj.Dose = datastr.Dose;
            if   obj.data_time(end)<48
                obj.num_RroundsOfTreatedcells=  2;
                obj.num_GenOfHealthycells = 1; % should be larger than 1
            elseif  obj.data_time(end) > 48
                    obj.num_RroundsOfTreatedcells= 5;
                    obj.num_GenOfHealthycells =  1;
            end
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
                            AbsolutePhase = mod(temp_G1_convergence,3);
                            if obj.TreatmentEffectOnPhases == 3 %G2M phase drug
                                idx_CyclingG2M = find(AbsolutePhase == 0);  % relative indexes
                                % leftout_Cycling =  idx_CyclingG2M( locs_col(idx_CyclingG2M)< idx  );
                                % LastTreatedState = temp_G1_convergence(idx_CyclingG2M(end)); % absolute indexes
                            elseif  obj.TreatmentEffectOnPhases == 2 % S phase drug
                                idx_CyclingG1S = find(AbsolutePhase == 2|AbsolutePhase == 1 );  % relative indexes
                                % leftout_Cycling =  find( locs_col< idx ) ;
                                % LastTreatedState = temp_G1_convergence(idx_CyclingG1S(end)); % absolute indexes
                            end
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
            %         States_convergence_RemainingTimeDist = zeros( num_CyclingState ,  obj.num_t  );
            M_aftertreatment_col =  cell(num_CyclingState,1 );
            if obj.TreatmentEffectOnPhases == 3 % G2M phase cells
                    lambda_d_max_G1 = obj.para_unknown(10);
                    lambda_G2arrest_UR_max = obj.para_unknown(11);
                    lambda_G2arrest_FR_max = obj.para_unknown(12);
                    lambda_G2arrest_MS_max = obj.para_unknown(13);
                    lambda_G1arrest_max = obj.para_unknown(14);
                    lambda_d_max_S = obj.para_unknown(15);
                    lambda_Sarrest_max = obj.para_unknown(16);
                    
                    tilde_b =  obj.para_unknown(17);
                    tilde_b_d =  obj.para_unknown(18);
                    n_G2M =  obj.para_unknown(19);
                    EC50_d =  obj.para_unknown(20);
                    EC50 =  obj.para_unknown(21);
                    %Dose = 20*EC50_d;
                    mat_m = G2Mphasedrug_BuildMatrixm_AfterTreatment( obj.Dose, obj, obj.num_GenOfHealthycells, obj.num_RroundsOfTreatedcells);
                    
                    
                    
                    %% create G2 arrest
                    lambda_fun = @(lambda_max) lambda_max * ((obj.Dose / EC50)^n_G2M) / (1 + (obj.Dose / EC50)^n_G2M);
                    
                    lambda_G2Marrest_UR  = lambda_fun(lambda_G2arrest_UR_max);
                    [~, ~, G2_arrest_UR ]  =   PdfcdfCal_treated('Weibull',  obj.t_span , lambda_G2Marrest_UR  , tilde_b  ); %A: scale, B , shape
                    lambda_G2Marrest_FR  = lambda_fun(lambda_G2arrest_FR_max);
                    [~, ~, G2_arrest_FR ]  =   PdfcdfCal_treated('Weibull',  obj.t_span , lambda_G2Marrest_FR  , tilde_b  ); %A: scale, B , shape
                    lambda_G2Marrest_MS = lambda_fun(lambda_G2arrest_MS_max);
                    [~, ~, G2_arrest_MS ]  =   PdfcdfCal_treated('Weibull',  obj.t_span , lambda_G2Marrest_MS  , tilde_b  ); %A: scale, B , shape
                    %% Create tilde_G1
                    lambda_G1arrest = lambda_fun(lambda_G1arrest_max);
                    lambda_d_max_G1 =  lambda_fun(lambda_d_max_G1);
                    [ m_D_G1arrest ,D_G1arrest_pdf , ~ ]  =   PdfcdfCal_treated('Weibull',  obj.t_span ,  lambda_G1arrest  , tilde_b  );
                    D_G1arrest_pdf( D_G1arrest_pdf == inf) = 0;
                    [ m_Dd , Dd_pdf , ~ ]  =   PdfcdfCal_treated('Weibull',  obj.t_span ,  lambda_d_max_G1  , tilde_b_d  );
                    Dd_pdf( Dd_pdf == inf) = 0;
                    ChoiceOfDist_treated_U = 'Exponential';
                    [m_U, U_pdf, ~]  =   PdfcdfCal_treated(ChoiceOfDist_treated_U, obj.t_span , 1 ,[] );
                    U_pdf( U_pdf == inf) = 0;
                    tilde_G1  = fun_treated_PDFWatingTime( m_U, m_Dd , m_D_G1arrest, G1_pdf , U_pdf, Dd_pdf, D_G1arrest_pdf ,  obj.num_t, obj.t_span , obj.StepSize   )  ;
                    
                    %% Create tilde_S
                    lambda_Sarrest = lambda_fun(lambda_Sarrest_max);
                    lambda_d_max_S =  lambda_fun(lambda_d_max_S);
                    [ m_D_Sarrest ,D_Sarrest_pdf , ~ ]  =   PdfcdfCal_treated('Weibull',  obj.t_span ,  lambda_Sarrest  , tilde_b  );
                    D_Sarrest_pdf( D_Sarrest_pdf == inf) = 0;
                    [ m_Dd , Dd_pdf , ~ ]  =   PdfcdfCal_treated('Weibull',  obj.t_span ,  lambda_d_max_S  , tilde_b_d  );
                    Dd_pdf( Dd_pdf == inf) = 0;
                    
                   
                    tilde_S = fun_treated_PDFWatingTime( m_U, m_Dd , m_D_Sarrest, S_pdf , U_pdf, Dd_pdf, D_Sarrest_pdf ,  obj.num_t, obj.t_span , obj.StepSize   )  ;
                    parpool('local', 11);
                    parfor i = 1:num_CyclingState
                        CyclingSate =  States_convergence(i);
                        phase = mod( CyclingSate ,3 );
                        %%
                        g_est2_cdf =  obj.States_convergence_RemainingTimeDist(i,:);
                        g_est2_cdf = adjust_vector_length(g_est2_cdf, obj);
                        g_est2_cdf(g_est2_cdf>0.95) = 1;
                        
                        if  phase == 0
                            mat_M =  G2Mphasedrug_AfterTreatmentBranchingProcess_G2M_2(  mat_m , G2_arrest_UR,  G2_arrest_FR,  G2_arrest_MS,  tilde_G1,  tilde_S,  g_est2_cdf,  G1_cdf, S_cdf, G2M_cdf,  obj.num_GenOfHealthycells ,obj.StepSize,obj.t_span_AfterTreatment,   obj.num_RroundsOfTreatedcells);
                        elseif phase == 1
                            mat_M =  G2Mphasedrug_AfterTreatmentBranchingProcess_G1_2( mat_m , G2_arrest_UR,  G2_arrest_FR,  G2_arrest_MS ,   tilde_G1,   tilde_S, g_est2_cdf,   G1_cdf, S_cdf, G2M_cdf,  obj);
                        elseif phase == 2
                            mat_M =  G2Mphasedrug_AfterTreatmentBranchingProcess_S_2(  mat_m , G2_arrest_UR,  G2_arrest_FR,  G2_arrest_MS ,  tilde_G1,   tilde_S, g_est2_cdf,   G1_cdf, S_cdf, G2M_cdf,  obj.num_GenOfHealthycells ,obj.StepSize,obj.t_span_AfterTreatment,   obj.num_RroundsOfTreatedcells);
                        end
                        M_aftertreatment_col{i} = mat_M(1,:,:);
                    end
                    delete(gcp('nocreate'))
                    
            end
            if obj.TreatmentEffectOnPhases == 2 % S phase cells
                %num_GenOfHealthycells = 2 ;% one round of healty cells is mandatory
                lambda_d_max_G2M = obj.para_unknown(8);
                lambda_Sarrest_UR_max = obj.para_unknown(9);
                lambda_Sarrest_FR_max = obj.para_unknown(10);
                lambda_G2Marrest_max = obj.para_unknown(11);
                lambda_G1block_max = obj.para_unknown(12);
                tilde_b =  obj.para_unknown(13);
                tilde_b_d =  obj.para_unknown(14);
                n_S =  obj.para_unknown(15);
                EC50_d =  obj.para_unknown(16);
                EC50 =  obj.para_unknown(17);
                %Dose = 20*EC50_d;
                mat_m = Sphasedrug_BuildMatrixm_AfterTreatment( obj.Dose, obj, obj.num_GenOfHealthycells, obj.num_RroundsOfTreatedcells);
                
                %% create S arrest
                lambda_fun = @(lambda_max) lambda_max * ((obj.Dose / EC50)^n_S) / (1 + (obj.Dose / EC50)^n_S);
                
                %% create G1 block
                lambda_G1block  = lambda_fun(lambda_G1block_max );
                [~, ~, G1block ]  =   PdfcdfCal_treated('Weibull',  obj.t_span ,  lambda_G1block  , tilde_b  ); %A: scale, B , shape
                          
                 %% create S arrest
                lambda_Sarrest_UR  = lambda_fun(lambda_Sarrest_UR_max);
                [~, ~, S_arrest_UR ]  =   PdfcdfCal_treated('Weibull',  obj.t_span , lambda_Sarrest_UR  , tilde_b  ); %A: scale, B , shape
                lambda_Sarrest_FR  = lambda_fun(lambda_Sarrest_FR_max);
                [~, ~, S_arrest_FR ]  =   PdfcdfCal_treated('Weibull',  obj.t_span , lambda_Sarrest_FR  , tilde_b  ); %A: scale, B , shape
                %% Create tilde_G2M
                lambda_G2Marrest = lambda_fun(lambda_G2Marrest_max);
                lambda_d_G2M =  lambda_fun(lambda_d_max_G2M);
                [ m_D_G2Marrest ,D_G2Marrest_pdf , ~ ]  =   PdfcdfCal_treated('Weibull',  obj.t_span ,   lambda_G2Marrest   , tilde_b  );
                D_G2Marrest_pdf(D_G2Marrest_pdf  == inf) = 0;
                [ m_Dd , Dd_pdf , ~ ]  =   PdfcdfCal_treated('Weibull',  obj.t_span ,  lambda_d_G2M  , tilde_b_d  );
                Dd_pdf( Dd_pdf == inf) = 0;
                ChoiceOfDist_treated_U = 'Exponential';
                [m_U, U_pdf, ~]  =   PdfcdfCal_treated(ChoiceOfDist_treated_U, obj.t_span , 1 ,[] );
                U_pdf( U_pdf == inf) = 0;
                tilde_G2M  = fun_treated_PDFWatingTime( m_U, m_Dd , m_D_G2Marrest, G1_pdf , U_pdf, Dd_pdf, D_G2Marrest_pdf ,  obj.num_t, obj.t_span , obj.StepSize   )  ;
                for i = 1:num_CyclingState
                    CyclingSate =  States_convergence(i);
                    phase = mod( CyclingSate ,3 );
                    %%
                    g_est2_cdf =  obj.States_convergence_RemainingTimeDist(i,:);
                    g_est2_cdf = adjust_vector_length(g_est2_cdf, obj);
                    if  phase == 0
                        mat_M =  Sphasedrug_AfterTreatmentBranchingProcess_G2M_2(G1block,   mat_m , S_arrest_UR,  S_arrest_FR,  tilde_G2M, g_est2_cdf, G1_cdf, S_cdf, G2M_cdf,  obj.num_GenOfHealthycells ,obj.StepSize,obj.t_span_AfterTreatment , obj.num_RroundsOfTreatedcells);
                    elseif phase == 1
                        mat_M =  Sphasedrug_AfterTreatmentBranchingProcess_G1_2( G1block, mat_m , S_arrest_UR,  S_arrest_FR,   tilde_G2M, g_est2_cdf, G1_cdf, S_cdf, G2M_cdf,   obj.num_GenOfHealthycells ,obj.StepSize,obj.t_span_AfterTreatment , obj.num_RroundsOfTreatedcells );
                    elseif phase == 2
                        mat_M =  Sphasedrug_AfterTreatmentBranchingProcess_S_2( G1block ,  mat_m , S_arrest_UR,  S_arrest_FR,  tilde_G2M, g_est2_cdf, G1_cdf, S_cdf, G2M_cdf,   obj);
                    end
                    M_aftertreatment_col{i} = mat_M(1,:,:);
                end
                
            end
            
        end
     %%
     
     
     function [exception, SimData, CellFractions, N_d, t_span_new] =  M_matrix_plot(obj)
         %accross generation
        % full_mat_M_original =  obj.M_beforetreatement();
         [exception, SimData, CellFractions, N_d, t_span_new, mat_M_convergence_BeforeTreatment_complete,   mat_M_G1_BeforeTreatment,   M_aftertreatment_col, States_convergence   ]  = model_simulation(obj );
         % cells in different generations after the treatment 
         fun_reshape = @(x,n) reshape(x, n, length(obj.t_span_AfterTreatment) );
         if obj.TreatmentEffectOnPhases == 3
             TotalGens = obj.num_RroundsOfTreatedcells+obj.num_GenOfHealthycells; 
             NumberofCellsInDifferentGensAfterTreatment = zeros( 3*TotalGens, length(obj.t_span_AfterTreatment)  );
             
             for  k  = 1: length(States_convergence)
                 CurrentState =  States_convergence(k);
                 AbsolutePhase = mod(  CurrentState ,3);
                 M_aftertreatment_CurrentState  =    M_aftertreatment_col{k};
                 TotalState = size(M_aftertreatment_CurrentState, 2 );
                 [a_G1, a_S, a_G2M]  = G2Mphasedrug_CollectIdx_CrossGens( AbsolutePhase, obj);
                 
                 for i = 1:  TotalGens
                       NumberofCellsInDifferentGensAfterTreatment(3*(i-1) +1,:) =  NumberofCellsInDifferentGensAfterTreatment(3*(i-1) +1,:) + sum(fun_reshape(M_aftertreatment_CurrentState(1, a_G1{i},: ),  length(a_G1{i})      ),1) ;
                       NumberofCellsInDifferentGensAfterTreatment(3*(i-1) +2,:) =  NumberofCellsInDifferentGensAfterTreatment(3*(i-1) +2 ,:) + sum(fun_reshape(M_aftertreatment_CurrentState(1, a_S{i},: ),   length(a_S{i} )     ),1) ;
                       NumberofCellsInDifferentGensAfterTreatment(3*(i-1) +3,:) =  NumberofCellsInDifferentGensAfterTreatment(3*(i-1) +3 ,:) + sum(fun_reshape(M_aftertreatment_CurrentState(1, a_G2M{i},: ),  length( a_G2M{i}  )    ),1) ;
                       
                 end
             end
%                  if AbsolutePhase== 0 %The branching process is intiated by the G2 phase
%                      
%                      
%                      NumberofCellsInDifferentGensAfterTreatment(3,:) =  NumberofCellsInDifferentGensAfterTreatment(3,:) + sum(fun_reshape(M_aftertreatment_CurrentState(1, [1,3,4,5],: ), 4),1) ;
%                      NumberofCellsInDifferentGensAfterTreatment(4,:) =   NumberofCellsInDifferentGensAfterTreatment(4,:)  +  sum(fun_reshape( M_aftertreatment_CurrentState(1, [2,6],: ),2),1) ;
%                      NumberofCellsInDifferentGensAfterTreatment(5,:) =  NumberofCellsInDifferentGensAfterTreatment(5,:) +  sum(fun_reshape(M_aftertreatment_CurrentState(1, [7,8],: ),2),1) ;
%                       
%                      ctr = 1;
%                      while   ctr<= obj.num_GenOfHealthycells-1
%                          NumberofCellsInDifferentGensAfterTreatment(6+ (ctr-1)*3+1,:) =  NumberofCellsInDifferentGensAfterTreatment(6+ (ctr-1)*3+1,:) + fun_reshape(M_aftertreatment_CurrentState(1,  9 + (ctr-1)*3+1,:),1) ;
%                          NumberofCellsInDifferentGensAfterTreatment(6+ (ctr-1)*3+2,:) =  NumberofCellsInDifferentGensAfterTreatment(6+ (ctr-1)*3+2,:) + fun_reshape(M_aftertreatment_CurrentState(1,  9 + (ctr-1)*3+2,:),1) ;
%                          NumberofCellsInDifferentGensAfterTreatment(6+ (ctr-1)*3+3,:) =  NumberofCellsInDifferentGensAfterTreatment(6+ (ctr-1)*3+3,:) + fun_reshape(M_aftertreatment_CurrentState(1,  9 + (ctr-1)*3+3,:),1) ;
%                          ctr  = ctr + 1;
%                      end
%                  elseif  AbsolutePhase== 1 %The branching process is intiated by the  G1 phase
%                      NumberofCellsInDifferentGensAfterTreatment(1,:) =  NumberofCellsInDifferentGensAfterTreatment(1,:) + fun_reshape(M_aftertreatment_CurrentState(1, 1,: ),1) ;
%                      NumberofCellsInDifferentGensAfterTreatment(2,:) =  NumberofCellsInDifferentGensAfterTreatment(2,:) + fun_reshape(M_aftertreatment_CurrentState(1,2 ,: ),1) ;
%                      NumberofCellsInDifferentGensAfterTreatment(3,:) =  NumberofCellsInDifferentGensAfterTreatment(3,:) +  sum(fun_reshape(M_aftertreatment_CurrentState(1, [3,5,6,7],: ),4),1) ;
%                      NumberofCellsInDifferentGensAfterTreatment(4,:) =   NumberofCellsInDifferentGensAfterTreatment(4,:)  +  sum(fun_reshape(M_aftertreatment_CurrentState(1, [4,8],: ),2) ,1);
%                      NumberofCellsInDifferentGensAfterTreatment(5,:) =  NumberofCellsInDifferentGensAfterTreatment(5,:) +  sum(fun_reshape(M_aftertreatment_CurrentState(1, [9,10],: ),2) ,1);
%                      NumberofCellsInDifferentGensAfterTreatment(6,:) =  NumberofCellsInDifferentGensAfterTreatment(6,:)+  fun_reshape(M_aftertreatment_CurrentState(1, 11,: ),1) ;
%                      ctr = 1;
%                      while   ctr<= obj.num_GenOfHealthycells-1
%                          NumberofCellsInDifferentGensAfterTreatment(6+ (ctr-1)*3+1,:) =  NumberofCellsInDifferentGensAfterTreatment(6+ (ctr-1)*3+1,:) + fun_reshape(M_aftertreatment_CurrentState(1,   11 + (ctr-1)*3+1,:),1) ;
%                          NumberofCellsInDifferentGensAfterTreatment(6+ (ctr-1)*3+2,:) =  NumberofCellsInDifferentGensAfterTreatment(6+ (ctr-1)*3+2,:) + fun_reshape(M_aftertreatment_CurrentState(1,   11 + (ctr-1)*3+2,:),1) ;
%                          NumberofCellsInDifferentGensAfterTreatment(6+ (ctr-1)*3+3,:) =  NumberofCellsInDifferentGensAfterTreatment(6+ (ctr-1)*3+3,:) + fun_reshape(M_aftertreatment_CurrentState(1,  11 + (ctr-1)*3+3,:),1) ;
%                          ctr  = ctr + 1;
%                      end
%                   elseif  AbsolutePhase==2 %The branching process is intiated by the  S phase
%                      NumberofCellsInDifferentGensAfterTreatment(2,:) =  NumberofCellsInDifferentGensAfterTreatment(2,:) + fun_reshape(M_aftertreatment_CurrentState(1,1 ,: ),1) ;
%                      NumberofCellsInDifferentGensAfterTreatment(3,:) =  NumberofCellsInDifferentGensAfterTreatment(3,:) +  sum(fun_reshape(M_aftertreatment_CurrentState(1, [2,4,5,6],: ),4) ,1);
%                      NumberofCellsInDifferentGensAfterTreatment(4,:) =   NumberofCellsInDifferentGensAfterTreatment(4,:)  +  sum( fun_reshape(M_aftertreatment_CurrentState(1, [3,7],: ),2) ,1);
%                      NumberofCellsInDifferentGensAfterTreatment(5,:) =  NumberofCellsInDifferentGensAfterTreatment(5,:) +  sum( fun_reshape(M_aftertreatment_CurrentState(1, [8,9], : ) ,2),1) ;
%                      NumberofCellsInDifferentGensAfterTreatment(6,:) =  NumberofCellsInDifferentGensAfterTreatment(6,:)+ fun_reshape(M_aftertreatment_CurrentState(1, 10,: ),1) ;
%                      ctr = 1;
%                      while   ctr<= obj.num_GenOfHealthycells-1
%                          NumberofCellsInDifferentGensAfterTreatment(6+ (ctr-1)*3+1,:) =  NumberofCellsInDifferentGensAfterTreatment(6+ (ctr-1)*3+1,:) + fun_reshape(M_aftertreatment_CurrentState(1,   10+ (ctr-1)*3+1,:),1) ;
%                          NumberofCellsInDifferentGensAfterTreatment(6+ (ctr-1)*3+2,:) =  NumberofCellsInDifferentGensAfterTreatment(6+ (ctr-1)*3+2,:) + fun_reshape(M_aftertreatment_CurrentState(1,   10 + (ctr-1)*3+2,:),1) ;
%                          NumberofCellsInDifferentGensAfterTreatment(6+ (ctr-1)*3+3,:) =  NumberofCellsInDifferentGensAfterTreatment(6+ (ctr-1)*3+3,:) + fun_reshape(M_aftertreatment_CurrentState(1,  10 + (ctr-1)*3+3,:),1) ;
%                          ctr  = ctr + 1;
%                      end
%                  end
             f1 = figure;
             
             fontName = 'Arial';
             fontSize = 14;
             
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
             legend(legendEntries, 'FontSize', fontSize-1, 'Location', 'bestoutside'); % Adjusted font size and location
             xlabel('Time(hr)', 'FontName', fontName, 'FontSize', fontSize);
             ylabel('Number of Cells', 'FontName', fontName, 'FontSize', fontSize); % Modified label
             %title('Your Title Here', 'FontName', fontName, 'FontSize', fontSize+1); % Add your title
             
             grid on; % Add grid for better readability
             
             % Adjust axis properties
             set(gca, 'FontName', fontName, 'FontSize', fontSize, 'LineWidth', 1.5);
             exportgraphics( f1 ,  fullfile(obj.outputfile_dir,'Gen.png' ) ,'Resolution' , 300   );
         
             
         end
         
         if obj.TreatmentEffectOnPhases == 2
             TotalGens =  obj.num_RroundsOfTreatedcells+obj.num_GenOfHealthycells  ;
             NumberofCellsInDifferentGensAfterTreatment = zeros( 3* TotalGens, length(obj.t_span_AfterTreatment)    ); % The extra one state is for the additional G2M phase
             for  k  = 1: length(States_convergence)
                 CurrentState =  States_convergence(k);
                 AbsolutePhase = mod(  CurrentState ,3);
                 M_aftertreatment_CurrentState  =    M_aftertreatment_col{k};
                 TotalState = size(M_aftertreatment_CurrentState, 2 );
                 [a_G1, a_S, a_G2M]  = Sphasedrug_CollectIdx_CrossGens( AbsolutePhase, obj);
                 
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
             %                  if AbsolutePhase== 0 %The branching process is intiated by the G2 phase
             %                      NumberofCellsInDifferentGensAfterTreatment(1,:) =  NumberofCellsInDifferentGensAfterTreatment(1,:) + sum(fun_reshape(M_aftertreatment_CurrentState(1, 1,: ), 1),1) ; % first G2M
             %
             %                      NumberofCellsInDifferentGensAfterTreatment(2,:) =  NumberofCellsInDifferentGensAfterTreatment(2,:) + sum(fun_reshape(M_aftertreatment_CurrentState(1, [2,3],: ), 2),1) ; % first G1
             %                      NumberofCellsInDifferentGensAfterTreatment(3,:) =   NumberofCellsInDifferentGensAfterTreatment(3,:)  +  sum(fun_reshape( M_aftertreatment_CurrentState(1, [ 4, 6, 7 ],: ), 3 ),1) ; %S
             %                      NumberofCellsInDifferentGensAfterTreatment(4,:) =  NumberofCellsInDifferentGensAfterTreatment(4,:) +  sum(fun_reshape(M_aftertreatment_CurrentState(1, [5,8],: ),2),1) ; %G2M
             %
             %                      NumberofCellsInDifferentGensAfterTreatment(5,:) =  NumberofCellsInDifferentGensAfterTreatment(5,:)+  sum( fun_reshape( M_aftertreatment_CurrentState(1, [9,10],: ),2) ,1);% second G1'
             %                      NumberofCellsInDifferentGensAfterTreatment(6,:) =  NumberofCellsInDifferentGensAfterTreatment(6,:) + sum(fun_reshape(M_aftertreatment_CurrentState(1, [11, 13,14],: ), 3),1) ; % second S
             %                      NumberofCellsInDifferentGensAfterTreatment(7,:) =   NumberofCellsInDifferentGensAfterTreatment(7,:)  +  sum(fun_reshape( M_aftertreatment_CurrentState(1, [ 12, 15 ],: ), 2 ),1) ; %G2M
             %                      ctr = 1;
             %                      while   ctr<= obj.num_GenOfHealthycells-1
             %                          NumberofCellsInDifferentGensAfterTreatment(7+ (ctr-1)*3+1,:) =  NumberofCellsInDifferentGensAfterTreatment(7+ (ctr-1)*3+1,:) + fun_reshape(M_aftertreatment_CurrentState(1,  15 + (ctr-1)*3+1,:),1) ;
             %                          NumberofCellsInDifferentGensAfterTreatment(7+ (ctr-1)*3+2,:) =  NumberofCellsInDifferentGensAfterTreatment(7+ (ctr-1)*3+2,:) + fun_reshape(M_aftertreatment_CurrentState(1,  15 + (ctr-1)*3+2,:),1) ;
             %                          NumberofCellsInDifferentGensAfterTreatment(7+ (ctr-1)*3+3,:) =  NumberofCellsInDifferentGensAfterTreatment(7+ (ctr-1)*3+3,:) + fun_reshape(M_aftertreatment_CurrentState(1,  15 + (ctr-1)*3+3,:),1) ;
             %                          ctr  = ctr + 1;
             %                      end
             %                  elseif  AbsolutePhase== 1 %The branching process is intiated by the  G1 phase
             %                      NumberofCellsInDifferentGensAfterTreatment(2,:) =  NumberofCellsInDifferentGensAfterTreatment(2,:) + sum( fun_reshape(M_aftertreatment_CurrentState(1, [1,2 ],: ),2),1);
             %                      NumberofCellsInDifferentGensAfterTreatment(3,:) =  NumberofCellsInDifferentGensAfterTreatment(3,:) + sum(fun_reshape(M_aftertreatment_CurrentState(1,[3,5,6] ,: ),3) ,1) ;
             %                      NumberofCellsInDifferentGensAfterTreatment(4,:) =  NumberofCellsInDifferentGensAfterTreatment(4,:) +  sum(fun_reshape(M_aftertreatment_CurrentState(1, [4,7],: ),2),1) ;
             %                      NumberofCellsInDifferentGensAfterTreatment(5,:) =   NumberofCellsInDifferentGensAfterTreatment(5,:)  +  sum(fun_reshape(M_aftertreatment_CurrentState(1, [8,9],: ),2) ,1);
             %                      NumberofCellsInDifferentGensAfterTreatment(6,:) =  NumberofCellsInDifferentGensAfterTreatment(6,:) +  sum(fun_reshape(M_aftertreatment_CurrentState(1, [10,12,13],: ),3) ,1);
             %                      NumberofCellsInDifferentGensAfterTreatment(7,:) =  NumberofCellsInDifferentGensAfterTreatment(7,:)+  sum(fun_reshape(M_aftertreatment_CurrentState(1, [11,14],: ),2),1) ;
             %                      ctr = 1;
             %                      while   ctr<= obj.num_GenOfHealthycells-1
             %                          NumberofCellsInDifferentGensAfterTreatment(7+ (ctr-1)*3+1,:) =  NumberofCellsInDifferentGensAfterTreatment(7+ (ctr-1)*3+1,:) + fun_reshape(M_aftertreatment_CurrentState(1,   14 + (ctr-1)*3+1,:),1) ;
             %                          NumberofCellsInDifferentGensAfterTreatment(7+ (ctr-1)*3+2,:) =  NumberofCellsInDifferentGensAfterTreatment(7+ (ctr-1)*3+2,:) + fun_reshape(M_aftertreatment_CurrentState(1,   14 + (ctr-1)*3+2,:),1) ;
             %                          NumberofCellsInDifferentGensAfterTreatment(7+ (ctr-1)*3+3,:) =  NumberofCellsInDifferentGensAfterTreatment(7+ (ctr-1)*3+3,:) + fun_reshape(M_aftertreatment_CurrentState(1,  14 + (ctr-1)*3+3,:),1) ;
             %                          ctr  = ctr + 1;
             %                      end
             %                   elseif  AbsolutePhase==2 %The branching process is intiated by the S phase
             %                      NumberofCellsInDifferentGensAfterTreatment(3,:) =  NumberofCellsInDifferentGensAfterTreatment(3,:) + sum(fun_reshape(M_aftertreatment_CurrentState(1,[1,3,4] ,: ),3),1) ;%S
             %                      NumberofCellsInDifferentGensAfterTreatment(4,:) =  NumberofCellsInDifferentGensAfterTreatment(4,:) +  sum(fun_reshape(M_aftertreatment_CurrentState(1, [2,5],: ),2) ,1);%G2M
             %                      NumberofCellsInDifferentGensAfterTreatment(5,:) =   NumberofCellsInDifferentGensAfterTreatment(5,:)  +  sum( fun_reshape(M_aftertreatment_CurrentState(1, [6,7],: ),2) ,1); % second G1
             %                      NumberofCellsInDifferentGensAfterTreatment(6,:) =  NumberofCellsInDifferentGensAfterTreatment(6,:) +  sum( fun_reshape(M_aftertreatment_CurrentState(1, [8,10,11], : ) ,3),1) ; %second S
             %                      NumberofCellsInDifferentGensAfterTreatment(7,:) =  NumberofCellsInDifferentGensAfterTreatment(7,:)+ sum( fun_reshape(M_aftertreatment_CurrentState(1, [9, 12],: ),2),1 );%second G2M
             %                      ctr = 1;
             %                      while   ctr<= obj.num_GenOfHealthycells-1
             %                          NumberofCellsInDifferentGensAfterTreatment(7+ (ctr-1)*3+1,:) =  NumberofCellsInDifferentGensAfterTreatment(7+ (ctr-1)*3+1,:) + fun_reshape(M_aftertreatment_CurrentState(1,   12+ (ctr-1)*3+1,:),1) ;
             %                          NumberofCellsInDifferentGensAfterTreatment(7+ (ctr-1)*3+2,:) =  NumberofCellsInDifferentGensAfterTreatment(7+ (ctr-1)*3+2,:) + fun_reshape(M_aftertreatment_CurrentState(1,   12 + (ctr-1)*3+2,:),1) ;
             %                          NumberofCellsInDifferentGensAfterTreatment(7+ (ctr-1)*3+3,:) =  NumberofCellsInDifferentGensAfterTreatment(7+ (ctr-1)*3+3,:) + fun_reshape(M_aftertreatment_CurrentState(1,  12 + (ctr-1)*3+3,:),1) ;
             %                          ctr  = ctr + 1;
             %                      end
             %                  end
             %               end
             f2 = figure('Position', [100, 100, 1200, 800]); % Increase figure size
             
             numGenerations = TotalGens;
             colors = lines(numGenerations);
             lineStyles = {'-', ':', '--'};
             legendEntries = {};
             
             hold on; % Move the hold command up to ensure all subsequent plots are held
             
             for gen = 1: numGenerations
                 for cellType = 1:3  % 1 for G1, 2 for S, 3 for G2M
                     if gen == 1 && (cellType == 1 || cellType == 2)
                         continue
                     end
                     
                     xData = t_span_new;
                     yData = NumberofCellsInDifferentGensAfterTreatment((gen-1)*3 + cellType, :);
                     
                     plot(xData, yData, 'LineStyle', lineStyles{cellType}, 'Color', colors(gen,:), 'LineWidth', 2); % Increase linewidth
                     
                     cellTypeStr = {'G1', 'S', 'G2M'};
                     legendEntries{end + 1} = [cellTypeStr{cellType} ' - Gen ' num2str(gen)];
                 end
             end
             
             xlabel('Time (hr)', 'FontSize', 16);
             ylabel('Number of Cells', 'FontSize', 16); % Adjust as needed
             title('Distribution of Cells Across Generations', 'FontSize', 20); % Add a title
             
             set(gca, 'FontSize', 14, 'LineWidth', 1.5); % Adjust fonts and linewidth
             grid on; % Add a grid
             
             legend(legendEntries, 'Location', 'northeastoutside', 'FontSize', 12); % Adjust legend location
             %tight_layout(); % Adjust layout to prevent overlap
             
         end
         %          p = 1;
         %          f1 = figure('Position', [1440,748,1122,589]);
         %          p1 = plot(obj.t_span,  reshape(full_mat_M_original(p ,States_convergence(mod(States_convergence ,3) ==0 ),: ), length( States_convergence(mod(States_convergence ,3) ==0 ) ), obj.num_t)' , 'LineWidth', 1.5 );
         %          hold on
         %          p2 = plot(obj.t_span,  M_afterTreatment_complete( mod(States_convergence ,3) ==0 ,:)' , '--', 'LineWidth', 1.5 );
         %          for i = 1:length(p1)
         %              p2(i).Color = p1(i).Color;
         %          end
         %          xline( 160 ,'--', ['treatment time:'  num2str(obj.clinical_time),'h']  ,'LineWidth',3,'LabelOrientation','horizontal' );
         %          legend( append('M_{1,', string(States_convergence(mod(States_convergence ,3) ==0)  ),'}'  ),'FontSize', 15 );
         %          exportgraphics( f1 ,  fullfile(obj.outputfile_dir,'CyclingG2M.png' ) ,'Resolution' , 300   );
         %
         %
         %          f4 = figure('Position', [1440,748,1122,589]);
         %          p1 = plot(obj.t_span,  reshape(full_mat_M_original(p , States_convergence(mod(States_convergence ,3) ~=0 ),: ), length(States_convergence(mod(States_convergence ,3) ~=0 ) ), obj.num_t)' , 'LineWidth', 1.5 );
         %          hold on
         %          p2 = plot(obj.t_span,  M_afterTreatment_complete( mod(States_convergence ,3) ~=0 ,:)' , '--', 'LineWidth', 1.5 );
         %          for i = 1:length(p1)
         %              p2(i).Color = p1(i).Color;
         %          end
         %          xline( 160 ,'--', ['treatment time:'  num2str(obj.clinical_time),'h']  ,'LineWidth',3,'LabelOrientation','horizontal' );
         %          legend( append('M_{1,', string(States_convergence(mod(States_convergence ,3) ~=0)  ),'}'  ),'FontSize', 15 );
         %          exportgraphics( f4 ,  fullfile(obj.outputfile_dir,'CyclingNonG2M.png' ) ,'Resolution' , 300   );
         %
         %
         %          f3 = figure('Position', [1440,748,1122,589]);
         %          p1 = plot(obj.t_span,  reshape(full_mat_M_original(p , States_convergence(end) +[1,2,3,4,5],: ), 5 , obj.num_t)' , 'LineWidth', 1.5 );
         %          hold on
         %          p2 = plot(obj.t_span, reshape(next_mat_M_G1(p , States_convergence(end) +[1,2,3,4,5],: ), 5 ,obj.num_t)'    , '--', 'LineWidth', 1.5 );
         %          for i = 1:length(p1)
         %              p2(i).Color = p1(i).Color;
         %          end
         %          xline( 160 ,'--', ['treatment time:'  num2str(obj.clinical_time),'h']  ,'LineWidth',3,'LabelOrientation','horizontal' );
         %          legend( append('M_{1,', string(  States_convergence(end) +[1,2,3,4,5] ),'}'  ),'FontSize', 15 );
         %          exportgraphics( f3 ,  fullfile(obj.outputfile_dir,'Newborn.png' ) ,'Resolution' , 300   );
         
         
         %%
         colorMatrix = [[0 0.4470 0.7410]; [0.8500 0.3250 0.0980]; [0.4940 0.1840 0.5560]];
         f4 = figure('Position', [100, 100, 1200, 800]); % Larger figure size
         hold on;
         
         % Plot simulation lines
         for i = 1:3
             plot(t_span_new, CellFractions(i,:), 'Color', colorMatrix(i,:) ,'LineWidth', 2.5);
         end
         
         % Plot experimental data points
         markerSymbols = {'o', 's', 'd'};
         for i = 1:3
             s = scatter(obj.data_time, obj.current_data(:,i), 100, 'filled', markerSymbols{i}); % Increase marker size
             s.MarkerEdgeColor = colorMatrix(i,:);
             s.MarkerFaceColor = colorMatrix(i,:);
         end
         
         % Set Labels, Legend, Title
         xlabel('Hours', 'FontSize', 18);
         ylabel('Cell Fraction', 'FontSize', 18); % Assuming it's a fraction; adjust as needed
         title('Model Simulation vs. Experimental Data', 'FontSize', 20);
         legend({'G1', 'S', 'G2M', 'data G1', 'data S', 'data G2M'}, 'Location', 'northwest', 'FontSize', 16);
         
         % Styling
         set(gca, 'FontSize', 16, 'LineWidth', 1.5); % Set font size and axes line width
         grid on; % Turn on the grid
         
         % Save the figure
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
         
         
         a = (obj.current_data  - SimData) ;
         vecnorm(a(:,1),1)
         vecnorm(a(:,2),1)
         vecnorm(a(:,3),1)
         
     end
     
     %%
     function   [exception, SimData, CellFractions, N_d, t_span_new, mat_M_convergence_BeforeTreatment_complete,   mat_M_G1_BeforeTreatment,   M_aftertreatment_col, States_convergence   ]  = model_simulation(obj )
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
         if obj.TreatmentEffectOnPhases == 3 %G2M phase drug
             for k = 1: length(States_convergence)
                 num_cell_AtTreatment = mat_M_convergence_BeforeTreatment_complete(k,  idx);
                 CurrentState =  States_convergence(k);
                 AbsolutePhase = mod(  CurrentState ,3);
                 M_aftertreatment_CurrentState  =    M_aftertreatment_col{k};
                 TotalState = size(M_aftertreatment_CurrentState, 2 );
                 
                 G2M_idx_col = G2Mphasedrug_CollectIdx( AbsolutePhase, obj, 'G2M');
                 G1_idx_col = G2Mphasedrug_CollectIdx( AbsolutePhase, obj, 'G1');
                 S_idx_col = G2Mphasedrug_CollectIdx( AbsolutePhase, obj, 'S');
                 
                 M_aftertreatment_S = fun_reshape(  M_aftertreatment_CurrentState(1,  S_idx_col , : ), length( S_idx_col  )  );
                 M_aftertreatment_G2M = fun_reshape(  M_aftertreatment_CurrentState(1,  G2M_idx_col , : ), length( G2M_idx_col  )  );
                 M_aftertreatment_G1 = fun_reshape(  M_aftertreatment_CurrentState(1,  G1_idx_col , : ), length( G1_idx_col  )  );
                 N_col_S_AcrossGen = N_col_S_AcrossGen +   num_cell_AtTreatment* sum(M_aftertreatment_S,1 ) ;
                 N_col_G2M_AcrossGen = N_col_G2M_AcrossGen +   num_cell_AtTreatment* sum(M_aftertreatment_G2M,1) ;
                 N_col_G1_AcrossGen = N_col_G1_AcrossGen +   num_cell_AtTreatment * sum( M_aftertreatment_G1,1) ;
                 N_d =  N_d + reshape(M_aftertreatment_CurrentState(1,end,:), 1, num_t_new);  %+  N_0(2)*reshape(mat_M_S(2,end,:), 1,num_t_new) + N_0(3)*reshape(mat_M_G2M(3,end,:), 1,num_t_new);
                
             end
         end
         
         if obj.TreatmentEffectOnPhases == 2 %S phase drug
             for k = 1: length(States_convergence)
                 num_cell_AtTreatment = mat_M_convergence_BeforeTreatment_complete(k,  idx);
                 CurrentState =  States_convergence(k);
                 AbsolutePhase = mod(  CurrentState ,3);
                 M_aftertreatment_CurrentState  =    M_aftertreatment_col{k};
                 TotalState = size(M_aftertreatment_CurrentState, 2 );
                 G2M_idx_col = Sphasedrug_CollectIdx( AbsolutePhase, obj, 'G2M');
                 G1_idx_col = Sphasedrug_CollectIdx( AbsolutePhase, obj, 'G1');
                 S_idx_col = Sphasedrug_CollectIdx( AbsolutePhase, obj, 'S');
                 
                 
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
         
         
          try
             SimData  = interp1( t_span_new , CellFractions',  obj.data_time ,'pchip') ;
             if ~isempty(SimData)
                 exception = [];
             end
         catch exception
             SimData  = NaN;
         end
         
     end


  end
    
    
end
