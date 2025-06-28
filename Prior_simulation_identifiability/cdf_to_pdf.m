clear
close all

addpath('../Sobol_3_G2M') 
addpath('../Sobol_4_G2M/ReadYAML-master/ReadYAML-master/')

ymlflie  = 'prior_identifiability_config_41.yaml';
datastr = ReadYaml(ymlflie);
Model_Index        = datastr.Model_Index;
TreatmentTarget    = datastr.TreatmentTarget;
batch_number       = datastr.batch_number;
num_batches        = datastr.num_batches;
num_worker         = datastr.num_worker;
load_history_flag  = datastr.load_history_flag;
%history_folder     = datastr.history_folder;
T =  datastr.T; 
Total_sample_limit = datastr.Total_sample_limit;


fprintf('batch_number %d',  batch_number)

history_folder  = get_latest_batch_file(batch_number, Model_Index)
fprintf('history_folder %s',  history_folder)



rootdir  = fullfile( '../Predictive_check' , sprintf( 'samples_%s_EC50_fixed', TreatmentTarget ) ) ;
%rootdir = 'GSA_outputs2';
allFiles = dir(rootdir);
dirFlags = [allFiles.isdir];
subFolders = allFiles(dirFlags);
subFolders  = subFolders(3:end);
M_A_col = [];
L_A_col = [];
for i = 1:length(subFolders  )
    load(fullfile(rootdir, subFolders(i).name , 'predictive_check.mat'),  'M_A', 'L_A');%% calculate the overall response variance
    M_A_col = [M_A_col ;  M_A];
    L_A_col = [L_A_col; L_A];
end

M_A_col = M_A_col(1: Total_sample_limit);
L_A_col  =  L_A_col(1: Total_sample_limit, :);
Total_number_points = size( L_A_col,1);
disp(Total_number_points)
batchSize = ceil(Total_number_points / num_batches);

para_label = splitYamlString(datastr.para_label);
para_exclude = splitYamlString(datastr.para_exclude);


num_samples = length( M_A_col );
[num_species, num_timepoints] = size(M_A_col{1});

if Model_Index == 1
    lb =  [
        1e-6*ones(9,1)
        0.1*ones(9,1)
        0.05*ones(3,1)
        0.1*ones(2,1) 
        7e-4*ones(3,1)
        ];
    
    ub = [
        ones(9,1)
        100*ones(7,1)
        25*ones(2,1)
        5
        5
        5
        10*ones(2,1)  
        9e-4*ones(3,1)
        ];

elseif Model_Index == 41
    lb =  [
          1e-6*ones(9,1)
          0.1*ones(9,1)
          0.1
          50
          50
          0.1*ones(2,1) 
          5e-4*ones(3,1)
        ];
        
    ub = [
            ones(9,1)
            100*ones(7,1)
            25*ones(2,1)
            5 %n_i
            700
            700
            10*ones(2,1)  
            9e-4*ones(3,1)  
     ];
end


all_outputfile_dir = [sprintf('Batch%d', batch_number), '/prior_simulation_identifiability_batch', num2str(batch_number), 'Model', num2str(Model_Index), '_', datestr(now,'yyyy_mm_dd_HH_MM')]; 

fprintf('new folder is %s', all_outputfile_dir)

startIdx = (batch_number - 1) * batchSize + 1;
endIdx = min(batch_number * batchSize, Total_number_points);

sample_numbers = startIdx : endIdx

% Create cell array of sample names
expected_samples = arrayfun(@(x) sprintf('sample%d', x), sample_numbers, 'UniformOutput', false);


para_change_exclude_idx = find(contains( para_label, para_exclude));
para_fixed = [1.1475 , 1.1489 , 1.1489 ];

hist_MCMC_samples_folder_coll = {};
if load_history_flag == 1
    [prioritized_samples, hist_MCMC_samples_folder_coll]   = get_index(history_folder, expected_samples );
    disp(prioritized_samples)
    prioritized_samples_idx = sample_numbers(prioritized_samples)
else
    prioritized_samples_idx = sample_numbers;
end



%for i = startIdx: endIdx
% Get local cluster object
localCluster = parcluster('local');

localCluster.NumWorkers = num_worker;

% Save the updated profile
saveProfile(localCluster);


parpool('local', num_worker)
parfor k = 1:length(prioritized_samples_idx) 
    j = prioritized_samples_idx(k);
    find_hist_flag = 0
    if load_history_flag == 1
        % Create the search pattern
        search_str = ['/sample' num2str(j) '/'];
    
        % Find indices of matches
        hist_set_idx = find(contains(hist_MCMC_samples_folder_coll, search_str));
    
        if ~isempty(hist_set_idx)
            find_hist_flag = 1;
            file_path = hist_MCMC_samples_folder_coll{hist_set_idx(1)};
            s = load(file_path, 'Para_col_accepted', 'S_curr' ,'Para_col_all' , 'Diff_norm_col' ); 
            Para_col_accepted = s.Para_col_accepted;
            S_curr = s.S_curr;
            Para_col_all = s.Para_col_all;
            Diff_norm_col = s.Diff_norm_col;
        end
    end
    
    task = getCurrentTask;
    taskid = task.ID;
    true_para  =  L_A_col(j,:); 
    synthetic_data = M_A_col{j};
    
   
    outputfile_dir= fullfile(all_outputfile_dir, "log_"+ num2str( taskid) );

    if ~isfolder(outputfile_dir)
        mkdir(outputfile_dir)
    end
    
    
    outputfile_dir_indivi_sample =  fullfile(all_outputfile_dir, "log_"+ num2str( taskid),  "sample"+num2str(j) ); 
    if ~isfolder(outputfile_dir_indivi_sample )
        mkdir(outputfile_dir_indivi_sample )
    end
    
    outfile = fopen(  fullfile(outputfile_dir,  'outlog.log' ),'a');
    fprintf(outfile, 'Start sample %d \n', j );


    if find_hist_flag == 1
        fprintf(outfile, 'load: %s',  file_path )
        Para_col_accepted_hist =  Para_col_accepted;
        nonzeroRows = all(Para_col_accepted_hist,2);
        Para_col_accepted_hist = Para_col_accepted_hist(nonzeroRows, :);  
        Diff_norm_col = Diff_norm_col(nonzeroRows(2:end));  
        
        Para_col_all_hist = Para_col_all(nonzeroRows, :);  
         
        num_HistPoints = size(Para_col_accepted_hist ,1);
        fprintf(outfile, 'num of hist points: %d', num_HistPoints)
        num_para = size(Para_col_accepted_hist  ,2);
        Para_col =   [ Para_col_accepted_hist; zeros( T, num_para ) ];
        Diff_norm_col_new  =   [ Diff_norm_col; zeros( T,  1) ];
        Para_col_all_new  =   [ Para_col_all_hist; zeros( T,    num_para )  ];
        
        fprintf(outfile, 'size of para_col: %d, %d', size(Para_col,1 ),  size(Para_col,2 ) );
        
        curr_T = num_HistPoints;
    else
        fprintf(outfile, 'start from 0')
        num_HistPoints = 0;

        init_para =  unifrnd(lb,ub);
        sum_first_three = sum( init_para(1:3));
        init_para(1:3) =init_para(1:3) / sum_first_three;
    
        sum_first_three = sum( init_para([4,8]));
        init_para([4,8]) =init_para([4,8]) / sum_first_three;
    
        init_para = init_para';
        num_para = length(init_para);

        Para_col = zeros( T+1,    num_para - length(para_fixed) ) ;
        init_para_exclude = init_para ;
        init_para_exclude(para_change_exclude_idx) = [];
        Para_col(1,:) = init_para_exclude ; %Para_col_hist(end,:);
        
        Para_col_all_new = zeros(size(Para_col ));
        Diff_norm_col_new = zeros(T,1);
        curr_T =  1;
    end

 
    if find_hist_flag  == 0
        fprintf(outfile, 'create S_curr' );
        %a lower diagonal matrix with positive diagonal elements
        S_curr = zeros(num_para);
        %subset1 = find(contains( para_label    ,["$"d_{G1}$"" "$"d_S$""  "$"d_{G2M}$""   "$"p1_{G1}$""  "$"p1_{S}$""  "$"p3_{G2M}$""  "$"p2_{G2M}$""  "$"p4_{G2M}$"" ] ));
        if Model_Index == 1
            subset2 = find(contains( para_label    ,["$q_1$"
                "$q_2$"
                "$q_3$"
                "$m_{d,UR}^{max}$"
                "$m_{d,FR}^{max}$"
                "$m_{d,MS}^{max}$"
                "$m_{d,UR,G1}^{max}$"
                "$m_{d,UR,S}^{max}$"
                "$n_{G2M}$"
                "$EC_{50,d}$"
                "$EC_{50}$"] ));
        elseif Model_Index ==41
            subset2 = find(contains( para_label    ,["$q_1$"
                "$q_2$"
                "$q_3$"
                "$q_4$"
                "$m_{d,UR}^{max}$"
                "$m_{d,FR}^{max}$"
                "$m_{d,MS}^{max}$"
                "$m_{d,UR,G1}^{max}$"
                "$m_{d,UR,S}^{max}$"
                "$n_{G2M}$"
                ] ));
        end
    
        subset4 = find(contains( para_label , [ "$\theta_1$"   "$\theta_2$"   "$\theta_3$" ] ));
        % subset3 = 3:8;
        for i = 1: num_para
            S_curr(i,1:i) =  randn*5;
            %         if ismember(i, subset1  )
            %             S_curr(i,1:i) =  abs(randn)*0.01;
            %         end
    
            if ismember(i, subset2  )
                temp = rand;
                if temp > 0.5
                    temp =    temp*0.5;
                end
                S_curr(i,1:i) = temp;
            end
            %          if ismember(i,  subset3  )
            %             S_curr(i,1:i) =  abs(randn)*0.5;
            %         end
    
            if ismember(i,  subset4  )
                S_curr(i,1:i) =  abs(randn)*0.0005;
            end
        end
        S_curr  = S_curr( 1: num_para, 1:num_para );
        S_curr(para_change_exclude_idx ,  :  ) = [];
        S_curr(:,  para_change_exclude_idx    ) = []; 
    end
    lb_local = lb; 
    ub_local = ub; 
    lb_local(para_change_exclude_idx ) = [];
    ub_local(para_change_exclude_idx ) = [];
    myfun = @Target_fun_Treated_prior_identifiability;

    [Para_col_accepted, acc_rate_col, S_curr, Para_col_all, Diff_norm_col] = ...
        MyRAM_treated_prior_identifiability(...
            T, size(Para_col,2), Para_col, S_curr, myfun, para_label, outputfile_dir_indivi_sample , Model_Index, ...
            num_HistPoints,  lb_local, ub_local , synthetic_data, para_change_exclude_idx, para_fixed,  outfile, j , curr_T, Diff_norm_col_new, Para_col_all_new );

    %if load_history_flag == 1
     %   Para_col_accepted= [Para_col_accepted_hist; Para_col_accepted(2:end,:)   ];
     %   Para_col_all = [Para_col_all_hist ; Para_col_all(2:end,:)   ];
   % end
   
    %save( fullfile( outputfile_dir_indivi_sample,  strcat('Para_treated_', num2str(Model_Index), 'sample', num2str(j),   '.mat' ) )  ,'Para_col_accepted'  , 'acc_rate_col',  'S_curr' , 'Para_col_all','Diff_norm_col' ,'-v7.3');
       % Store results in a struct
    result_struct = struct( ...
        "Para_col_accepted", Para_col_accepted, ...
        "acc_rate_col", acc_rate_col, ...
        "S_curr", S_curr, ...
        "Para_col_all", Para_col_all, ...
        "Diff_norm_col", Diff_norm_col ...
    );

    % Construct filename for saving
    filename = fullfile(outputfile_dir_indivi_sample, ...
        sprintf("Para_treated_%d_sample%d.mat", Model_Index, j));

    % Save using -struct (equivalent to -fromstruct in the code pattern you referenced)
    save(filename, "-fromstruct", result_struct, '-v7.3');

    fprintf(outfile, 'Saved sample %d \n', j  ); 
    fprintf(outfile, '====================' ); 
end
