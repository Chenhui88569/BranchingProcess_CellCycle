% Global sensitivity analysis
% Model output is the kinectics of the tumor proliferating cells + damaging cells(100 mg/kg weekly )
% step1. choose the most influential parameter
% step2. Variance-based sensitivity analysis
% purturb the parameter with the range (determine the range)
% create 1000 parameter sets using Latin hypercube sample
close all
clear 
%dbstop if error

addpath('../Sobol_2_S')
num_worker = 15;
num_sample = 200; %Total number of samples is 500
Model_Index = 42;



outputfile_dir = ['samples_S_EC50_fixed/Predictive_Check_CellCycle_Model_with_parameters', num2str(Model_Index) ,'_' datestr(now,'yyyy_mm_dd_HH_MM')];
if ~isfolder(outputfile_dir)
    mkdir(outputfile_dir)
else 
    rmdir(outputfile_dir,'s');
     mkdir(outputfile_dir)
end


addpath('../Sobol_2_S/ReadYAML-master/ReadYAML-master/')
ymlfile = ['DatasetsInfo'  num2str( Model_Index )  '.yml'];
datastr = ReadYaml(ymlfile);

%% load all the parameters 
%MCMC_obj = MCMC_sample_analysis(Model_Index);
% [Para_trun, Para_mean_col ] = MCMC_obj.sample_trun( );
% para_all_label = MCMC_obj.para_label;
% para_all_label = para_all_label(1:end-3);

para_all_label =  ["q_{1, S}"
        "q_{2, S}"
        "q_{3, S}"
        "m_{d,FR, S}^{max}"
        "m_{d,UR, S}^{max}"
        "m_{d,G2/Marrest}^{max}"
        "m_{d,G1block}^{max}"
        "\lambda_{d,G2/M}^{max}"
        "\lambda_{UR,S}^{max}"
        "\lambda_{FR,S}^{max}"
        "\lambda_{G2/Marrest}^{max}"
        "\lambda_{G1block}^{max}"
        "b_{S}"
        "b_{d,S}"
        "n_{S}"
        "EC_{50,d, S}"
        "EC_{50, S}"
        "a_{G2/M}"
        "b_{G2/M}"
        "theta_1"
        "theta_2"
        "theta_3"
];


Para_mean_col = [0.794640761865667	0.234743947792112	0.711938836096859	0.140777190181406	0.521746170025760	0.520124560016750	0.343708610726482	25.0406526783250	14.9822738069382	46.1891808219113	11.1718270222620	83.1671628832383	12.6855972819988	19.4549767021315	2.77992552387977	169.966388987478	64.4414618246546	4.19702043214536	3.88132867758518	0.000672619165833518	0.000623650760682560	0.000659379806498318];
Para_Value =  Para_mean_col(1:end-3);

para_change =  para_all_label(1:end-3);

para_exclude =  [ "n_{S}"
        "EC_{50,d, S}"
        "EC_{50, S}" ];

para_change_exclude_idx = find(contains( para_all_label, para_exclude));
disp(para_change_exclude_idx)
para_fixed = [1.1475 , 1.1489 , 1.1489 ];

num_para_change = length(para_change);
para_change_idx =  find(contains( para_all_label, para_change));

%%
A = fnc_getSobolSetMatlab(num_para_change , num_sample);%lhsdesign(num_sample,num_para_change);
%B = fnc_getSobolSetMatlab(num_para_change , num_sample); %lhsdesign(num_sample,num_para_change);

model_obj = create_model(Model_Index, [], outputfile_dir );
t_span = model_obj.t_span_AfterTreatment;

%%

mu_para_change = Para_Value(para_change_idx);

% ub_col = zeros(1, num_para_change);
% lb_col = zeros(1, num_para_change);

if Model_Index == 1
    lb_col =  [
        1e-6*ones(9,1)
        0.1*ones(9,1)
        0.05*ones(3,1)
        0.1*ones(2,1) 
        ]';
    
    ub_col = [
        ones(9,1)
        100*ones(7,1)
        25*ones(2,1)
        5
        5
        5
        10*ones(2,1)  
        ]';
elseif Model_Index == 41
    lb_col =  [
              1e-6*ones(9,1)
              0.1*ones(9,1)
              0.1
              50
              50
              0.1*ones(2,1) 
            ]';
            
    ub_col = [
            ones(9,1)
            100*ones(7,1)
            25*ones(2,1)
            5 %n_i
            700
            700
            20*ones(2,1)  
     ]';

elseif Model_Index == 42

    lb_col =  [
        0
        0
        0
        0
        0
        0
        0
        0.1
        0.1
        0.1
        0.1
        0.1
        0.05
        0.05
        0.05
        4
        4
        0.1
        0.1
        ]';
    
    ub_col = [
        1
        1
        1
        1
        1
        1
        1
        100
        100
        100
        100
        500
        30
        30  
        5
        200
        200
        10*ones(2,1)
     ]';

end

L_A =  repmat(lb_col, num_sample,1)  +  A.*repmat(ub_col- lb_col, num_sample,1 )  ;
%L_B =  repmat(lb_col, num_sample,1)  +  B.*repmat(ub_col- lb_col, num_sample,1 )  ;
% Modify L_A and L_B to make sure that the sum of the first three elements is <= 1
for i = 1:num_sample
    sum_first_three_A = sum(L_A(i, 2:3));

    if sum_first_three_A > 1
        L_A(i, 2:3) = L_A(i, 2:3) / sum_first_three_A;
    end

end


M_A =  cell(num_sample,1); % matrix A 
M_A_N_d = cell(num_sample,1); % matrix A 



parpool('local', num_worker)
parfor j = 1:num_sample
    tic; % Start timing
    task = getCurrentTask;
    taskid = task.ID;
    outfile = fopen( fullfile(outputfile_dir,  ['logfile_' num2str(taskid) '.log'] ), 'a');

    local_model_obj = model_obj;
    theta_A = L_A(j,:);
    theta_A(para_change_exclude_idx) = para_fixed;
   
    local_model_obj.para_unknown =  theta_A;
    [~, ~, CellFractions, N_d,~, ~,   ~ , ~,~]  =   local_model_obj.model_simulation(  );
                           
  
    M_A{j}=    CellFractions ;
    M_A_N_d{j}  = N_d;
    %M_B_G2M(j,:) =   SimData_B(:,3)';
    
    elapsedTime = toc; % End timing and get elapsed time
    fprintf(outfile, 'successful : %d/ %d | Time taken: %.2f seconds \n',j, num_sample,  elapsedTime);
   

end
delete(gcp('nocreate'))


save(fullfile(outputfile_dir , 'predictive_check.mat' ),  'L_A', 'M_A', 'M_A_N_d', 't_span',  '-v7.3');




