% Global sensitivity analysis
% Model output is the kinectics of the tumor proliferating cells + damaging cells(100 mg/kg weekly )
% step1. choose the most influential parameter
% step2. Variance-based sensitivity analysis
% purturb the parameter with the range (determine the range)
% create 1000 parameter sets using Latin hypercube sample
close all
clear 
%dbstop if error


num_worker = 18;
num_sample =  250; %Total number of samples is 500
Model_Index = 42;


outputfile_dir = ['GSA_CellCycle_Model', num2str(Model_Index) ,'_' datestr(now,'yyyy_mm_dd_HH_MM')];
if ~isfolder(outputfile_dir)
    mkdir(outputfile_dir)
else 
    rmdir(outputfile_dir,'s');
     mkdir(outputfile_dir)
end


addpath('ReadYAML-master/ReadYAML-master/')
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


%% load selected parameters
para_change=  para_all_label(1:end-3);


num_para_change = length(para_change);
para_change_idx =  find(contains( para_all_label, para_change));

%%
A = fnc_getSobolSetMatlab(num_para_change , num_sample);%lhsdesign(num_sample,num_para_change);
B = fnc_getSobolSetMatlab(num_para_change , num_sample); %lhsdesign(num_sample,num_para_change);

model_obj = create_model(Model_Index, [], outputfile_dir );
t_span = model_obj.data_time; 
T_num = length(t_span);

%%

mu_para_change = Para_Value(para_change_idx);

% ub_col = zeros(1, num_para_change);
% lb_col = zeros(1, num_para_change);

lb_col =  [
          1e-2*ones(7,1)
          Para_Value(8:end)'*0.1
        ]';
        
ub_col = [
        0.9*ones(7,1)
        Para_Value(8:end)'*10
 ]';

L_A =  repmat(lb_col, num_sample,1)  +  A.*repmat(ub_col- lb_col, num_sample,1 )  ;
L_B =  repmat(lb_col, num_sample,1)  +  B.*repmat(ub_col- lb_col, num_sample,1 )  ;
% Modify L_A and L_B to make sure that the sum of the first three elements is <= 1
for i = 1:num_sample
    sum_first_three_A = sum(L_A(i, 2:3));
    sum_first_three_B = sum(L_B(i, 2:3));

    if sum_first_three_A > 1
        L_A(i, 2:3) = L_A(i, 2:3) / sum_first_three_A;
    end

    if sum_first_three_B > 1
        L_B(i, 2:3) = L_B(i, 2:3) / sum_first_three_B;
    end
end


%loop over the parameters
%[Tsim, Ysim,f_A_prime] =  kinetics_TGI(mu,1,T);


M_A_G2M = zeros(num_sample,T_num); % matrix A 
M_B_G2M = zeros(num_sample,T_num); % matrix B



parpool('local', num_worker)
parfor j = 1:num_sample
    task = getCurrentTask;
    taskid = task.ID;
    outfile = fopen( fullfile(outputfile_dir,  ['logfile_' num2str(taskid) '.log'] ), 'a');

    local_model_obj = model_obj;
    theta_A = L_A(j,:);
    theta_A_complete = Para_Value;
    theta_A_complete_var = theta_A_complete ;   theta_A_complete_var(para_change_idx) = theta_A;
    theta_B = L_B(j,:);
    theta_B_complete = Para_Value;
    theta_B_complete_var = theta_B_complete ;   theta_B_complete_var(para_change_idx) = theta_B;
    local_model_obj.para_unknown =  theta_A_complete_var;
    [~, SimData_A, ~, ~,~, ~,   ~ , ~,~]  =   local_model_obj.model_simulation(  );
                           
    local_model_obj.para_unknown =  theta_B_complete_var ;
    [~, SimData_B, ~, ~,~, ~,   ~ , ~,~]   =   local_model_obj.model_simulation(  );
    
    %disp(size(SimData_A));

    M_A_G2M(j,:) =   SimData_A(:,3)';
    M_B_G2M(j,:) =   SimData_B(:,3)';
    fprintf(outfile,' loop1 : successful : %d/ %d\n', j,num_sample );

end
delete(gcp('nocreate'))

%%

M_A_prime_col_G2M = cell( num_sample ,1);

parpool('local', num_worker)
parfor j = 1 : num_sample  
    task = getCurrentTask;
    taskid = task.ID;
    outfile = fopen( fullfile(outputfile_dir,  ['logfile_' num2str(taskid) '.log'] ), 'a');
    fprintf(outfile,'loop2 [%d/ %d] \n', j,num_sample );

    M_A_prime_G2M = zeros(num_para_change, T_num);
    for i = 1: num_para_change
        % To prevent the error 'unclassified variable.', create a temporary array in each paraloop and then fill in in the nested for loop
        theta_A = L_A(j,:);
        theta_A_complete = Para_Value;
        theta_A_complete_var = theta_A_complete ;   theta_A_complete_var(para_change_idx) = theta_A;
        theta_B = L_B(j,:);
        theta_B_complete = Para_Value;
        theta_B_complete_var = theta_B_complete ;   theta_B_complete_var(para_change_idx) = theta_B;
        para_id = para_change_idx(i);
        if    para_id   == 1
            theta_A_prime = [  theta_B_complete_var( para_id )    theta_A_complete_var(2:end)];
        else
            theta_A_prime = [theta_A_complete_var(1: para_id-1)  theta_B_complete_var( para_id)   theta_A_complete_var( para_id+1:end)];
        end
        local_model_obj = model_obj;
        local_model_obj.para_unknown = theta_A_prime ;
         [~, SimData_A_prime, ~, ~,~, ~,   ~ , ~,~]   =   local_model_obj.model_simulation(  );

        M_A_prime_G2M(i,:) = SimData_A_prime(:,3)';
        fprintf(outfile,'\t loop2 : success for para : %d/ %d\n',i,num_para_change );
    end
    M_A_prime_col_G2M{j} = M_A_prime_G2M;
    fprintf(outfile,'loop2 : successful : %d/ %d\n', j,num_sample );
    fprintf(outfile,'--------------------------------------------\n' );

end
save(fullfile(outputfile_dir , 'SobolInfo.mat' ),  'M_B_G2M', 'M_A_G2M', 'M_A_prime_col_G2M', '-v7.3');
delete(gcp('nocreate'))



