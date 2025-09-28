% For individual TGl
close all
clear
delete(gcp('nocreate'))
Model_Index = 3;
switch Model_Index
    case 1
        %load(['Para_TGI_col_Model' num2str(Model_Index) '.mat'])
        load('Output_CellCycle_Model1_2023_05_17_08_14/Para_untreated_1.mat')
    case 3
        load('Output_CellCycle_Model3_2023_05_17_08_19/Para_untreated_3.mat')
    case 4
        load('Output_CellCycle_Model4_2023_05_15_19_53/Para_untreated_4.mat')
end
%portion_coll = [0.2884;0; 0.2 ;0.6];
portion_coll = [0.9;0; 0.5;0.45]; %0.5
Para_col_for_analysis = Para_col;

nonzeroRows = all(Para_col_for_analysis,2);
Para_col_for_analysis = Para_col_for_analysis(nonzeroRows, :);  


outputfile_dir= ['CI_CellCycle_Model', num2str(Model_Index) ,'_' datestr(now,'yyyy_mm_dd_HH_MM')];
if ~isfolder(outputfile_dir)
    mkdir(outputfile_dir)
else 
    rmdir(outputfile_dir,'s');
     mkdir(outputfile_dir)
end

percent = portion_coll(Model_Index);%0.3 96
num_sample = size(Para_col_for_analysis,1);
num_para =  size(Para_col_for_analysis,2);
para_trun = Para_col_for_analysis( floor(num_sample*percent)+1:end,1:num_para);
num_sample_trun = size(para_trun,1 );
Para_mean = mean(para_trun , 1);
theta =  Para_mean;
[SteadyState,  ClinicalTime, MeanCellFractionsTimeCourse , t_span] = MCMC_myfun_untreated_Multivariate_Pgf_v2(theta);

TimeCourse_col_G1 = zeros(num_sample_trun, length( t_span ));
TimeCourse_col_S = zeros(num_sample_trun, length( t_span ));
TimeCourse_col_G2M = zeros(num_sample_trun, length( t_span ));

parpool(20)
parfor i = 1:num_sample_trun
    task = getCurrentTask;
    taskid = task.ID;
    outfile = fopen( fullfile(outputfile_dir,  ['logfile_' num2str(taskid) '.log'] ), 'a');
    
    CurrPara = para_trun(i,:);
    [SteadyState,  ClinicalTime, CurrCellFractionsTimeCourse , t_span] = MCMC_myfun_untreated_Multivariate_Pgf_v2(CurrPara);
    TimeCourse_col_G1(i,:) = CurrCellFractionsTimeCourse(1,:);
    TimeCourse_col_S(i,:) = CurrCellFractionsTimeCourse(2,:);
    TimeCourse_col_G2M(i,:) = CurrCellFractionsTimeCourse(3,:);
    fprintf(outfile,'successful : %d\n', i);
end

delete(gcp('nocreate'))
save( fullfile( outputfile_dir ,  strcat('CI_untreated_', num2str(Model_Index),  '.mat' ) )  ,'TimeCourse_col_G1'  , 'TimeCourse_col_S',  'TimeCourse_col_G2M',  'Para_mean', 'MeanCellFractionsTimeCourse','t_span' ,'-v7.3');


