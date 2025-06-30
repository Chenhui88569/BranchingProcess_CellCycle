clear
close all

Model_Index = 1;
dir = 'Output_CellCycle_Model1_2025_06_23_16_41';
str = sprintf('Para_treated_%d.mat', Model_Index);
switch Model_Index
    case 1
        DataTime = [ 0
            4
            16
            20
            24
            ];
        Current_Data  = [39.78    35.02    25.20
            28.1102   26.6131   45.2766
            18.1365   20.6562   61.2073
            17.4534   16.9032   65.6434
            14.9088   13.8842   71.2070
            ]   *0.01  ;  % Time * 3
        Current_Data  = Current_Data - (sum(Current_Data ,2)-1)/3;
        load( fullfile( dir, 'Para_treated_1.mat' ));
        para_label =  ["q_{1,G2/M}"
        "q_{2,G2/M}"
        "q_{3,G2/M}"
        "q_{4,G2/M}"
        "m_{d,FR,G2/M}^{max}"
        "m_{d,UR,G2/M}^{max}"
        "m_{d,MS,G2/M}^{max}"
        "m_{d,G1arrest}^{max}"
        "m_{d,Sarrest}^{max}"
        "\lambda_{d,G1}^{max}"
        "\lambda_{UR,G2/M}^{max}"
        "\lambda_{FR,G2/M}^{max}"
        "\lambda_{MS,G2/M}^{max}"
        "\lambda_{G1arrest}^{max}"
        "\lambda_{d,S}^{max}"
        "\lambda_{Sarrest}^{max}"
        "b_{G2/M}"
        "b_{d,G2/M}"
        "n_{G2/M}"
        "EC_{50,d, G2/M}"
        "EC_{50,G2/M}"
         "a_{S}"
        "b_{S}"
        "\sigma_{1, treated}^2"
        "\sigma_{2, treated}^2"
        "\sigma_{3, treated}^2"
        ];
        para_exclude = ["n_{G2/M}", "EC_{50,d, G2/M}","EC_{50,G2/M}"];
        para_fixed = [1 ,10.295, 2.301];
        percent = 0.15; % 0.7;
        Para_col_for_analysis = Para_col_accepted;
    case 41
        DataTime = [ 0
            4
            24
            48
            72
            ];
        Current_Data  = [   40.5645161290322    31.209677419354865    28.70967741935491
            32.5	18.79032258	47.82258065
            16.37096774	12.82258065	70.24193548
            25.32258065	30.88709677	44.43548387
            17.5	51.69354839	29.91935484 ] *0.01;
        Current_Data  = Current_Data - (sum(Current_Data ,2)-1)/3;
        load( fullfile( dir, 'Para_treated_41.mat' ));
        para_label =  ["q_{1,G2M}"
        "q_{2,G2M}"
        "q_{3,G2M}"
        "q_{4,G2M}"
        "m_{d,FR,G2M}^{max}"
        "m_{d,UR,G2M}^{max}"
        "m_{d,MS,G2M}^{max}"
        "m_{d,UR,G1}^{max}"
        "m_{d,UR,S}^{max}"
        "\lambda_{d,G1}^{max}"
        "\lambda_{G2arrest,UR}^{max}"
        "\lambda_{G2arrest,FR}^{max}"
        "\lambda_{G2arrest,MS}^{max}"
        "\lambda_{G1arrest}^{max}"
        "\lambda_{d,S}^{max}"
        "\lambda_{Sarrest}^{max}"
        "b_{G2M}"
        "b_{d,G2M}"
        "n_{G2M}"
        "EC_{50,d, G2M}"
        "EC_{50,G2M}"
         "a_{S}"
        "b_{S}"
        "\theta_1"
        "\theta_2"
        "\theta_3"
        ]; 
        percent  = 0.1;
        Para_col_for_analysis = Para_col_accepted;
    case 42
        DataTime = [ 0
            4
            24
            48
            72
            ];
        Current_Data  = [  40.40322581	31.0483871	28.5483871
            59.43548387	38.9516129	2.5
            63.30645161	29.91935484	4.758064516
            72.5	21.69354839	4.435483871
            86.69354839	10.16129032	3.306451613 ] *0.01;
        Current_Data  = Current_Data - (sum(Current_Data ,2)-1)/3;
        load( fullfile( dir , 'Para_treated_42.mat'));
        percent  = 0.6;
        para_label = ["q_{1, S}"
        "q_{2, S}"
        "q_{3, S}"
        "m_{d,FR, S}^{max}"
        "m_{d,UR, S}^{max}"
        "m_{d,UR,G2M}^{max}"
        "m_{d,G1,block}"
        "\lambda_{d,G2M}^{max}"
        "\lambda_{Sarrest, UR}^{max}"
        "\lambda_{Sarrest, FR}^{max}"
        "\lambda_{G2Marrest}^{max}"
        "\lambda_{G1block}^{max}"
        "b_{S}"
        "b_{d,S}"
        "n_{S}"
        "EC_{50,d, S}"
        "EC_{50, S}"
        "a_{G2M}"
        "b_{G2M}"
        "\theta_1"
        "\theta_2"
        "\theta_3"
        ];
        num_sample = size(Para_col_accepted,1 );
        Para_col_for_analysis = Para_col_accepted( num_sample - 15000: end, :);


        
        
end

para_change_exclude_idx = find(contains( para_label, para_exclude));
para_label(para_change_exclude_idx )  = [];

nonzeroRows = all(Para_col_for_analysis,2);
Para_col_for_analysis = Para_col_for_analysis(nonzeroRows, :);

num_sample = size(Para_col_for_analysis,1);
Para_trun = Para_col_for_analysis( num_sample*percent +1:end    ,:);
Para_trun_uniq = unique( Para_trun,'rows' );
new_num_sample = size(Para_trun_uniq,1);


outputfile_dir = ['CredibleInterval_', dir  ];
if ~isfolder(outputfile_dir)
    mkdir(outputfile_dir)
    CellFraction_col = cell(new_num_sample,1);
    N_d_col = cell(new_num_sample,1);
%else 
    %load( fullfile(outputfile_dir,   strcat('CredibleInterval_', num2str(Model_Index),  '.mat' ) )); 
end



model_obj = create_model(Model_Index, [], outputfile_dir );

%num_worker = 30;
%localCluster = parcluster('local');
%localCluster.NumWorkers = num_worker;
%parpool(localCluster , 30);
for i = 1 :  new_num_sample
    %task = getCurrentTask;
    taskid =1;
    fileid = fopen( fullfile(outputfile_dir, ['log_' num2str(taskid ) '.log']), 'a');
    theta = Para_trun_uniq(i,:);
    theta_full = assemble_theta(theta, para_fixed, para_change_exclude_idx); 
    local_obj = model_obj;
    local_obj.para_unknown =  theta_full ;
    [exception , SimData,CellFractions, N_d ,~, ~,   ~ , ~,~]  =   local_obj.model_simulation(  );
   % [exception,~, CellFractions, N_d, t_span_new]= local_obj.M_matrix_plot(  );
    CellFraction_col{i} =  CellFractions;
    N_d_col{i} = N_d;
    fprintf( fileid, 'Complete for sample %d\n', i);
   % save( fullfile( outputfile_dir ,  strcat('CredibleInterval_', num2str(Model_Index),  '.mat' ) )  , 'CellFraction_col'  , 'N_d_col' ,'-v7.3');
end
save( fullfile( outputfile_dir ,  strcat('CredibleInterval_', num2str(Model_Index),  '.mat' ) )  , 'CellFraction_col'  , 'N_d_col' ,'-v7.3');



