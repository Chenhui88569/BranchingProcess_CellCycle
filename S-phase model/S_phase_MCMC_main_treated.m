clear
close all
T = 2000;

init_para = [   
    0.7500
    0.0319
    0.8570
    0.4335
    0.2615
    0.6946
    0.3206
   51.7330
   10.3197
   42.4213
   30.3945
   95.4084
   15.5446
    4.0785
    2.4288
  102.5251
   38.9257
    3
   1
   8e-4
   8e-4
   8e-4
];

lb =  [
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
    7e-4*ones(3,1)
    ];

ub = [
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
    100
    25
    25
    5
    200
    200
    10*ones(2,1)
    9e-4*ones(3,1)
 ];

init_para = init_para';
num_para = length(init_para);

para_label =["$q_{1, S}$"
"$q_{2, S}$"
"$q_{3, S}$"
"$m_{d,FR, S}^{max}$"
"$m_{d,UR, S}^{max}$" 
"$m_{d,UR,G2M}^{max}$"
"$m_{d,G1,block}$"
"$\tilde{\lambda}_{d,max}$"
"$\tilde{\lambda}_{Sarrest, UR, max}$"
"$\tilde{\lambda}_{Sarrest, FR, max}$"
"$\tilde{\lambda}_{G2Marrest, max}$"
"$\tilde{\lambda}_{G1block, max}$"
"$\tilde{b}$"
"$\tilde{b}_d$"
"$n_{S}$"
"$EC_{50,d, S}$"
"$EC_{50, S}$"
"$a_G2M$"
"$b_G2M$"
"$\theta_1$"
"$\theta_2$"
"$\theta_3$"
];



Model_Index =  42; 

history_flag = 0;
if  history_flag 
    history_folder  = get_latest_batch_file(Model_Index)
    fprintf('history_folder %s',  history_folder)
    load( fullfile(history_folder,  sprintf('Para_treated_%d.mat', Model_Index) ))
end
outputfile_dir= ['Model_Calibration_Results/'  'Output_CellCycle_Model', num2str(Model_Index) ,'_' datestr(now,'yyyy_mm_dd_HH_MM')];


para_exclude = ["$n_{S}$","$EC_{50,d, S}$","$EC_{50, S}$"];
para_change_exclude_idx = find(contains( para_label, para_exclude));
para_fixed = [1 ,10.066, 26.073 ];

if ~isfolder(outputfile_dir)
    mkdir(outputfile_dir)
end
    
outfile = fopen(  fullfile(outputfile_dir,  'outlog.log' ),'a');

%Output_CellCycle_Model41_2023_12_05_10_30
%load('Output_CellCycle_Model1_2024_10_24_18_52/Para_treated_1.mat')
%Output_CellCycle_Model1_2023_11_26_13_2
if exist('Para_col_accepted','var')

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
    
    curr_T = num_HistPoints
        
else
    fprintf(outfile, 'start from 0')
    num_HistPoints = 0;

    num_para = length(init_para);

    Para_col = zeros( T,    num_para - length(para_fixed) ) ;
    init_para_exclude = init_para ;
    init_para_exclude(para_change_exclude_idx) = [];
    disp(para_change_exclude_idx)
    Para_col(1,:) = init_para_exclude ; %Para_col_hist(end,:);
        
    Para_col_all_new = zeros(size(Para_col ));
    Diff_norm_col_new = zeros(T,1);
    curr_T =  1
end

%a lower diagonal matrix with positive diagonal elements
if ~exist('Para_col_all_hist','var')
    S_curr = zeros(num_para);
    %subset1 = find(contains( para_label    ,[""$"d_{G1}$""" ""$"d_S"$""  ""$"d_{G2M}$""" ""$"p1_{G1}$"""  ""$"p1_{S}$"""  ""$"p3_{G2M}$"""  ""$"p2_{G2M}$"""  ""$"p2_{G1}$""" ""$"p2_{S}$"""  ""$"theta_1"$""  ""$"theta_2"$""  ""$"theta_3"$"" ] ));
    %subset1 = find(contains( para_label    ,[""$"d_{G1}$""" ""$"d_S"$""  ""$"d_{G2M}$"""   ""$"p1_{G1}$"""  ""$"p1_{S}$"""  ""$"p3_{G2M}$"""  ""$"p2_{G2M}$"""  ""$"p4_{G2M}$""" ] ));
    subset2 = find(contains( para_label    ,["$q_{1, S}$"
        "$q_{2, S}$"
        "$q_{3, S}$"
        "$m_{d,UR, S}^{max}$"
        "$m_{d,FR, S}^{max}$"
        "$m_{d,UR,G2M}^{max}$"
        "$m_{d,G1,block}$"
        "$n_{S}$"
        ] ));
    subset4 = find(contains( para_label    ,[ "$\theta_1$"
        "$\theta_2$"
        "$\theta_3$"] ));
    % subset3 = 3:8;
    for i = 1: num_para
        S_curr(i,1:i) =  0.5*randn;
        %         if ismember(i, subset1  )
        %             S_curr(i,1:i) =  abs(randn)*0.01;
        %         end
        
         if ismember(i, subset2  )
            S_curr(i,1:i) =  abs(randn)*0.005;
         end
%          if ismember(i,  subset3  )
%             S_curr(i,1:i) =  abs(randn)*0.5;
%         end
        
         if ismember(i,  subset4  )
            S_curr(i,1:i) =  abs(randn)*0.00001;
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
    
myfun = @Target_fun_Treated;

[Para_col_accepted, acc_rate_col, S_curr, Para_col_all, Diff_norm_col] = ...
        MyRAM_treated(...
            T, size(Para_col,2), Para_col, S_curr, myfun, para_label, outputfile_dir , Model_Index, ...
            num_HistPoints,  lb_local, ub_local, para_change_exclude_idx, para_fixed,  outfile, j , curr_T, Diff_norm_col_new, Para_col_all_new );


%if exist('Para_col_accepted_hist','var')
  %  Para_col_accepted= [Para_col_accepted_hist; Para_col_accepted(2:end,:)   ];
 %   Para_col_all = [Para_col_all_hist ; Para_col_all(1:end,:)   ];
%else
  %  Para_col_accepted =Para_col_accepted;
%end

save( fullfile( outputfile_dir ,  strcat('Para_treated_', num2str(Model_Index),  '.mat' ) )  ,'Para_col_accepted'  , 'acc_rate_col',  'S_curr' , 'Para_col_all','Diff_norm_col' ,'-v7.3');
   