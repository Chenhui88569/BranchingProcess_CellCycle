clear
close all
T =   2000;

init_para = [   
       0.5420
            0.1040
            0.2060
            0.4700
            0.6430
            0.8000
            0.3040
            0.3700
            0.8530
            36.3680
            41.4260
            10.9410
            8.0180
            59.9780
            37.4030
            47.0860
            6.2900
            12.5460
            1.0660
            0.1700
            0.5200
            0.3550
            1.8310
    6e-4
    6e-4
    6e-4 
];

init_para = init_para';
num_para = length(init_para);

para_label =["$q_1$"
"$q_2$"
"$q_3$"
"$q_4$"
"$m_{d,FR}^{max}$"
"$m_{d,UR}^{max}$"
"$m_{d,MS}^{max}$"
"$m_{d,UR,G1}^{max}$"
"$m_{d,UR,S}^{max}$"
"$\tilde{\lambda}_{d,max,G1}$"
"$\tilde{\lambda}_{G2arrest,UR, max}$"
"$\tilde{\lambda}_{G2arrest,FR, max}$"
"$\tilde{\lambda}_{G2arrest,MS, max}$"
"$\tilde{\lambda}_{G1arrest, max}$"
"$\tilde{\lambda}_{d,max,S}$"
"$\tilde{\lambda}_{Sarrest, max}$"
"$\tilde{b}$"
"$\tilde{b}_d$"
"$n_{G2M}$"
"$EC_{50,d}$"
"$EC_{50}$"
"$a_S$"
"$b_S$"
"$\theta_1$"
"$\theta_2$"
"$\theta_3$"
];




Model_Index = 1;
history_flag = 1;
if  history_flag 
    history_folder  = get_latest_batch_file(Model_Index)
    fprintf('history_folder %s',  history_folder)
    load( fullfile(history_folder,  sprintf('Para_treated_%d.mat', Model_Index) ))
end

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
          7e-4*ones(3,1)
        ];
        
ub = [
        ones(9,1)
        100*ones(7,1)
        25*ones(2,1)
        5 %n_i
        700
        700
        20*ones(2,1)  
        9e-4*ones(3,1)  
 ];
end

outputfile_dir = ['Model_Calibration_Results/' 'Output_CellCycle_Model', num2str(Model_Index) ,'_' datestr(now,'yyyy_mm_dd_HH_MM')];

para_exclude = ["$n_{G2M}$","$EC_{50,d}$","$EC_{50}$"];
para_change_exclude_idx = find(contains( para_label, para_exclude));
para_fixed = [1 ,10.295, 2.301];

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
    Para_col(1,:) = init_para_exclude ; %Para_col_hist(end,:);
        
    Para_col_all_new = zeros(size(Para_col ));
    Diff_norm_col_new = zeros(T,1);
    curr_T =  1
end

%a lower diagonal matrix with positive diagonal elements
if ~exist('Para_col_all_hist','var')
    fprintf('create S_curr');
    S_curr = zeros(num_para);
    %subset1 = find(contains( para_label    ,["$"d_{G1}$"" "$"d_S$""  "$"d_{G2M}$"" "$"p1_{G1}$""  "$"p1_{S}$""  "$"p3_{G2M}$""  "$"p2_{G2M}$""  "$"p2_{G1}$"" "$"p2_{S}$""  "$"theta_1$""  "$"theta_2$""  "$"theta_3$"" ] ));
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
    
    subset4 = find(contains( para_label    ,[ "$\theta_1$"   "$\theta_2$"   "$\theta_3$" ] ));
   % subset3 = 3:8;
    for i = 1: num_para
        S_curr(i,1:i) =  randn;
%         if ismember(i, subset1  )
%             S_curr(i,1:i) =  abs(randn)*0.01;
%         end
        
         if ismember(i, subset2  )
            S_curr(i,1:i) =  abs(randn)*0.03;
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

    
   