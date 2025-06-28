clear
close all
T =  3000;

init_para = [   
   0.5
   0.2
   0.2
   0.8
   0.2
   0.2
   0.8
   30
   60
   60
   84
   30
   10
   10
   0.3
   0.5
   0.5
   1e-3
   1e-3
   1e-3
];

init_para = init_para';
num_para = length(init_para);

para_label =["$q_1$"
"$q_2$"
"$q_3$"
"$m_{d,UR}^{max}$"
"$m_{d,FR}^{max}$"
"$m_{d,MS}^{max}$"
"$m_{d,UR,G1}^{max}$"
"$\tilde{\lambda}_{d,max}$"
"$\tilde{\lambda}_{G2arrest,UR, max}$"
"$\tilde{\lambda}_{G2arrest,FR, max}$"
"$\tilde{\lambda}_{G2arrest,MS, max}$"
"$\tilde{\lambda}_{G1arrest, max}$"
"$\tilde{b}$"
"$\tilde{b}_d$"
"$n_{G2M}$"
"$EC_{50,d}$"
"$EC_{50}$"
"$\theta_1$"
"$\theta_2$"
"$\theta_3$"
];



Model_Index =  1; 
outputfile_dir= ['Output_CellCycle_Model', num2str(Model_Index) ,'_' datestr(now,'yyyy_mm_dd_HH_MM')];


load('Output_CellCycle_Model1_2023_09_12_14_28/Para_treated_1.mat')
if exist('Para_col_accepted','var')
    Para_col_accepted_hist =  Para_col_accepted;
    nonzeroRows = all(Para_col_accepted_hist,2);
    Para_col_accepted_hist = Para_col_accepted_hist(nonzeroRows, :);  
    
    Para_col_all_hist = Para_col_all(nonzeroRows, :);  
     
    num_HistPoints = size(Para_col_accepted_hist ,1);
    fprintf('num of hist points: %d', num_HistPoints)
    Para_col = zeros( T,    num_para ) ;
    Para_col(1,:) =   Para_col_accepted_hist(end,:);
else
    num_HistPoints = [];
    Para_col = zeros( T,    num_para ) ;
    Para_col(1,:) =  init_para ; %Para_col_hist(end,:);
end

%a lower diagonal matrix with positive diagonal elements
if ~exist('Para_col_hist','var')
    S_curr = zeros(num_para);
    %subset1 = find(contains( para_label    ,["$"d_{G1}$"" "$"d_S$""  "$"d_{G2M}$"" "$"p1_{G1}$""  "$"p1_{S}$""  "$"p3_{G2M}$""  "$"p2_{G2M}$""  "$"p2_{G1}$"" "$"p2_{S}$""  "$"theta_1$""  "$"theta_2$""  "$"theta_3$"" ] ));
    %subset1 = find(contains( para_label    ,["$"d_{G1}$"" "$"d_S$""  "$"d_{G2M}$""   "$"p1_{G1}$""  "$"p1_{S}$""  "$"p3_{G2M}$""  "$"p2_{G2M}$""  "$"p4_{G2M}$"" ] ));
    subset2 = find(contains( para_label    ,["$q_1$"
        "$q_2$"
        "$q_3$"
        "$m_{d,UR}^{max}$"
        "$m_{d,FR}^{max}$"
        "$m_{d,MS}^{max}$"
        "$m_{d,UR,G1}^{max}$"
        "$EC_{50,d}$"
        "$EC_{50}$"] ));
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
            S_curr(i,1:i) =  abs(randn)*0.001;
        end
    end
end

S_curr  = S_curr( 1: num_para, 1:num_para );
myfun = @Target_fun_Treated;
[Para_col_accepted,  acc_rate_col,    S_curr, Para_col_all, Diff_norm_col] = MyRAM_treated(  T, num_para, Para_col , S_curr,myfun,para_label ,  outputfile_dir, Model_Index , num_HistPoints );
if exist('Para_col_accepted_hist','var')
    Para_col_accepted = [Para_col_accepted_hist; Para_col_accepted(2:end,:)   ];
    Para_col_all = [Para_col_all_hist ; Para_col_all(1:end,:)   ];
else
    Para_col_accepted =Para_col_accepted;
end
save( fullfile( outputfile_dir ,  strcat('Para_treated_', num2str(Model_Index),  '.mat' ) )  ,'Para_col_accepted'  , 'acc_rate_col',  'S_curr' , 'Para_col_all','Diff_norm_col' ,'-v7.3');

    
   