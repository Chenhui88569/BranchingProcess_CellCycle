clear
close all

Model_Index = 42;
%dir = 'Output_CellCycle_Model42_2023_10_15_14_35'; 
switch Model_Index
    case 1
        percent  = 0.1;
        dir = 'Output_CellCycle_Model1_2023_10_09_13_13';
        str = sprintf('Para_treated_%d.mat', Model_Index);
        load( fullfile( dir, str ));
        Dose_col = [0.001, 0.01, 0.1, 1, 10, 20];
        num_Dose = length(Dose_col);
    case 41
     
    case 42
        percent  = 0.2;
        dir = 'Output_CellCycle_Model42_2023_10_15_14_35'; 
        str = sprintf('Para_treated_%d.mat', Model_Index);
        load( fullfile( dir, str ));
        Dose_col = [5, 10,  100, 1000, 10000];
        num_Dose = length(Dose_col);
end


Para_col_for_analysis = Para_col_accepted;
num_sample = size(Para_col_for_analysis,1);
Para_trun = Para_col_for_analysis( num_sample*percent +1:end    ,:);
Para_trun_uniq = unique( Para_trun,'rows' );
new_num_sample = size(Para_trun_uniq,1);
para_mean = median(Para_trun, 1);


outputfile_dir = ['DoseSimulation_', dir  ];
if ~isfolder(outputfile_dir)
    mkdir(outputfile_dir)
else 
    rmdir(outputfile_dir,'s');
     mkdir(outputfile_dir)
end



model_obj = create_model(Model_Index, [], outputfile_dir );
model_obj.para_unknown = para_mean ;
%[exception , SimData,CellFractions, N_d , t_span_new , ~,   ~ , ~,~]  =   model_obj.model_simulation(  );

Dose_simulation = cell(num_Dose,1); 

f = figure('Position', [291,100,800,1237]);
tcl = tiledlayout('flow');
% Define colors and line styles
colorMatrix = [[0 0.4470 0.7410]; [0.8500 0.3250 0.0980]; [0.4940 0.1840 0.5560]];
lineStyles = {'-', '--', ':'};
last_precentage = zeros(num_Dose,3 );


for i = 1:num_Dose     
    fprintf('Complete %d',i)
    model_obj.Dose = Dose_col(i);
    if i ==1 && Model_Index == 42
        model_obj.StepSize = 0.02;
        t_span_full = 0: model_obj.StepSize: model_obj.MaxT;
        if t_span_full(end) ~=  model_obj.MaxT
            t_span_full  = [t_span_full model_obj.MaxT];
        end
        model_obj.t_span = t_span_full ;
        model_obj.t_span_AfterTreatment = 0: model_obj.StepSize: model_obj.data_time(end) + 4 ;
        model_obj.num_t = length(t_span_full );
    else
        model_obj.StepSize = 0.1;
        t_span_full = 0: model_obj.StepSize: model_obj.MaxT;
        if t_span_full(end) ~=  model_obj.MaxT
            t_span_full  = [t_span_full model_obj.MaxT];
        end
        model_obj.t_span = t_span_full ;
        model_obj.t_span_AfterTreatment = 0: model_obj.StepSize: model_obj.data_time(end) + 4 ;
        model_obj.num_t = length(t_span_full );
    end
    
    [exception , SimData,CellFractions, N_d , t_span_new , ~,   ~ , ~,~]  =   model_obj.model_simulation(  );
    Dose_simulation{i} =  CellFractions;
    
    current_CellFraction = Dose_simulation{i};
    nexttile(tcl)
    % Plot simulation lines with different line styles
    for p = 1:3
        plot(t_span_new, current_CellFraction(p, :), 'Color', colorMatrix(p, :), 'LineStyle', lineStyles{p}, 'LineWidth', 2.5);
        hold on
        last_precentage(i,p) = current_CellFraction(p, end);
    end
    
    % Set labels and adjust plot properties
    title(sprintf('Dose %d: %g nM', i, Dose_col(i)  ), 'FontSize', 20);
    set(gca, 'FontSize', 12);
    grid on;
    box off
    % Customize other plot properties as needed
end


% Add a legend
hL = legend({'G1M', 'S', 'G2M'}, 'FontSize', 20);
% Move the legend to the right side of the figure
hL.Layout.Tile = 'East';
xlabel(tcl, 'Hours', 'FontSize', 18);
ylabel(tcl, 'Cell Fraction', 'FontSize', 18);

exportgraphics( f ,  fullfile( outputfile_dir, ['Figs_Model' num2str(Model_Index) '.png' ]) ,'Resolution' , 800   );
table_data = table(     Dose_col (:), last_precentage(:,1),  last_precentage(:,2) ,  last_precentage(:,3) );
table_data.Properties.VariableNames = {'Dose level', 'G1', 'S', 'G2M'};
save(fullfile( outputfile_dir, 'Dose_simulation.mat'),  'Dose_simulation','table_data')


    



