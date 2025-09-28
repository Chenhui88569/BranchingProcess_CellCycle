% For individual TGI
close all
%clear
Model_Index = 1;
addpath('/Users/chenhuima/Matlab_code/BranchingProcess_CellCycle/MCMC_ModelCalibration_treated_model_1_dtx_GM3_new/ReadYAML-master/')
addpath('/Users/chenhuima/Matlab_code/BranchingProcess_CellCycle/MCMC_ModelCalibration_untreated_plot/github_repo/cbrewer2/')
addpath('/Users/chenhuima/Matlab_code/BranchingProcess_CellCycle/MCMC_ModelCalibration_untreated_plot/github_boundedline/boundedline/'  )
addpath('/Users/chenhuima/Matlab_code/BranchingProcess_CellCycle/MCMC_ModelCalibration_untreated_plot/github_boundedline/Inpaint_nans/')

ymlflie = ['Datasetsinfo'  num2str( Model_Index )  '.yml'];
datastr = ReadYaml(ymlflie);
CI_DIR = datastr.CI_DIR;
Data = datastr.data_untreated;
load(  CI_DIR)
endpoint = datastr. untreated_plot_endpoint;


outputfile_dir= ['Plot_CI_CellCycle_Model', num2str(Model_Index)];
if ~isfolder(outputfile_dir)
    mkdir(outputfile_dir)
else 
    rmdir(outputfile_dir,'s');
     mkdir(outputfile_dir)
end


idx = find(t_span>endpoint,1);
ClinicalTime = 250;


% The dimensions of CI are (num_t,2)
CI_G1 = TimeCourses_ConfidenceInterval(TimeCourse_col_G1) ;  CI_G1  = CI_G1(1: idx,:); 
CI_G1_dist =  abs(CI_G1- repmat( MeanCellFractionsTimeCourse(1,1:idx)'  , 1,2) ) ;

CI_S = TimeCourses_ConfidenceInterval(TimeCourse_col_S) ;  CI_S  = CI_S(1: idx,:);
CI_S_dist =  abs( CI_S-  repmat(MeanCellFractionsTimeCourse(2,1:idx)' ,1,2) );

CI_G2M = TimeCourses_ConfidenceInterval(TimeCourse_col_G2M) ; CI_G2M  = CI_G2M(1: idx,:);
CI_G2M_dist = abs( CI_G2M- repmat( MeanCellFractionsTimeCourse(3,1:idx)', 1,2) );
 

CI_col = {CI_G1_dist; CI_S_dist ; CI_G2M_dist    };

f = figure('Units', 'pixels', ...
    'Position', [634,577,930,745]);
hold  on;

colors = cbrewer2('qual', 'Set2', 3);
legnames = {'G1', 'S', 'G2M'};


for i = 1:3
    subplot(3,1,i)  
    [l, bl ] = boundedline(t_span(1:idx), MeanCellFractionsTimeCourse(i,1:idx), CI_col{i} ,...
          'cmap', colors(i,:),  'transparency', 0.2);
    outlinebounds(l,bl); 
    %xline( ClinicalTime ,'--', ['treatment time:'  num2str(ClinicalTime ),'h']  ,'LineWidth',3,'LabelOrientation','horizontal', 'LabelHorizontalAlignment', 'center', 'FontSize' , 15 );
    yline(   Data(i), '-.k', 'Observed steady state' , 'LineWidth', 3 ,'LabelHorizontalAlignment', 'center' ,'FontSize' , 15 )
    drawnow 
    xlabel('Time (h)'); ylabel('Percentage');
    lh = legend(l,'FontSize', 15);
    str = ['\' sprintf('color[rgb]{%f,%f,%f} %s', colors(i, 1), colors(i, 2), colors(i, 3), legnames{i})];
    lh.String = str;
    lh.Box = 'off';
    
    % move a bit closer
    lpos = lh.Position;
    lpos(1) = lpos(1) + 0.1;
    lh.Position = lpos;
    
    
    set(groot, ...
        'DefaultFigureColor', 'w', ...
        'DefaultAxesLineWidth', 0.5, ...
        'DefaultAxesXColor', 'k', ...
        'DefaultAxesYColor', 'k', ...
        'DefaultAxesFontUnits', 'points', ...
        'DefaultAxesFontSize', 15, ...
        'DefaultAxesFontName', 'Helvetica', ...
        'DefaultLineLineWidth', 2, ...
        'DefaultTextFontUnits', 'Points', ...
        'DefaultTextFontSize', 20, ...
        'DefaultTextFontName', 'Helvetica', ...
        'DefaultAxesBox', 'off', ...
        'DefaultAxesTickLength', [0.02 0.025]);
    
    
    set(gca,...
        'Box'         , 'off'     , ...
        'TickDir'     , 'out'     , ...
        'TickLength'  , [.02 .02] , ...
        'XMinorTick'  , 'on'      , ...
        'YMinorTick'  , 'on'      , ...
        'YGrid'       , 'on'      , ...
        'XColor'      , [.3 .3 .3], ...
        'YColor'      , [.3 .3 .3], ...
        'YTick'       , 0:0.3:1, ...
        'XTick'       , 0:20:endpoint, ...
        'Xlim'  , [0, endpoint],...
        'LineWidth'   , 1         );
    set(gcf, 'PaperPositionMode', 'auto');
    
    
end
exportgraphics( f ,  fullfile(outputfile_dir,'Sim.png' ) ,'Resolution' , 300   );
%bl = boundedline(t_span(1:idx), CellFractions(1,1:idx), [.25 .25] , ...
  %  t_span(1:idx), CellFractions(2,1:idx), [.25 .25] , ...
  %  t_span(1:idx), CellFractions(3,1:idx), [.25 .25] , ...
  %  'cmap', colors,  'transparency', 0.2 );


%plot( t_span(1:idx), CellFractions(:,1:idx)' ,'LineWidth', 1.5 )
%legend({'G1', 'S', 'G2M'}, 'Location', 'best')
% instead of a legend, show colored text

% lh = legend;
% legnames = {'G1', 'S', 'G2M'};
% for i = 1:length(legnames)
%     str{i} = ['\' sprintf('color[rgb]{%f,%f,%f} %s', colors(i, 1), colors(i, 2), colors(i, 3), legnames{i})];
% end
% lh.String = str;
% lh.Box = 'off';




