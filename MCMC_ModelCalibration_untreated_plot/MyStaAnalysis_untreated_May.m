close all
clear
Model_Index = 1;
addpath('/Users/chenhuima/MatlabWorkShop/ReadYAML-master/ReadYAML-master/')

switch Model_Index
    case 1
        %load(['Para_TGI_col_Model' num2str(Model_Index) '.mat'])
        load('Output_CellCycle_Model1_2023_04_05_09_43/Para_untreated_1.mat')
        priorcol = {'Weibull' , 'Weibull',  'Weibull', 'Gamma', 'Gamma', 'Gamma',  'Gamma' } ;
        priorpara_a = [ 9
            9
            9
            1
            1
            1
            1];
        priorpara_b = [ 2
            2
            2
            1
            1
            1
            1];
    case 3
        load('Output_CellCycle_Model3_2023_05_17_08_19/Para_untreated_3.mat')
    case 4
        priorcol = {'Weibull' , 'Weibull',  'Weibull', 'Gamma', 'Gamma', 'Gamma',  'Gamma' } ;
        load('Output_CellCycle_Model4_2023_05_15_19_53/Para_untreated_4.mat')
        priorpara_a = [ 10
            10
            10
            2
            1
            1
            1];
        priorpara_b = [ 3
            3
            3
            1
            0.1
            0.1
            0.1];
end
%portion_coll = [0.2884;0; 0.2 ;0.6];
portion_coll = [0.2884;0; 0.5;0.45];
Para_col_for_analysis = Para_col;

nonzeroRows = all(Para_col_for_analysis,2);
Para_col_for_analysis = Para_col_for_analysis(nonzeroRows, :);  


outputfile_dir = ['Figs_', num2str(Model_Index )];
if ~isfolder(outputfile_dir)
    mkdir(outputfile_dir)
end

para_label = ["$\hat{\alpha}_1$"   
    "$\hat{\alpha}_2$" 
     "$\hat{\alpha}_3$"   
    "$\hat{\beta}_0$"
    "$\sigma_{1,  untreated }^2$"
    "$\sigma_{2,  untreated }^2$"
    "$\sigma_{3,  untreated}^2$"
];

percent = portion_coll(Model_Index);%0.3 96
num_sample = size(Para_col_for_analysis,1);
num_para = length(para_label);
Para_trun = Para_col_for_analysis( floor(num_sample*percent)+1:end,1:num_para);
Para_mean_col = zeros(num_para,1);

%cmap = colormap('cool');
binwidth = 0.5;
alpha = 0.5;

factor = 9;
T_col =  cell(1,num_para);
for iFig = 1: ceil(num_para/factor )
    f= figure;
    set(0, 'CurrentFigure', f);
    set(0, 'DefaultFigurePosition', [1440,728,679,609]);
    for i = 1:factor 
        if ismember(i,[5,6,7])
            binwidth = 0.01;
        end
        subplot(3,3, i)
        idx = factor*(iFig-1)+i;
        dbstop if error
        histogram(Para_trun(:,idx  ), 'FaceColor', [0.7 0.7 0.7], 'FaceAlpha', alpha, 'Normalization', 'pdf' ) %'BinWidth', binwidth,
        hold on
                % Calculate the x values to cover the data range
        xi = linspace(min(Para_trun(:, idx)), max(Para_trun(:, idx)), 500);
        
        ChoiceOfDist =  priorcol{idx}; 
        curr_a = priorpara_a(idx);
        curr_b = priorpara_b(idx);
        pd = makedist(ChoiceOfDist,'A',curr_a,'B',curr_b);%distribution object
        % Calculate Weibull pdf
        yi = pdf(pd, xi);
        
        % Plot Weibull pdf
        plot(xi, yi, 'k-', 'LineWidth', 2);

        
        
       
        [counts, edges] = histcounts(Para_trun(:,idx ),'BinWidth', binwidth);
        locs = movmean(edges, 2, 'Endpoints', 'discard');
        %plot(locs, counts, 'LineWidth', 3);
        ax = gca;
        xlabel('Sampled Values', 'FontSize', 10);
        ylabel('Frequency', 'FontSize', 10);
        title(para_label(idx ),'FontSize', 15, 'Interpreter','Latex')
        Para_mean = median(Para_trun(:,idx ));
        %Para_mean = round(Para_mean);
        Para_mean_col(idx ) =  Para_mean;
        hold on
        xline(Para_mean,'r', 'LineWidth',3)
        x = Para_trun(:,idx );
        CI  = quantile(x, [0.25, 0.75]);        %CI = round(CI,3);
        %bar_label = strcat( ' bar', '{', para_label(i),'}'   );
        bar_label = strcat( '(',para_label(idx ),')','_{avg}' );
        temp = strcat(num2str( Para_mean) ,'(', num2str(CI(1)), ',', num2str(CI(2)),')');
        T_col{idx} = temp;

        %text( ax.XTick(1),  (ax.YTick(end)+ax.YTick(end-1))/2  ,temp,'tex'  )
        hold off
        if  idx == num_para
            break
        end
    end
    exportgraphics( f ,  fullfile(outputfile_dir, ['hit', num2str( iFig), '.png'] ) ,'Resolution' , 1200   );

end

T = table(T_col);

for  i = 1: floor(num_para/6)
       f =  figure;
       set(0, 'CurrentFigure', f);
      % set(f, 'Position', get(0, 'Screensize'));
        for j = 1:6
            subplot(2,3, j)
            idx = 6*(i-1)+j;
            plot(Para_trun(:, idx ),'k-','LineWidth',1.5)
            xlabel('iteration')
            %axis([1+T*fix(curr_T/T),T+T*fix(curr_T/T),0,15]);
            grid off
            title(para_label(idx ),'FontSize', 15, 'Interpreter','Latex')
            if  idx == num_para
                break
            end
        end
end


Data = [39.819277108433674      35.060240963855364    25.240963855421647
    0 0 0 
    55.89041095890411    24.883720930232457   19.310344827586523
    40.5645161290322    31.209677419354865    28.70967741935491 ] *0.01   ;
Data = Data - (sum(Data,2)-1)/3;

SimulatedPercentage = fun_untreated_lambdaMalObj(Para_mean_col);
TruePercentage = Data(Model_Index,:)';
std = zeros(3,1);
mean = zeros(3,1);
ChoiceOfDist = 'Gamma';
for i = 1:3
    a =  Para_mean_col(i);
    b = 1/Para_mean_col(4);
    [s, m]  = PdfcdfCal_MeanStd(ChoiceOfDist, a ,b);
     std (i) = s;
     mean (i) = m;
     
end
T_durations = table(std, mean, SimulatedPercentage,TruePercentage )

theta =  Para_mean_col;
[SteadyState,  ClinicalTime, CellFractions , t_span] = MCMC_myfun_untreated_Multivariate_Pgf_v2(theta);
ClinicalTime = t_span(ClinicalTime);
f2 = figure;
plot( t_span, CellFractions' ,'LineWidth', 1.5 )
legend({'G1', 'S', 'G2M'}, 'Location', 'best')
yline(   Data(Model_Index,:), '--k', 'LineWidth', 3  )

xline( ClinicalTime ,'--', ['treatment time:'  num2str(ClinicalTime ),'h']  ,'LineWidth',3,'LabelOrientation','horizontal' );
exportgraphics( f2 ,  fullfile(outputfile_dir,'Simulation.png' ) ,'Resolution' , 300   );
