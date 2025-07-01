
close all
Cotreatment = 1;

plotflag = 1;

order = 2;
switch order
    case 1
        CellLine_Model_Index =1 ;
        para_unknown = [
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
            0.7950
            0.2350
            0.7120
            0.1410
            0.5220
            0.5200
            0.3440
            25.0410
            14.9820
            46.1890
            11.1720
            83.1670
            12.6860
            19.4550
            2.7800
            169.9660
            64.4420
            4.1970
            3.8810
            ]   ;
        
        
        Dose_1thDrug_col =  [ 0.5, 1.3, 2.6, 5, 10, 20 ]; %G2M phase drug, Dtx
        Dose_2thDrug_col  =  [3.3, 33, 170, 330, 1700, 3300]; %Gem
        numdrugcombo = length(Dose_1thDrug_col );
    case 2
        CellLine_Model_Index = 41 ;
        % para_unknown =    ;
        para_unknown =  [
            0.5920
            0.0260
            0.0710
            0.0150
            0.6280
            0.1340
            0.9330
            0.2890
            0.8170
            15.8860
            38.2940
            93.8220
            59.7440
            94.9020
            91.4040
            97.2710
            3.0640
            11.7460
            0.6350
            141.3970
            225.2360
            0.6190
            5.0030
            0.7950
            0.2350
            0.7120
            0.1410
            0.5220
            0.5200
            0.3440
            25.0410
            14.9820
            46.1890
            11.1720
            83.1670
            12.6860
            19.4550
            2.7800
            169.9660
            64.4420
            4.1970
            3.8810
            
            
            
            ];
        Dose_1thDrug_col =  [ 3.3, 33, 170, 330, 1700, 3300 ]; %G2M phase drug, Pac
        Dose_2thDrug_col  =  [10, 31, 63, 122, 1000, 10000 ] ; %Gem
        numdrugcombo = length(Dose_1thDrug_col );
end
if plotflag == 0
    if Cotreatment == 1
        
        
        output_dir =  ['Output_Cotreatment', num2str(order) ,'_' datestr(now,'yyyy_mm_dd_HH_MM')];
        if ~isfolder(   output_dir )
            mkdir(output_dir)
        end
        
        CellFractions_col = cell(numdrugcombo*numdrugcombo,1);
        for i =  1 : numdrugcombo
            for j =1 : numdrugcombo
                Dose_1thDrug =   Dose_1thDrug_col(i);
                Dose_2thDrug =   Dose_2thDrug_col(j);
                model_obj = create_model_CombinationTreatment(CellLine_Model_Index, para_unknown, output_dir, Dose_1thDrug,Dose_2thDrug  );
                [ CellFractions, N_d, t_span_new] =  model_obj.M_matrix_plot()  ;
                CellFractions_col{ (i-1)*numdrugcombo+j } = CellFractions;
            end
        end
        save(fullfile(  output_dir, 'Cotreatment.mat'),  'CellFractions_col', 't_span_new')
        %elseif  Cotreatment == 0
    end
    
elseif plotflag ==1
    output_dir = 'Output_Cotreatment2_2023_11_23_17_04';
    load(fullfile(  output_dir, 'Cotreatment.mat'),  'CellFractions_col', 'N_d_col',  't_span_new')
    
    f = figure('Position', [462,211,1461,1055]);
    tcl = tiledlayout(numdrugcombo, numdrugcombo); % 1 row, num_Dose columns
    tclcol = [];
    colorMatrix = [
        %[0, 107, 164]/255; % Blue
        [255, 128, 14]/255; % Orange
        [44, 160, 44]/255; % Green
        [148, 103, 189]/255; % Purple
        ];
    lineStyles = {'-', '--', ':'};
    
    
    for i =   1: numdrugcombo
        for j = 1: numdrugcombo
            Dose_1thDrug =   Dose_1thDrug_col(i);
               Dose_2thDrug =   Dose_2thDrug_col(j); 
                t = nexttile;
                tclcol = [t, tclcol];
                current_CellFractions = CellFractions_col{ (i-1)*numdrugcombo+j };
                current_N_d = N_d_col{ (i-1)*numdrugcombo+j  }; 
                for p = 1:3
                    plot(t_span_new, current_CellFractions(p, :), 'Color', colorMatrix(p, :), 'LineStyle', lineStyles{p}, 'LineWidth', 2.5);
                    hold on
                end
               
                % Store the position of the subplot
                subplotPosition = get( t, 'Position');
                subplot_pos =  subplotPosition;
                % Create an inset axes within the current subplot
                % The position vector is [left bottom width height]
                % Define the position for the inset axes based on the subplot
                insetWidth = 0.25 * subplotPosition(3); % 25% of the subplot's width
                insetHeight = 0.25 * subplotPosition(4); % 25% of the subplot's height
                insetLeft = subplotPosition(1)+0.015 ;%subplotPosition(1)*0.1;
                insetBottom = subplotPosition(2) + subplotPosition(4) - insetHeight;
                insetPosition = [insetLeft insetBottom insetWidth insetHeight];
                % Create an inset axes at the calculated position
                
                inset_axes = axes('Position',insetPosition);
                % create smaller axes in top right, and plot on it
              
                plot( t_span_new  , current_N_d, 'LineWidth', 2 )
                ylim([0 1])
                axis tight;
                axes(t);

                xticks(0:20:t_span_new(end));
                title(sprintf('Pac %d: %g nM, Gem %d: %g nM', i, Dose_1thDrug ,j,  Dose_2thDrug ), 'FontSize', 14);
                %title(sprintf('Dose %d: %g nM', i, Dose_col(i)), 'FontSize', 14);
                set(gca, 'FontSize', 12);
                set(gca, 'FontName', 'Arial')
                grid on;
                box off;
                
        end
    end
    linkaxes(tclcol, 'xy');
    ax1 = tclcol(1);
    ax1.YLim = [0 1.2];
    
    % Add a legend
    hL = legend({'G1', 'S', 'G2M'}, 'FontSize', 20);
    % Move the legend to the right side of the figure
    hL.Layout.Tile = 'East';
    
    xlabel(tcl, 'Hours', 'FontSize', 18);
    ylabel(tcl, 'Cell Fraction', 'FontSize', 18);
    set(0, 'DefaultAxesFontName', 'Arial');
    exportgraphics( f ,  fullfile( output_dir  , 'Cotreatment.png' ) ,'Resolution' , 1200   );

end
