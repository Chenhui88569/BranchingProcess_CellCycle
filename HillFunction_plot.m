clear all
close all
% Parameters
outputfile_dir = "HillFunction_plot";
if ~isfolder(outputfile_dir)
    mkdir(outputfile_dir)
end
n_G2M = 1;               % Hill coefficient
E= @(Dose_col, EC50) ((Dose_col ./ EC50).^n_G2M) ./ (1 + (Dose_col ./ EC50).^n_G2M); 

EC50_G2M = 2.301;
% Given discrete drug concentrations
Dose_col = [0.5, 1.3, 2.6, 5, 10, 20];

% Compute E_G2M at each dose
E_doc = E(Dose_col, EC50_G2M  );
% Plot
f = figure;
plot(Dose_col, E_doc , 'o-', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Drug Concentration $\mathcal{D}_{\mathrm{Docetaxel}}$ (nM)', 'Interpreter', 'latex');
ylabel('$E_{\mathrm{G2/M}}$', 'Interpreter', 'latex');
title('$E_{\mathrm{G2/M}}$ vs. $\mathcal{D}_{\mathrm{Docetaxel}}$', 'Interpreter', 'latex');
set(gca, 'FontSize', 18);
grid on;
exportgraphics( f ,  fullfile( outputfile_dir, ['E_doc' '.png' ]) ,'Resolution' , 1200   );

EC50_G2M_d  = 10.295;
E_d_doc = E(Dose_col, EC50_G2M_d  );
% Plot
f = figure;
plot(Dose_col, E_d_doc , 'o-', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Drug Concentration $\mathcal{D}_{\mathrm{Docetaxel}}$ (nM)', 'Interpreter', 'latex');
ylabel('$E_{d, \mathrm{G2/M}}$', 'Interpreter', 'latex');
title('$E_{d, \mathrm{G2/M}}$ vs. $\mathcal{D}_{\mathrm{Docetaxel}}$', 'Interpreter', 'latex');
set(gca, 'FontSize', 18);
grid on;
exportgraphics( f ,  fullfile( outputfile_dir, ['E_d_doc' '.png' ]) ,'Resolution' , 1200   );


%%%%%%%%%%%%  Paclitaxel
Dose_col =  [ 3.3, 33, 170, 330, 1700, 3300 ];
EC50_G2M = 12.023;
E_pac = E(Dose_col, EC50_G2M  );
% Plot
f = figure;
plot(Dose_col, E_pac , 'o-', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Drug Concentration $\mathcal{D}_{\mathrm{Paclitaxel}}$ (nM)', 'Interpreter', 'latex');
ylabel('$E_{\mathrm{G2/M}}$', 'Interpreter', 'latex');
title('$E_{\mathrm{G2/M}}$ vs. $\mathcal{D}_{\mathrm{Paclitaxel}}$', 'Interpreter', 'latex');
set(gca, 'FontSize', 18);
grid on;
exportgraphics( f ,  fullfile( outputfile_dir, ['E_pac' '.png' ]) ,'Resolution' , 1200   );

EC50_G2M_d  = 5.802;
E_d_pac = E(Dose_col, EC50_G2M_d  );
% Plot
f = figure;
plot(Dose_col, E_d_pac , 'o-', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Drug Concentration $\mathcal{D}_{\mathrm{Paclitaxel}}$ (nM)', 'Interpreter', 'latex');
ylabel('$E_{d, \mathrm{G2/M}}$', 'Interpreter', 'latex');
title('$E_{d, \mathrm{G2/M}}$ vs. $\mathcal{D}_{\mathrm{Paclitaxel}}$', 'Interpreter', 'latex');
set(gca, 'FontSize', 18);
grid on;
exportgraphics( f ,  fullfile( outputfile_dir, ['E_d_pac' '.png' ]) ,'Resolution' , 1200   );



%%%%%%%%%%%%  Gemcitabine
Dose_col = [ 10, 31, 63, 122, 1000, 10000];
EC50_G2M = 26.073;
E_gem = E(Dose_col, EC50_G2M  );
% Plot
f = figure;
plot(Dose_col, E_gem , 'o-', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Drug Concentration $\mathcal{D}_{\mathrm{Gemcitabine}}$ (nM)', 'Interpreter', 'latex');
ylabel('$E_{\mathrm{G2/M}}$', 'Interpreter', 'latex');
title('$E_{\mathrm{G2/M}}$ vs. $\mathcal{D}_{\mathrm{Gemcitabine}}$', 'Interpreter', 'latex');
set(gca, 'FontSize', 18);
grid on;
exportgraphics( f ,  fullfile( outputfile_dir, ['E_gem' '.png' ]) ,'Resolution' , 1200   );

EC50_G2M_d  = 10.066;
E_d_gem = E(Dose_col, EC50_G2M_d  );
% Plot
f = figure;
plot(Dose_col, E_d_gem , 'o-', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Drug Concentration $\mathcal{D}_{\mathrm{Gemcitabine}}$ (nM)', 'Interpreter', 'latex');
ylabel('$E_{d, \mathrm{G2/M}}$', 'Interpreter', 'latex');
title('$E_{d, \mathrm{G2/M}}$ vs. $\mathcal{D}_{\mathrm{Gemcitabine}}$', 'Interpreter', 'latex');
set(gca, 'FontSize', 18);
grid on;
exportgraphics( f ,  fullfile( outputfile_dir, ['E_d_gem' '.png' ]) ,'Resolution' , 1200   );


