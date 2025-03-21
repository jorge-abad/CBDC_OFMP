% -------------------------------------------------------------------------
%
% Abad, Nuño and Thomas (2025): "CBDC and the operational framework of 
% monetary policy"
%
% This script is the main file that solves and calibrates the baseline
% model in Abad, Nuño and Thomas (2025). It should be run first.
%
% This version: Nov 2024. Contact: jorgeabad@bde.es
%
% -------------------------------------------------------------------------

% Housekeeping
close all; clear; clc; tic;

% Add Dynare path (check that the version C:\dynare\x.y\matlab is correct)
addpath /Applications/Dynare/6.1-arm64/matlab

% Add paths
addpath functions

% -------------------------------------------------------------------------

% Run baseline exercise (fig. 2)
run('A1_baseline.m');

% Store value of excess reserves in the baseline without CBDC
reserves = Phil(1).*(1-Gammal(1))./Y(1)./4*100;
save('baseline_eqm.mat','reserves'); clear reserves;

% Calculate welfare
L = ( P.eta_D.*D.^((P.varepsilon-1)/P.varepsilon) + ...
    P.eta_M.*M.^((P.varepsilon-1)/P.varepsilon) + ...
    eta_DC.*DDC.^((P.varepsilon-1)/P.varepsilon) ) .^ (P.varepsilon/(P.varepsilon-1));
W = (log(C) + P.vartheta*log(L) + H.^(1+P.kappa)/(1+P.kappa))./(1-P.beta);
W_C = log(C)./(1-P.beta);
W_L = P.vartheta*log(L)./(1-P.beta);
W_H = -(H.^(1+P.kappa)/(1+P.kappa))./(1-P.beta);
Cbar = exp(W*(1-P.beta));
Cbar_C = exp(W_C*(1-P.beta));
Cbar_L = exp(W_L*(1-P.beta));
Cbar_H = exp(W_H*(1-P.beta));
save('welfare_baseline.mat','Cbar','Cbar_C','Cbar_L','Cbar_H','DDC','Y','eta_DC')
Cbar_baseline = Cbar(1);
clear Cbar welfare;

% -------------------------------------------------------------------------

% Run the Dynare scripts (transitional dynamics)
cd transitions;

figure_index = 1;

for demand = [0.04,0.08,0.12]
    for rho_DC = [0.80,0.0001]
        dynare transition.mod;  
    end
    j=0;
    for rho_DC = [0.80,0.0001]
        j=j+1;
        figure_transition; 
    end
    saveas(gcf,strcat('../graphs/fig_transition_',num2str(demand*100,'%0.2f'),'.eps'),'epsc');
    close;
end

figure_index = 2;
rho_DC = 0.80; j = 0;

for demand = [0.04,0.08,0.12]
    j=j+1;
    figure_transition;
end
saveas(gcf,'../graphs/fig_transition.eps','epsc');

close;
cd ..

% -------------------------------------------------------------------------

% Run the policy exercises with compensating policies (keeping agg.
% reserves constant)
run('A2_app.m');
run('A2_tlo.m');
close;

% Calculate welfare
L = ( P.eta_D.*D.^((P.varepsilon-1)/P.varepsilon) + ...
    P.eta_M.*M.^((P.varepsilon-1)/P.varepsilon) + ...
    eta_DC.*DDC.^((P.varepsilon-1)/P.varepsilon) ) .^ (P.varepsilon/(P.varepsilon-1));
W = (log(C) + P.vartheta*log(L) + H.^(1+P.kappa)/(1+P.kappa))./(1-P.beta);
Cbar = exp(W*(1-P.beta));
save('welfare_floor.mat','Cbar','DDC','Y','eta_DC')

% -------------------------------------------------------------------------

clear Rib;

% Neutral rate exercise
neutral = 1; p = P.lambda;
run('calibration/main.m');
run('A3_neutral.m');

% Calculate welfare
L = ( P.eta_D.*D.^((P.varepsilon-1)/P.varepsilon) + ...
    P.eta_M.*M.^((P.varepsilon-1)/P.varepsilon) + ...
    eta_DC.*DDC.^((P.varepsilon-1)/P.varepsilon) ) .^ (P.varepsilon/(P.varepsilon-1));
W = (log(C) + P.vartheta*log(L) + H.^(1+P.kappa)/(1+P.kappa))./(1-P.beta);
Cbar = exp(W*(1-P.beta));
save('welfare_neutral.mat','Cbar','DDC','Y','eta_DC')
reserves = Phil(1).*(1-Gammal(1))./Y(1)./4*100;
K = K(1);
save('neutral_eqm.mat','reserves'); clear reserves;

% Neutral rate + floor
run('A3_neutral_floor.m');
% Calculate welfare
L = ( P.eta_D.*D.^((P.varepsilon-1)/P.varepsilon) + ...
    P.eta_M.*M.^((P.varepsilon-1)/P.varepsilon) + ...
    eta_DC.*DDC.^((P.varepsilon-1)/P.varepsilon) ) .^ (P.varepsilon/(P.varepsilon-1));
W = (log(C) + P.vartheta*log(L) + H.^(1+P.kappa)/(1+P.kappa))./(1-P.beta);
W_C = log(C)./(1-P.beta);
W_L = P.vartheta*log(L)./(1-P.beta);
W_H = -(H.^(1+P.kappa)/(1+P.kappa))./(1-P.beta);
Cbar = exp(W*(1-P.beta));
Cbar_C = exp(W_C*(1-P.beta));
Cbar_L = exp(W_L*(1-P.beta));
Cbar_H = exp(W_H*(1-P.beta));
save('welfare_neutral_floor.mat','Cbar','Cbar_C','Cbar_L','Cbar_H','DDC','Y','eta_DC')

% -------------------------------------------------------------------------

close all; set(0,'DefaultAxesFontSize',18)

set(0,'DefaultFigurePosition',[100 50 900 600])

load 'welfare_baseline.mat'
plot(DDC./Y*25,(Cbar./Cbar(1)-1)*100,'LineWidth',3,'LineStyle','-'); hold on;
xlim([DDC(1)./Y(1)*25 12.5])

load 'welfare_floor.mat'
plot(DDC./Y*25,(Cbar./Cbar(1)-1)*100,'LineWidth',3,'LineStyle','--'); hold on;
xlim([DDC(1)./Y(1)*25 12.5])

load 'welfare_neutral.mat'
plot(DDC./Y*25,(Cbar./Cbar(1)-1)*100,'LineWidth',3,'LineStyle',':'); hold on; grid on;
xlim([DDC(1)./Y(1)*25 12.5])

load 'welfare_neutral_floor.mat'
plot(DDC./Y*25,(Cbar./Cbar(1)-1)*100,'LineWidth',3,'LineStyle','-.'); hold on; grid on;
xlim([DDC(1)./Y(1)*25 12.5])

ylabel('Equivalent consumption units ($\Delta$\% wrt baseline without CBDC)','interpreter','latex','FontSize',18)
xlabel('CBDC in circulation (\% of GDP)','interpreter','latex','FontSize',18)
legend({'Baseline scenario','Floor preserving policies','Wealth-neutral remuneration','Floor preserving policies and wealth-neutral remuneration'},...
    'Interpreter','Latex','FontSize',18,'Location','NorthWest')

saveas(gcf,'graphs\welfare','epsc');

% -------------------------------------------------------------------------

close all; set(0,'DefaultAxesFontSize',18)

set(0,'DefaultFigurePosition',[100 50 900 600])

load 'welfare_baseline.mat'

plot(DDC./Y*25,(Cbar./Cbar(1)-1)*100,'LineWidth',3,'LineStyle','-'); hold on; grid on;
xlim([DDC(1)./Y(1)*25 12.5])

plot(DDC./Y*25,(Cbar_C./Cbar_C(1)-1)*100,'LineWidth',3,'LineStyle','--'); hold on;
xlim([DDC(1)./Y(1)*25 12.5])

plot(DDC./Y*25,(Cbar_H./Cbar_H(1)-1)*100,'LineWidth',3,'LineStyle',':'); hold on;
xlim([DDC(1)./Y(1)*25 12.5])

plot(DDC./Y*25,(Cbar_L./Cbar_L(1)-1)*100,'LineWidth',3,'LineStyle','-.'); hold on;
xlim([DDC(1)./Y(1)*25 12.5])

ylabel('Equivalent consumption units ($\Delta$\% wrt baseline without CBDC)','interpreter','latex','FontSize',18)
xlabel('CBDC in circulation (\% of GDP)','interpreter','latex','FontSize',18)
legend({'Total welfare effect','Utility of consumption','Disutility of labor','Utility of liquidity services'},...
    'Interpreter','Latex','FontSize',18,'Location','NorthWest')

saveas(gcf,'graphs\fig_app_welfare_decomp','epsc');

% -------------------------------------------------------------------------