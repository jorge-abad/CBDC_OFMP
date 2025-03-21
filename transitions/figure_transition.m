set(0,'DefaultFigurePosition',[100 50 1100 900])
set(0,'DefaultLineLineWidth',2)
set(0,'DefaultLineMarkerSize',5)
set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize',12)
set(0,'DefaultFigureUnits','pixels')
set(0,'DefaultFigureRenderer','painters')


if j == 1
    line = '-';
    line2 = '-o';
elseif j==2
    line = '--';
    line2 = '--o';
else
    line = ':';
    line2 = ':o';
end

load(strcat('simul_transition_',num2str(demand*100,'%0.2f'),'_',num2str(rho_DC,'%0.4f'),'.mat'))

T = 60;

t = 1:1:length(Simulated_time_series.data(1:T,1));

V = struct();

for i = 1:length(Simulated_time_series.name)
    V.(Simulated_time_series.name{i}) = i;
end

figure(1);
subplot(3,3,1); grid on; hold on; box on;
p=plot(t, 100*(Simulated_time_series.data(t,V.DDC)./Simulated_time_series.data(t,V.Y)/4-Simulated_time_series.data(1,V.DDC)./Simulated_time_series.data(1,V.Y)/4),'LineStyle',line);
title('(a) CBDC in circulation  (\% of GDP)','interpreter','latex','fontsize',16)
xlabel('$t$','interpreter','latex','FontSize',14)
ylabel('pp diff.','interpreter','latex','FontSize',14)
xlim([t(1) t(end)])
c = get(p,'Color');

subplot(3,3,2); grid on; hold on; box on;
plot(t, 100*(Simulated_time_series.data(t,V.M)./Simulated_time_series.data(t,V.Y)/4-Simulated_time_series.data(1,V.M)./Simulated_time_series.data(1,V.Y)/4),'LineStyle',line);
title('(b) Banknotes in circulation (\% of GDP)','interpreter','latex','FontSize',16)
xlabel('$t$','interpreter','latex','FontSize',14)
ylabel('pp diff.','interpreter','latex','FontSize',14)
xlim([t(1) t(end)])

subplot(3,3,3); grid on; hold on; box on;
plot(t, 100*(Simulated_time_series.data(t,V.D)./Simulated_time_series.data(t,V.Y)/4-Simulated_time_series.data(1,V.D)./Simulated_time_series.data(1,V.Y)/4),'LineStyle',line);
title('(c) Retail deposits (\% of GDP)','interpreter','latex','FontSize',16)
xlabel('$t$','interpreter','latex','FontSize',14)
ylabel('pp diff.','interpreter','latex','FontSize',14)
xlim([t(1) t(end)])
subplot(3,3,4); grid on; hold on; box on;
plot(t, 100*(Simulated_time_series.data(t,V.Phil).*(1-Simulated_time_series.data(t,V.Gammal))./Simulated_time_series.data(t,V.Y)/4-Simulated_time_series.data(1,V.Phil).*(1-Simulated_time_series.data(1,V.Gammal))./Simulated_time_series.data(1,V.Y)/4),'LineStyle',line);
title('(d) Excess reserves (\% of GDP)','interpreter','latex','FontSize',16)
xlabel('$t$','interpreter','latex','FontSize',14)
ylabel('pp diff.','interpreter','latex','FontSize',14)
xlim([t(1) t(end)])

subplot(3,3,5); grid on; hold on; box on;
plot(t, 400*Simulated_time_series.data(t,V.Rd)-400,line2,'MarkerIndices',[1:5:60],'color',c);
plot(t, 400*Simulated_time_series.data(t,V.Rdf)-400,'LineStyle',line,'color',c);
title('(e) Policy and interbank rates','interpreter','latex','FontSize',16)
legend({'$R^{IB}$', '$R^{DF}$'},'interpreter','latex','Location','east','FontSize',14)
xlabel('$t$','interpreter','latex','FontSize',14)
ylabel('Annualized \%','interpreter','latex','FontSize',14)
xlim([t(1) t(end)])

subplot(3,3,6); grid on; hold on; box on;
plot(t, 400*Simulated_time_series.data(t,V.pi),'LineStyle',line);
title('(f) Inflation (ann. \%)','interpreter','latex','FontSize',16)
xlabel('$t$','interpreter','latex','FontSize',14)
ylabel('pp diff.','interpreter','latex','FontSize',14)
xlim([t(1) t(end)])
if figure_index == 1
    legend({'Gradual transition', 'Fast transition'},'interpreter','latex','Location','east')
else
    legend({'Low demand', 'Interm. demand', 'High demand'},'interpreter','latex','Location','east','FontSize',14)
end

subplot(3,3,7); grid on; hold on; box on;
plot(t, 400*Simulated_time_series.data(t,V.real_rate)-400*Simulated_time_series.data(1,V.real_rate),'LineStyle',line);
title('(g) Real rate (ann. \%)','interpreter','latex','FontSize',16)
xlabel('$t$','interpreter','latex','FontSize',14)
ylabel('pp diff.','interpreter','latex','FontSize',14)
xlim([t(1) t(end)])

subplot(3,3,8); grid on; hold on; box on;
plot(t, Simulated_time_series.data(t,V.Y)./Simulated_time_series.data(1,V.Y)*100-100,'LineStyle',line);
title('(h) Output','interpreter','latex','FontSize',16)
xlabel('$t$','interpreter','latex','FontSize',14)
ylabel('\% diff.','interpreter','latex','FontSize',14)
xlim([t(1) t(end)])

subplot(3,3,9); grid on; hold on; box on;
plot(t, Simulated_time_series.data(t,V.C)./Simulated_time_series.data(1,V.C)*100-100,'LineStyle',line);
title('(i) Consumption ','interpreter','latex','FontSize',16)
xlabel('$t$','interpreter','latex','FontSize',14)
ylabel('\% diff.','interpreter','latex','FontSize',14)
xlim([t(1) t(end)])
