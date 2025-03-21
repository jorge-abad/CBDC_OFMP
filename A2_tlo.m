% -------------------------------------------------------------------------
%
% This runs the baseline exercise: introduction of CBDC for different
% demand levels
%
% -------------------------------------------------------------------------

% Load calibrated parameters
load('calibration/param_init'); clear x0;

load('baseline_eqm.mat')
P.reserves = reserves; clear reserves K;

% Grid of CBDC demand as a % of annualized GDP
cbdcy_grid = linspace(1.1054e-06,0.16,O.grid);

load('A1_guesses.mat')

for i=1:length(cbdcy_grid)
    
    % Set initial guess
    x0(1) = X0.omegal(i);
    x0(2) = X0.omegab(i);
    x0(3) = X0.K(i);
    x0(4) = X0.C(i);
    x0(5) = X0.D(i);
    x0(6) = X0.M(i);
    x0(7) = X0.DDC(i);
    x0(8) = X0.Omega(i);
    x0(9) = X0.eta_DC(i);
    x0(10) = P.psi;

    P.cbdcy = cbdcy_grid(i);
    
    % Solve steady state
    f = @(x0) steadystate_tlo(x0,P);
    [x,R_val] = fsolve(f,x0,O.options);
    
    omegal(i)  = x(1);
    omegab(i)  = x(2);
    K(i)       = x(3);
    C(i)       = x(4);
    D(i)       = x(5);
    M(i)       = x(6);
    DDC(i)     = x(7);
    Omega(i)   = x(8);
    eta_DC(i)  = x(9);
    x(14)      = x(10);
    x(10)      = P.vartheta;
    x(11)      = P.eta_M;
    x(12)      = P.mudist;
    x(13)      = P.sigdist;
    x(15)      = P.phi;
    x(16)      = P.Bgcbbar;

    f = model(x,P);

    L(i) = f(1);
    Rd(i) = f(2);
    I(i) = f(3);
    Ra(i) = f(4);
    Rb(i) = f(5);
    Rl(i) = f(6);
    Rk(i) = f(7);
    W(i) = f(8);
    H(i) = f(9);
    Y(i) = f(10);
    Bgcb(i) = f(11);
    Bg(i) = f(12);
    N(i) = f(13);
    Phib(i) = f(14);
    Phil(i) = f(15);
    Gammab(i) = f(16);
    Gammal(i) = f(17);
    varphi(i) = f(18);
    Bcb(i) = f(19);
    Rdf(i) = f(20);
    Rlf(i) = f(21);
    Rdc(i) = f(22);
    Rcb(i) = f(23);
    Rg(i) = f(24);
    Qg(i) = f(25);
    tassets(i) = f(26);
    Fomegab(i) = f(29);
    Fomegal(i) = f(30);
    Rib(i) = Rdf(i) + (1-varphi(i))*P.chi;

end

P.options = optimoptions('fsolve','display','off','FinDiffType','central','FunctionTolerance',1e-30,'Algorithm','trust-region-dogleg');

for i=1:length(cbdcy_grid)
    
    x0 = [omegal(i);omegab(i);K(i);C(i);D(i);M(i);DDC(i);Omega(i);eta_DC(i);x(14)];
    P.cbdcy = cbdcy_grid(i);
    
    % Solve steady state
    f = @(x0) steadystate_tlo(x0,P);
    [x,R_val] = fsolve(f,x0,P.options);    
    
    omegal(i)  = x(1);
    omegab(i)  = x(2);
    K(i)       = x(3);
    C(i)       = x(4);
    D(i)       = x(5);
    M(i)       = x(6);
    DDC(i)     = x(7);
    Omega(i)   = x(8);
    eta_DC(i)  = x(9);
    x(14)      = x(10);
    x(10)      = P.vartheta;
    x(11)      = P.eta_M;
    x(12)      = P.mudist;
    x(13)      = P.sigdist;
    x(15)      = P.phi;
    x(16)      = P.Bgcbbar;

    f = model(x,P);

    L(i) = f(1);
    Rd(i) = f(2);
    I(i) = f(3);
    Ra(i) = f(4);
    Rb(i) = f(5);
    Rl(i) = f(6);
    Rk(i) = f(7);
    W(i) = f(8);
    H(i) = f(9);
    Y(i) = f(10);
    Bgcb(i) = f(11);
    Bg(i) = f(12);
    N(i) = f(13);
    Phib(i) = f(14);
    Phil(i) = f(15);
    Gammab(i) = f(16);
    Gammal(i) = f(17);
    varphi(i) = f(18);
    Bcb(i) = f(19);
    Rdf(i) = f(20);
    Rlf(i) = f(21);
    Rdc(i) = f(22);
    Rcb(i) = f(23);
    Rg(i) = f(24);
    Qg(i) = f(25);
    tassets(i) = f(26);
    Fomegab(i) = f(29);
    Fomegal(i) = f(30);
    Rib(i) = Rdf(i) + (1-varphi(i))*P.chi;

end

set(0,'DefaultLineLineWidth',1.5)
set(0,'DefaultLineMarkerSize',5)
set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize',20)
set(0,'DefaultFigureUnits','pixels')
set(0,'DefaultFigureRenderer','painters')

fig = 1;

if fig
    set(0,'DefaultAxesFontSize',18)
    set(0,'DefaultFigurePosition',[100 50 1200 400])

    figure(1);
    subplot(1,3,1); grid on; hold on; box on;
    plot(DDC./Y*25, Bcb./Y./4*100,'--','LineWidth',3);
    title('(a) Central bank balance sheet','interpreter','latex','FontSize',22)
    xline(6.3,'--r','€3,000 p.a.','LabelVerticalAlignment','middle')
    legend({'Govt. bond holdings', 'Targeted loans'},'interpreter','latex','Location','southeast','FontSize',16)
    xlabel('CBDC in circulation (\% of GDP)','interpreter','latex','FontSize',18)
    ylabel('\% of GDP','interpreter','latex','FontSize',18)
    xlim([DDC(1)./Y(1)*25 14])
        
    subplot(1,3,2); grid on; hold on; box on; 
    plot(DDC./Y*25, Bcb./tassets*100,'--','LineWidth',3);
    title('(b) Banks'' balance sheet','interpreter','latex','FontSize',22)
     xline(6.3,'--r','€3,000 p.a.','LabelVerticalAlignment','middle')
   legend({'Govt. bond holdings', 'Targeted loans'},'interpreter','latex','Location','southeast','FontSize',16)
    xlabel('CBDC in circulation (\% of GDP)','interpreter','latex','FontSize',18)
    ylabel('\% of total bank assets','interpreter','latex','FontSize',18)
    xlim([DDC(1)./Y(1)*25 14])
    
    subplot(1,3,3); grid on; hold on; box on; 
    plot(DDC./Y*25, K./K(1).*100,'-','LineWidth',3);
    plot(DDC./Y*25, N./N(1).*100,'--','LineWidth',3);
    plot(DDC./Y*25, Y./Y(1).*100,':','LineWidth',3);
    title('(c) Bank credit, equity and output','interpreter','latex','FontSize',22)
    xline(6.3,'--r','€3,000 p.a.','LabelVerticalAlignment','middle')
    legend({'Bank credit ($K$)', 'Bank equity ($N$)', 'Output ($Y$)'},'interpreter','latex','Location','southeast','FontSize',16)
    xlabel('CBDC in circulation (\% of GDP)','interpreter','latex','FontSize',18)
    ylabel('\% of baseline','interpreter','latex','FontSize',18)
    xlim([DDC(1)./Y(1)*25 14])

    saveas(gcf,'graphs\fig_policies_slides','epsc');
end
