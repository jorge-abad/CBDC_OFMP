% -------------------------------------------------------------------------
%
% This runs the baseline exercise: introduction of CBDC for different
% demand levels
%
% -------------------------------------------------------------------------

% Load calibrated parameters
load('calibration/param_init'); clear x0;

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

    P.cbdcy = cbdcy_grid(i);

    % Solve steady state
    f = @(x0) steadystate_baseline(x0,P);
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
    x(10)      = P.vartheta;
    x(11)      = P.eta_M;
    x(12)      = P.mudist;
    x(13)      = P.sigdist;
    x(14)      = P.psi;
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

    x0 = [omegal(i);omegab(i);K(i);C(i);D(i);M(i);DDC(i);Omega(i);eta_DC(i)];
    P.cbdcy = cbdcy_grid(i);

    % Solve steady state
    f = @(x0) steadystate_baseline(x0,P);
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
    x(10)      = P.vartheta;
    x(11)      = P.eta_M;
    x(12)      = P.mudist;
    x(13)      = P.sigdist;
    x(14)      = P.psi;
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

X0.omegal = omegal;
X0.omegab = omegab;
X0.K = K;
X0.C = C;
X0.D = D;
X0.M = M;
X0.DDC = DDC;
X0.Omega = Omega;
X0.eta_DC = eta_DC;

save('A1_guesses.mat','X0')

set(0,'DefaultFigurePosition',[100 50 1800 1000])
set(0,'DefaultLineLineWidth',1.5)
set(0,'DefaultLineMarkerSize',5)
set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize',20)
set(0,'DefaultFigureUnits','pixels')
set(0,'DefaultFigureRenderer','painters')

fig = 1; alt_fig=0;

if alt_fig

    set(0,'DefaultAxesFontSize',12)
    set(0,'DefaultFigurePosition',[100 50 900 1000])

    figure(1);
    subplot(3,2,1); grid on; hold on; box on;
    area(DDC./Y*25, [K'-K(1), Bg'-Bg(1), ((1-Gammal).*Phil)'-((1-Gammal(1)).*Phil(1))]./(K(1)+Bg(1)+((1-Gammal(1)).*Phil(1))));
    plot(DDC./Y*25, [K'-K(1)+Bg'-Bg(1)+((1-Gammal).*Phil)'-((1-Gammal(1)).*Phil(1))]./(K(1)+Bg(1)+((1-Gammal(1)).*Phil(1))),'LineWidth',3,'Color',[1 0 0]);
    legend({'Loans ($Q K$)', 'Govt. bonds ($B^G$)', 'CB reserves ($(1-\Gamma^L)\Phi^L$)'},'interpreter','latex','Location','southwest','FontSize',10)
    title('(a) Aggregate banks'' assets','interpreter','latex','FontSize',14)
    xlabel('CBDC in circulation (\% of GDP)','interpreter','latex','FontSize',12)
    xlim([DDC(1)./Y(1)*25 13.5])
    ylim([-0.06, 0.02])

    subplot(3,2,2); grid on; hold on; box on;
    area(DDC./Y*25, [D'-D(1), N'-N(1)]./(N(1)+D(1)+(1-Gammab(1)).*Phib(1)));
    area(DDC./Y*25, ((1-Gammab').*Phib'-((1-Gammab(1)).*Phib(1)))./(N(1)+D(1)+(1-Gammab(1)).*Phib(1)));
    plot(DDC./Y*25, [D'-D(1)+N'-N(1)]./(N(1)+D(1)+((1-Gammab(1)).*Phib(1)))+((1-Gammab').*Phib'-((1-Gammab(1)).*Phib(1)))./(N(1)+D(1)+(1-Gammab(1)).*Phib(1)),'LineWidth',3,'Color',[1 0 0]);
    legend({'Deposits ($D$)', 'Equity ($N$)', 'CB borrowing ($(1-\Gamma^B)\Phi^B$)'},'interpreter','latex','Location','southwest','FontSize',10)
    title('(b) Aggregate banks'' liabilities','interpreter','latex','FontSize',14)
    xlabel('CBDC in circulation (\% of GDP)','interpreter','latex','FontSize',12)
    xlim([DDC(1)./Y(1)*25 13.5])


    % -------------------------------------------------------------------------

    % Borrowing bank
    subplot(3,2,3); grid on; hold on; box on;
    area(DDC./Y*25, ((1-Fomegab').*P.phi.*N'-(1-Fomegab(1)).*P.phi.*N(1))./((1-Fomegab(1)).*P.phi.*N(1)));
    plot(DDC./Y*25, ((1-Fomegab').*P.phi.*N'-(1-Fomegab(1)).*P.phi.*N(1))./((1-Fomegab(1)).*P.phi.*N(1)),'LineWidth',3,'Color',[1 0 0]);
    legend({'Loans ($(1-F(\omega^B))\phi N$)'},'interpreter','latex','Location','southwest','FontSize',10)
    title('(c) Borrowing banks'' assets','interpreter','latex','FontSize',14)
    xlabel('CBDC in circulation (\% of GDP)','interpreter','latex','FontSize',12)
    xlim([DDC(1)./Y(1)*25 13.5])
    ylim([-0.04, 0.04])

    subplot(3,2,4); grid on; hold on; box on;
    area(DDC./Y*25, [(1-Fomegab').*D'-(1-Fomegab(1)).*D(1), (1-Fomegab').*N'-(1-Fomegab(1)).*N(1)]./((1-Fomegab(1))*N(1)+(1-Fomegab(1))*D(1)+(Phib(1))));
    area(DDC./Y*25, [((1-Gammab).*Phib)'-((1-Gammab(1)).*Phib(1)),((Gammab).*Phib)'-((Gammab(1)).*Phib(1))]./((1-Fomegab(1))*N(1)+(1-Fomegab(1))*D(1)+(Phib(1))));
    plot(DDC./Y*25, ((1-Fomegab').*N'-(1-Fomegab(1)).*N(1) + (1-Fomegab').*D'-(1-Fomegab(1)).*D(1) + ((Gammab).*Phib)'-((Gammab(1)).*Phib(1)) + ((1-Gammab).*Phib)'-((1-Gammab(1)).*Phib(1)))./((1-Fomegab(1))*N(1)+(1-Fomegab(1))*D(1)+(Phib(1))),'LineWidth',3,'Color',[1 0 0]);
    legend({'Deposits ($(1-F(\omega^B))D$)', 'Equity ($(1-F(\omega^B))N$)', 'CB borrowing ($(1-\Gamma^B)\Phi^B$)', 'IB borrowing ($\Gamma^B \Phi^B$)'},'interpreter','latex','Location','southwest','FontSize',10)
    title('(d) Borrowing banks'' liabilities','interpreter','latex','FontSize',14)
    xlabel('CBDC in circulation (\% of GDP)','interpreter','latex','FontSize',12)
    xlim([DDC(1)./Y(1)*25 13.5])

    % -------------------------------------------------------------------------

    % Lending bank
    subplot(3,2,5); grid on; hold on; box on;
    area(DDC./Y*25, [((Gammal).*Phil)'-((Gammal(1)).*Phil(1))]./(Phil(1)+Bg(1)));
    area(DDC./Y*25, [Bg'-Bg(1),((1-Gammal).*Phil)'-((1-Gammal(1)).*Phil(1))]./(Phil(1)+Bg(1)));
    plot(DDC./Y*25, [((Gammal).*Phil)'-((Gammal(1)).*Phil(1)) + ((1-Gammal).*Phil)'-((1-Gammal(1)).*Phil(1)) + Bg'-Bg(1)]./(Phil(1)+Bg(1)),'LineWidth',3,'Color',[1 0 0]);
    legend({'IB lending ($\Gamma^L \Phi^L$)','Govt. bonds ($B^G$)','CB reserves $((1-\Gamma^L)\Phi^L)$'},'interpreter','latex','Location','southwest','FontSize',10)
    title('(e) Lending banks'' assets','interpreter','latex','FontSize',14)
    xlabel('CBDC in circulation (\% of GDP)','interpreter','latex','FontSize',12)
    xlim([DDC(1)./Y(1)*25 13.5])
    ylim([-0.05, 0.02])

    subplot(3,2,6); grid on; hold on; box on;
    area(DDC./Y*25, [(Fomegal').*D'-(Fomegal(1)).*D(1)]./((Fomegal(1))*N(1)+(Fomegal(1))*D(1)));
    area(DDC./Y*25, [(Fomegal').*N'-(Fomegal(1)).*N(1)]./((Fomegal(1))*N(1)+(Fomegal(1))*D(1)));
    plot(DDC./Y*25, [(Fomegal').*N'-(Fomegal(1)).*N(1)+(Fomegal').*D'-(Fomegal(1)).*D(1)]./((Fomegal(1))*N(1)+(Fomegal(1))*D(1)),'LineWidth',3,'Color',[1 0 0]);
    legend({'Deposits ($F(\omega^L)D$)', 'Equity ($F(\omega^L)N$)'},'interpreter','latex','Location','southwest','FontSize',10)
    title('(f) Lending banks'' liabilities','interpreter','latex','FontSize',14)
    xlabel('CBDC in circulation (\% of GDP)','interpreter','latex','FontSize',12)
    xlim([DDC(1)./Y(1)*25 13.5])
    ylim([-0.05, 0.02])

    saveas(gcf,'graphs\baseline_figure_bs','epsc');

    close;

end

if fig

    set(0,'DefaultAxesFontSize',14)
    set(0,'DefaultFigurePosition',[100 50 1800 1000])

    figure(1);
    subplot(2,3,1); grid on; hold on; box on;
    yyaxis left
    plot(DDC./Y*25, 100*M./(4*Y),'-','LineWidth',3);
    ylabel('\% of GDP','interpreter','latex','FontSize',18)
    yyaxis right
    plot(DDC./Y*25, D./(4*Y)*100,'--','LineWidth',3);
    legend({'Cash ($M$)','Deposits ($D$, rhs)'},'interpreter','latex','Location','northeast','FontSize',12)
    title('(a) Cash and retail deposits','interpreter','latex','FontSize',22)
    xlabel('CBDC in circulation (\% of GDP)','interpreter','latex','FontSize',18)
    xlim([DDC(1)./Y(1)*25 14])
    ylim([222 234])

    subplot(2,3,2); grid on; hold on; box on;
    plot(DDC./Y*25, Phil.*(1-Gammal)./Y*25,'-','LineWidth',3);
    plot(DDC./Y*25, Phib.*(1-Gammab)./Y*25,'--','LineWidth',3);
    legend({'Deposit facility','Lending facility'},'interpreter','latex','Location','northwest','FontSize',12)
    title('(b) Central bank facilities','interpreter','latex','FontSize',22)
    xlabel('CBDC in circulation (\% of GDP)','interpreter','latex','FontSize',18)
    ylabel('\% of GDP','interpreter','latex','FontSize',18)
    xlim([DDC(1)./Y(1)*25 14])

    subplot(2,3,3); grid on; hold on; box on;
    plot(DDC./Y*25, Rlf*400-400,'-','LineWidth',3);
    plot(DDC./Y*25, Rdf*400-400,'--','LineWidth',3);
    plot(DDC./Y*25, Rib*400-400,':','LineWidth',3);
    title('(c) Policy and interbank rates','interpreter','latex','FontSize',22)
    legend({'$R^{LF}$', '$R^{DF}$', '$R^{IB}$'},'interpreter','latex','Location','northeast','FontSize',12)
    xlabel('CBDC in circulation (\% of GDP)','interpreter','latex','FontSize',18)
    ylabel('Annualized \%','interpreter','latex','FontSize',18)
    xlim([DDC(1)./Y(1)*25 14])

    subplot(2,3,4); grid on; hold on; box on;
    plot(DDC./Y*25, Rd*400-400,'LineWidth',3);
    title('(d) Retail deposit rate','interpreter','latex','FontSize',22)
    xlabel('CBDC in circulation (\% of GDP)','interpreter','latex','FontSize',18)
    ylabel('Annualized \%','interpreter','latex','FontSize',18)
    xlim([DDC(1)./Y(1)*25 14])

    subplot(2,3,5); grid on; hold on; box on;
    plot(DDC./Y*25, (D+DDC+M)./(D(1)+DDC(1)+M(1))*100,'LineWidth',3);
    title('(e) Households'' wealth','interpreter','latex','FontSize',22)
    xlabel('CBDC in circulation (\% of GDP)','interpreter','latex','FontSize',18)
    ylabel('\% of baseline','interpreter','latex','FontSize',18)
    xlim([DDC(1)./Y(1)*25 14])

    subplot(2,3,6); grid on; hold on; box on;
    plot(DDC./Y*25, K./K(1).*100,'-','LineWidth',3);
    plot(DDC./Y*25, N./N(1).*100,'--','LineWidth',3);
    plot(DDC./Y*25, Y./Y(1).*100,':','LineWidth',3);
    title('(f) Bank credit, equity and output','interpreter','latex','FontSize',22)
    legend({'Bank credit ($K$)', 'Bank equity ($N$)', 'Output ($Y$)'},'interpreter','latex','Location','southwest','FontSize',12)
    xlabel('CBDC in circulation (\% of GDP)','interpreter','latex','FontSize',18)
    ylabel('\% of baseline','interpreter','latex','FontSize',18)
    xlim([DDC(1)./Y(1)*25 14])

    saveas(gcf,'graphs\baseline_figure','epsc');
end

close;

% -------------------------------------------------------------------------