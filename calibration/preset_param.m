function [P,S,O] = preset_param

% Household
P.beta      = 0.995;    % Household discount factor
P.gamma     = 1;        % CRRA instantaneous utility function parameter
P.kappa     = 0.276;    % Inverse Frisch elasticity

% HH liquidity preferences
P.vartheta   = 0.0359;  % Liquidity preference parameter
P.eta_D     = 1.00;     % Relative weight of deposits (normalized to 1)
P.eta_M     = 1.197;    % Relative weight of cash
P.eta_DC    = 0.30;     % Relative weight of CBDC
P.varepsilon= 6.6;      % Elasticity of substitution


% Production
P.alpha     = 0.33;     % Capital share
P.delta     = 0.025;    % Depreciation rate
P.epsilon   = 4.167;    % Markup
P.theta     = 0.779;    % Calvo parameter
P.iota      = 1.728;    % Investment adjustment cost

% Banks
P.varsigma  = 0.975;    % Dividend ratio
P.phi       = 15;       % Leverage constraint 
P.mudist    = -0.0014;  % Mean of idiosyncratic shocks
P.sigdist   = 0.0025;   % Std of idiosyncratic shocks
P.lambda    = 165;      % IB mkt matching function parameter

% Monetary policy
P.nu        = 1.5;      % Taylor rule inflation
P.rho       = 0.8;      % Persistence
P.pibar     = 0;        % Steady state inflation
P.chi       = 0.0025;   % Corridor width
P.chiCB     = 0.0;      % TLO spread
P.psi       = 0.00;     % TLO allowance

% Government
P.Bbar      = 0.6234*4; % Debt to SS GDP
P.Bgcbbar   = 0.16/(P.Bbar/4);% SS share of Govt bond holdings by Central Bank 
P.zeta      = 0.05; % Bond maturity

P.Rdf_target = 1.0025;      % Deposit facility rate
P.M_target = 0.1054;        % Banknotes as % of GDP
P.IB_target = 0.1884;       % Interbank assets
P.equity_target = 0.0788;   % Equity as % of TA
P.CB_target = 0.16;         % CB assets as % of GDP


% Steady state guesses (to be used as initial guess for the solver)
S.omegal = 0.9974;
S.omegab = 0.9975;
S.K      = 16.54;
S.C      = 1.5767;
S.D      = 18.72;
S.M      = 0.9842;

% Remuneration of CBDC
P.Rdc   = 1;

% Other parameters
O.options = optimoptions('fsolve','display','off','FunctionTolerance',1e-30,'Algorithm','trust-region');

O.grid = 500;