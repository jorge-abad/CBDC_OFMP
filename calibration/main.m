% -------------------------------------------------------------------------
%
% This calibrates the baseline model before CBDC is introduced
%
% -------------------------------------------------------------------------

% Load preset parameters
[P,S,O] = preset_param;
% Set initial guess
x0 = [S.omegal; S.omegab; S.K; S.C; S.D; S.M; 0.000001; 1; 
    P.eta_DC; P.vartheta; P.eta_M; P.mudist; P.sigdist; P.psi; P.phi; P.Bgcbbar]; 

P.lambda = p;

% Solve steady state
f = @(x0) steadystate_calibration(x0,P);
[x,R_val] = fsolve(f,x0,O.options);

% Save steady state variables
omegal  = x(1);
omegab  = x(2);
K       = x(3);
C       = x(4);
D       = x(5);
M       = x(6);
DDC     = x(7);
Omega   = x(8);
eta_DC  = x(9);
vartheta= x(10);
eta_M   = x(11);
mudist  = x(12);
sigdist = x(13);
psi     = max(0,x(14));
phi     = x(15);
Bgcbbar = x(16);

f = model(x,P);

L = f(1);
Rd = f(2);
I = f(3);
Ra = f(4);
Rb = f(5);
Rl = f(6);
Rk = f(7);
W = f(8);
H = f(9);
Y = f(10);
Bgcb = f(11);
Bg = f(12);
N = f(13);
Phib = f(14);
Phil = f(15);
Gammab = f(16);
Gammal = f(17);
varphi = f(18);
Bcb = f(19);
Rdf = f(20);
Rlf = f(21);
Rdc = f(22);
Rcb = f(23);
Rg = f(24);
Qg = f(25);
tassets = f(26);
condexpwb = f(27);
condexpwl = f(28);
Fomegab = f(29);
Fomegal = f(30);
condexpwlb = f(31);

% -------------------------------------------------------------------------

% Save calibrated parameters
P.vartheta = vartheta;
P.eta_M    = eta_M;
P.mudist   = mudist;
P.sigdist  = sigdist;
P.psi      = psi;
P.phi      = phi;
P.Bgcbbar  = Bgcbbar;
P.Rbar     = Rdf;
P.cbdcy    = DDC/Y/4;
if neutral == 1
    P.Rdc = (Rd*D + DDC + M)/(D+DDC+M); % Neutral rate
end

if neutral == 0
    save('param_init.mat','P','x','O','-mat')
else
    save('param_init_neutral.mat','P','x','O','-mat')
end

% -------------------------------------------------------------------------