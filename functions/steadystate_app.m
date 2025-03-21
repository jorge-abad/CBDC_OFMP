function R = steadystate_app(x0,P)

omegal  = x0(1);
omegab  = x0(2);
K       = x0(3);
C       = x0(4);
D       = x0(5);
M       = x0(6);
DDC     = x0(7);
Omega   = x0(8);
eta_DC  = x0(9);
Bgcbbar = x0(10);

vartheta= P.vartheta;
eta_M   = P.eta_M;
mudist  = P.mudist;
sigdist = P.sigdist;
psi     = P.psi;
phi     = P.phi;

x0(10) = P.vartheta;
x0(11) = P.eta_M;
x0(12) = P.mudist;
x0(13) = P.sigdist;
x0(14) = P.psi;
x0(15) = P.phi;
x0(16) = Bgcbbar;

f = model(x0,P);

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


R(1) = N - (P.varsigma*(Ra*Omega*K-Rb*Phib+Rl*Phil+Rg*Bg-Rd*D-Rcb*Bcb));
R(2) = W - H^P.kappa/(C)^-P.gamma;
R(3) = Y - C - I;
R(4) = Omega - (phi*(1-Fomegab)*condexpwb + (N+D)/((1-psi)*N)*(Fomegab-Fomegal)*condexpwlb)/(phi*(1-Fomegab)+(N+D)/((1-psi)*N)*(Fomegab-Fomegal));
R(5) = (1 - vartheta*(1/L)/(C^-P.gamma)*(eta_M)*(L/(M))^(1/P.varepsilon)) - P.beta;
R(6) = (1 - vartheta*(1/L)/(C^-P.gamma)*(eta_DC)*(L/(DDC))^(1/P.varepsilon)) - P.beta*P.Rdc;
R(7) = Rl-Rlf+(1-Gammal*(1-varphi))*P.chi;
R(8) = D + N + Bcb + Phib*(1-Gammab) - K - Bg - Phil*(1-Gammal);
R(9) = DDC/Y/4 - P.cbdcy;
R(10) = Phil*(1-Gammal)/Y/4*100 - P.reserves;