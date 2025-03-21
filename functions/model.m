function f = model(x,P)

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
psi     = x(14);
phi     = x(15);
Bgcbbar = x(16);

pi = 0;
SDF = P.beta;
Q = 1;
pstar = 1;
Z = 1;

condexpwb  =exp(mudist+sigdist^2/2)*normcdf((-log(omegab)+mudist+sigdist^2)/sigdist)/(.5-.5*erf((log(omegab)-mudist)/(sqrt(2)*sigdist)));
condexpwl  =exp(mudist+sigdist^2/2)*normcdf((-log(omegal)+mudist+sigdist^2)/sigdist)/(.5-.5*erf((log(omegal)-mudist)/(sqrt(2)*sigdist)));
Fomegab  = normcdf((log(omegab)-mudist)/sigdist);
Fomegal  = normcdf((log(omegal)-mudist)/sigdist);
condexpwlb =(condexpwl*(1-Fomegal)-condexpwb*(1-Fomegab))/(Fomegab-Fomegal+(abs(omegab-omegal)<1e-7));


L = ( P.eta_D*D^((P.varepsilon-1)/P.varepsilon) + ...
      eta_M*M^((P.varepsilon-1)/P.varepsilon) + ...
      eta_DC*DDC^((P.varepsilon-1)/P.varepsilon) ) ^ (P.varepsilon/(P.varepsilon-1));

Rd  = (1+pi)/P.beta*(1 - C^(P.gamma)/L*vartheta*P.eta_D*(L/D)^(1/P.varepsilon));
I = K*(1-(1-P.delta)*Omega);

Ra = Rd/((1+pi)*((1-Fomegab)*omegab+Fomegal*omegal+(Fomegab-Fomegal)*condexpwlb));

Rb = omegab*(Ra*(1+pi));
Rl = omegal*(Ra*(1+pi));
Rk = Ra-(1-P.delta);
X  = (P.epsilon-1)/P.epsilon;

Delta = 1;

W     = (1-P.alpha)*X*(Rk/(P.alpha*X))^(P.alpha/(P.alpha-1));
H     = (((1-P.alpha)*Z*X)/W)^(1/P.alpha)*Omega*K;
Y     = Z/Delta*H^(1-P.alpha)*(K*Omega)^P.alpha;

Xi1    = P.epsilon/(P.epsilon-1)*X*Y/(1-P.theta*P.beta);
Xi2    = Xi1;

Bgcb = Y*P.Bbar*Bgcbbar;
Bg   = P.Bbar*Y - Bgcb;

N = (K-D/(1-psi)*(Fomegab-Fomegal))/(phi*(1- Fomegab)+(Fomegab-Fomegal)/(1-psi)) ;
Phib = (N*(phi*(1-psi)-1)-D)*(1-Fomegab);
Phil = (N+D)*Fomegal - Bg;

Gammab = Phil/((Phil^P.lambda+Phib^P.lambda)^(1/P.lambda));
Gammal = Phib/((Phil^P.lambda+Phib^P.lambda)^(1/P.lambda));
varphi    = 1/((Phib/Phil)^P.lambda+1);

Bcb     = - Phib*(1-Gammab) - Bgcb + Phil*(1-Gammal) + DDC + M;

Rdf = Rb-(1-Gammab*varphi)*P.chi;
Rlf =  Rdf + P.chi;
Rdc = P.Rdc;
Rcb = Rdf - P.chiCB;
Rg  = Rl;
Qg = P.zeta/(Rg+1-P.zeta);

tassets = Bg+K+Phil;

f(1) = L;
f(2) = Rd;
f(3) = I;
f(4) = Ra;
f(5) = Rb;
f(6) = Rl;
f(7) = Rk;
f(8) = W;
f(9) = H;
f(10) = Y;
f(11) = Bgcb;
f(12) = Bg;
f(13) = N;
f(14) = Phib;
f(15) = Phil;
f(16) = Gammab;
f(17) = Gammal;
f(18) = varphi;
f(19) = Bcb;
f(20) = Rdf;
f(21) = Rlf;
f(22) = Rdc;
f(23) = Rcb;
f(24) = Rg;
f(25) = Qg;
f(26) = tassets;
f(27) = condexpwb;
f(28) = condexpwl;
f(29) = Fomegab;
f(30) = Fomegal;
f(31) = condexpwlb;
