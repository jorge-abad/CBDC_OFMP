//CBDC MODEL
// Jorge Abad, Galo Nuño, Carlos Thomas
// October 2023

clc;

load_guess = 0;

if demand > 0.01
    load_guess = 0;
end

if load_guess
    load(strcat('simul_transition_',num2str(demand-increments,4),'_',num2str(rho_DC,'%0.4f'),'.mat'),'oo_');
    path = oo_.endo_simul;
end

load('../calibration/param_init.mat')
load('../A1_guesses.mat')

p = cell2mat(struct2cell(P));

cbdcy_grid = linspace(1.1054e-06,0.2,O.grid);

for i=1:(O.grid-length(cbdcy_grid(cbdcy_grid>demand)))

    P.cbdcy = cbdcy_grid(i);
        x0(1) = X0.omegal(i);
        x0(2) = X0.omegab(i);
        x0(3) = X0.K(i);
        x0(4) = X0.C(i);
        x0(5) = X0.D(i);
        x0(6) = X0.M(i);
        x0(7) = X0.DDC(i);
        x0(8) = X0.Omega(i);
        x0(9) = X0.eta_DC(i);

    f = @(x0) steadystate_baseline(x0,P);
    [x,R_val] = fsolve(f,x0,O.options);
    
    x(10)      = P.vartheta;
    x(11)      = P.eta_M;
    x(12)      = P.mudist;
    x(13)      = P.sigdist;
    x(14)      = P.psi;
    x(15)      = P.phi;
    x(16)      = P.Bgcbbar;
    
    f = model(x,P);

    p(2) = f(2);

    if i==1;p1=p;x1=x;f1=f;end;
end

rho_Rb = 0.8;
Rbar_init   = p1(2);
Rbar_fin    = p(2);
eta_DC_init = p(7);
eta_DC_fin  = x(9);

NN = 80;

eta_DC_path = eta_DC_init*ones(NN+2,1);
Rbar_path = Rbar_init*ones(NN+2,1);

for i=2:length(eta_DC_path)
    eta_DC_path(i) = rho_DC*eta_DC_path(i-1) + (1-rho_DC)*eta_DC_fin;
    Rbar_path(i) = rho_Rb*Rbar_path(i-1) + (1-rho_Rb)*Rbar_fin;
end

eta_DC_path=eta_DC_path(2:end-1);
Rbar_path=Rbar_path(2:end-1);

// VARIABLES
var Y Q I C K N W H SDF X pi pstar Xi1 Xi2 Delta Ra Rk Rl Rb Rdf Rlf Rg Rd Gammab Gammal Phib Phil;
var varphi omegab omegal Bgcb Bg D Omega;
var L M DDC Bcb Rcb u v g V Vt;
var Z Fomegab Fomegal condexpwb condexpwl condexpwlb Qg real_rate;

varexo Rbar eta_DC;

parameters beta gamma kappa vartheta eta_D eta_M varepsilon;
parameters alpha delta epsilon theta iota varsigma phi mudist sigdist;
parameters lambda nu rho pibar chi chiCB psi Bbar Bgcbbar;
parameters Rdc;


beta        = p(1);
gamma       = p(2);
kappa       = p(3);
vartheta    = p(4);
eta_D       = p(5);
eta_M       = p(6);
varepsilon  = p(8);
alpha       = p(9);
delta       = p(10);
epsilon     = p(11);
theta       = p(12);
iota        = p(13);
varsigma    = p(14);
phi         = p(15);
mudist      = p(16);
sigdist     = p(17);
lambda      = p(18);
nu          = p(19);
rho         = p(20);
pibar       = p(21);
chi         = p(22);
chiCB       = p(23);
psi         = p(24);
Bbar        = p(25);
Bgcbbar     = p(26);

Rdc         = p(28);

model(block,differentiate_forward_vars);

SDF*Rd/(1+pi(+1))  = 1 - C^(gamma)/L*vartheta*eta_D*(L/D)^(1/varepsilon);
SDF/(1+pi(+1))     = 1 - C^(gamma)/L*vartheta*eta_M*(L/M)^(1/varepsilon);
SDF*Rdc/(1+pi(+1)) = 1 - C^(gamma)/L*vartheta*eta_DC*(L/DDC)^(1/varepsilon);
W   = H^kappa / C^(-gamma);
SDF = beta*C(+1)^(-gamma)/C^(-gamma);
1   = Q* (1-iota*(I/I(-1)-1)^2 - iota*2*(I/I(-1)-1)*I/I(-1))+ SDF*Q(+1)*iota*2*(I(+1)/I-1)*(I(+1)/I)^2; 
K   =  (1-iota*(I/I(-1)-1)^2)*I+ (1-delta)*Omega(-1)*K(-1);
L = (eta_D*D^((varepsilon-1)/varepsilon)+eta_M*M^((varepsilon-1)/varepsilon)+eta_DC*DDC^((varepsilon-1)/varepsilon))^(varepsilon/(varepsilon-1));

Y     = Z/Delta*H^(1-alpha)*(K(-1)*Omega(-1))^alpha;
1     = theta*(1+pi)^(epsilon-1) + (1-theta)*pstar^(1-epsilon);
pstar = Xi1/Xi2;
Xi1   = epsilon/(epsilon-1)*X*Y   + theta*SDF*(1+pi(+1))^epsilon    *Xi1(+1);
Xi2   = Y                          + theta*SDF*(1+pi(+1))^(epsilon-1)*Xi2(+1);
Delta = (1-theta)*pstar^(-epsilon) + theta     *(1+pi)^epsilon       *Delta(-1);
Ra    = (Rk+(1-delta)*Q)/Q(-1);
Rk    = alpha*X*Z*(((1-alpha)*X*Z)/W)^((1-alpha)/alpha);
H     = (((1-alpha)*Z*X)/W)^(1/alpha)*Omega(-1)*K(-1);

Q*K     = N*phi*(1- Fomegab) + (N+D)/(1-psi)*(Fomegab-Fomegal);
Bcb     = psi*phi*N*(1- Fomegab) + psi/(1-psi)*(N+D)*(Fomegab-Fomegal);
N       = varsigma*( Ra*Q(-1)*Omega(-1)*K(-1) - Rb(-1)*Phib(-1)/(1+pi) - Rcb(-1)*Bcb(-1)/(1+pi) + Rl(-1)*Phil(-1)/(1+pi) + Rg*Bg(-1)/(1+pi) - Rd(-1)*D(-1)/(1+pi) );
omegab  = Rb/(Ra(+1)*(1+pi(+1)));
omegal  = Rl/(Ra(+1)*(1+pi(+1)));
1/Qg  = Rl;
Rd      = (1-Fomegab)*Rb + Fomegal*Rl + (Fomegab-Fomegal)*Ra(+1)*(1+pi(+1))*condexpwlb;

Phib = (N*(phi*(1-psi)-1)-D)*(1-Fomegab);
Phil = (N+D)*Fomegal - Bg;
Gammab = Phil/((Phil^lambda+Phib^lambda)^(1/lambda));
Gammal = Phib/((Phil^lambda+Phib^lambda)^(1/lambda));
Rb     = varphi*Gammab*Rdf + (1-varphi*Gammab)*Rlf;
Rl     = (1-varphi)*Gammal*Rlf + (1-(1-varphi)*Gammal)*Rdf;
varphi    = 1/((Phib/Phil)^lambda+1);

Rdf = rho*(Rdf(-1) + (1-varphi(-1))*chi) + (1-rho)*(Rbar+nu*(pi-pibar)) - (1-varphi)*chi;
% Rdf = rho*Rdf(-1) + (1-rho)*(Rbar+nu*(pi-pibar));

Rlf =  Rdf + chi;
Rcb = Rdf - chiCB;
Bgcb + Bcb + Phib*(1-Gammab) = Phil*(1-Gammal) + M + DDC;
Bgcb = steady_state(Y)*Bbar*Bgcbbar;

Bbar*steady_state(Y) = Bgcb + Bg;

Rg = 1/Qg(-1);

Omega = (phi*(1-Fomegab)*condexpwb + (N+D)/((1-psi)*N)*(Fomegab-Fomegal)*condexpwlb)/(phi*(1-Fomegab)+(N+D)/((1-psi)*N)*(Fomegab-Fomegal));
C = Y - I;

condexpwb  =exp(mudist+sigdist^2/2)*normcdf((-log(omegab)+mudist+sigdist^2)/sigdist)/(.5-.5*erf((log(omegab)-mudist)/(sqrt(2)*sigdist)));
condexpwl  =exp(mudist+sigdist^2/2)*normcdf((-log(omegal)+mudist+sigdist^2)/sigdist)/(.5-.5*erf((log(omegal)-mudist)/(sqrt(2)*sigdist)));
condexpwlb =(condexpwl*(1-Fomegal)-condexpwb*(1-Fomegab))/(Fomegab-Fomegal+(abs(omegab-omegal)<1e-7));
Fomegab  = normcdf((log(omegab)-mudist)/sigdist);
Fomegal  = normcdf((log(omegal)-mudist)/sigdist);

Z = 1;

real_rate = Rd*(1+pi(+1)) - 1;

u = log(C);
v = vartheta*log(L);
g = H^(1+kappa)/(1+kappa);
V = u + v + g + beta*V(+1);
Vt = u + v + g;

end;

initval;
SDF = beta;
Q = 1;
X = (epsilon-1)/epsilon;
pi = 0;
pstar = 1;
Z = 1;
Xi1    = epsilon/(epsilon-1)*X*f1(10)/(1-theta*SDF);
Xi2    = f1(10)/(1-theta*SDF);

Delta = 1;

condexpwb = f1(27);
condexpwl = f1(28);
condexpwlb = f1(31);

omegal  = x1(1);
omegab  = x1(2);
K       = x1(3);
C       = x1(4);
D       = x1(5);
M       = x1(6);
DDC     = x1(7);
Omega   = x1(8);
eta_DC  = x1(9);

L = f1(1);
Rd = f1(2);
I = f1(3);
Ra = f1(4);
Rb = f1(5);
Rl = f1(6);
Rk = f1(7);
W = f1(8);
H = f1(9);
Y = f1(10);
Bgcb = f1(11);
Bg = f1(12);
N = f1(13);
Phib = f1(14);
Phil = f1(15);
Gammab = f1(16);
Gammal = f1(17);
varphi = f1(18);
Bcb = f1(19);
Rdf = f1(20);
Rlf = f1(21);
Rcb = f1(23);
Rg = f1(24);
Fomegab = f1(29);
Fomegal = f1(30);

Qg = 1/Rg;

Rbar = Rbar_init;
eta_DC = 0.3;

real_rate = Rd-1;

u = log(C);
v = vartheta*log(L);
g = H^(1+kappa)/(1+kappa);
V = (u + v + g)/(1-beta);
Vt = (u + v + g);

end;

%steady(tolf=1e-4);

endval;
SDF = beta;
Q = 1;
X = (epsilon-1)/epsilon;
pi = 0;
pstar = 1;
Z = 1;
Xi1    = epsilon/(epsilon-1)*X*f(10)/(1-theta*SDF);
Xi2    = f(10)/(1-theta*SDF);

Delta = 1;

condexpwb = f(27);
condexpwl = f(28);
condexpwlb = f(31);

omegal  = x(1);
omegab  = x(2);
K       = x(3);
C       = x(4);
D       = x(5);
M       = x(6);
DDC     = x(7);
Omega   = x(8);
eta_DC  = x(9);

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
Rcb = f(23);
Rg = f(24);
Fomegab = f(29);
Fomegal = f(30);
Rbar = Rbar_fin;

Qg = 1/Rg;


eta_DC = eta_DC_fin;

u = log(C);
v = vartheta*log(L);
g = H^(1+kappa)/(1+kappa);
V = (u + v + g)/(1-beta);
Vt = (u + v + g);

end;

shocks;
var eta_DC;
periods 1:80;
values (eta_DC_path);

var Rbar;
periods 1:80;
% values (Rbar_path); % Gradual change in Rbar
values (Rbar_fin); % Abrupt change in Rbar
end;

%steady(tolf=1e-4);

perfect_foresight_setup(periods=400,endval_steady);

if load_guess
    oo_.endo_simul = path;
end

perfect_foresight_solver(stack_solve_algo=3);

save(strcat('simul_transition_',num2str(demand*100,'%1.2f'),'_',num2str(rho_DC,'%0.4f'),'.mat'),'Simulated_time_series','oo_','P','demand')
