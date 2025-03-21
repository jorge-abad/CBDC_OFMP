function T = static_g1_tt(T, y, x, params)
% function T = static_g1_tt(T, y, x, params)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T         [#temp variables by 1]  double   vector of temporary terms to be filled by function
%   y         [M_.endo_nbr by 1]      double   vector of endogenous variables in declaration order
%   x         [M_.exo_nbr by 1]       double   vector of exogenous variables in declaration order
%   params    [M_.param_nbr by 1]     double   vector of parameter values in declaration order
%
% Output:
%   T         [#temp variables by 1]  double   vector of temporary terms
%

assert(length(T) >= 61);

T = transition.static_resid_tt(T, y, x, params);

T(40) = params(4)*getPowerDeriv(y(4),params(2),1)/y(35);
T(41) = getPowerDeriv(y(4),(-params(2)),1);
T(42) = params(1)*getPowerDeriv(y(4)+y(54),(-params(2)),1);
T(43) = getPowerDeriv(y(5)*y(34),params(8),1);
T(44) = (y(46)-y(47))*(y(6)*(1-params(23))-(y(33)+y(6))*(1-params(23)))/(y(6)*(1-params(23))*y(6)*(1-params(23)));
T(45) = getPowerDeriv(T(20),(1-params(8))/params(8),1);
T(46) = getPowerDeriv(T(20),1/params(8),1);
T(47) = (-((y(13)+y(57))*y(9)*params(11)*getPowerDeriv(1+y(11)+y(53),params(10),1)));
T(48) = (-((y(14)+y(58))*y(9)*params(11)*getPowerDeriv(1+y(11)+y(53),params(10)-1,1)));
T(49) = (-((-(y(19)*(y(16)+y(59))))/((1+y(11)+y(53))*(y(16)+y(59))*(1+y(11)+y(53))*(y(16)+y(59)))));
T(50) = (-((-(y(18)*(y(16)+y(59))))/((1+y(11)+y(53))*(y(16)+y(59))*(1+y(11)+y(53))*(y(16)+y(59)))));
T(51) = (-((-((1+y(11)+y(53))*y(19)))/((1+y(11)+y(53))*(y(16)+y(59))*(1+y(11)+y(53))*(y(16)+y(59)))));
T(52) = (-((-((1+y(11)+y(53))*y(18)))/((1+y(11)+y(53))*(y(16)+y(59))*(1+y(11)+y(53))*(y(16)+y(59)))));
T(53) = getPowerDeriv(T(25),1/params(17),1);
T(54) = getPowerDeriv(y(26),params(17),1)*T(53);
T(55) = getPowerDeriv(y(26)/y(27),params(17),1);
T(56) = T(53)*getPowerDeriv(y(27),params(17),1);
T(57) = getPowerDeriv(y(35)/y(33),T(3),1);
T(58) = getPowerDeriv(T(13),params(7)/(params(7)-1),1);
T(59) = params(4)*(-T(1))/(y(35)*y(35));
T(60) = getPowerDeriv(y(35)/y(36),T(3),1);
T(61) = getPowerDeriv(y(35)/y(37),T(3),1);

end
