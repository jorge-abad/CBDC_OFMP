function T = static_resid_tt(T, y, x, params)
% function T = static_resid_tt(T, y, x, params)
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

assert(length(T) >= 39);

T(1) = y(4)^params(2);
T(2) = T(1)/y(35)*params(4);
T(3) = 1/params(7);
T(4) = (y(35)/y(33))^T(3);
T(5) = (y(35)/y(36))^T(3);
T(6) = (y(35)/y(37))^T(3);
T(7) = y(8)^params(3);
T(8) = y(4)^(-params(2));
T(9) = params(1)*(y(4)+y(54))^(-params(2));
T(10) = 2*params(12)*y(9)*(y(2)+y(55))*((y(3)+y(56))/y(3)-1);
T(11) = ((y(3)+y(56))/y(3))^2;
T(12) = (params(7)-1)/params(7);
T(13) = params(5)*y(33)^T(12)+params(6)*y(36)^T(12)+x(2)*y(37)^T(12);
T(14) = (y(5)*y(34))^params(8);
T(15) = y(8)^(1-params(8));
T(16) = T(15)*1/y(15);
T(17) = (1+y(11)+y(53))^params(10);
T(18) = (1+y(11)+y(53))^(params(10)-1);
T(19) = params(11)*(1+y(11))^params(10);
T(20) = (1-params(8))*y(10)/y(7);
T(21) = T(20)^((1-params(8))/params(8));
T(22) = T(20)^(1/params(8));
T(23) = y(34)*T(22);
T(24) = (y(33)+y(6))*params(23)/(1-params(23));
T(25) = y(27)^params(17)+y(26)^params(17);
T(26) = T(25)^(1/params(17));
T(27) = 1+(y(26)/y(27))^params(17);
T(28) = (y(33)+y(6))/(y(6)*(1-params(23)));
T(29) = params(14)*(1-y(46))*y(48)+y(50)*(y(46)-y(47))*T(28);
T(30) = params(14)*(1-y(46))+(y(46)-y(47))*T(28);
T(31) = params(16)^2;
T(32) = exp(params(15)+T(31)/2);
T(33) = (T(31)+params(15)-log(y(29)))/params(16);
T(34) = .5-.5*erf((log(y(29))-params(15))/(params(16)*1.414213562373095));
T(35) = (T(31)+params(15)-log(y(30)))/params(16);
T(36) = .5-.5*erf((log(y(30))-params(15))/(params(16)*1.414213562373095));
T(37) = y(46)-y(47)+(abs(y(29)-y(30))<1e-7);
T(38) = (log(y(29))-params(15))/params(16);
T(39) = (log(y(30))-params(15))/params(16);

end
