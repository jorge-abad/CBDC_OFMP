function T = dynamic_resid_tt(T, y, x, params, steady_state, it_)
% function T = dynamic_resid_tt(T, y, x, params, steady_state, it_)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T             [#temp variables by 1]     double  vector of temporary terms to be filled by function
%   y             [#dynamic variables by 1]  double  vector of endogenous variables in the order stored
%                                                    in M_.lead_lag_incidence; see the Manual
%   x             [nperiods by M_.exo_nbr]   double  matrix of exogenous variables (in declaration order)
%                                                    for all simulation periods
%   steady_state  [M_.endo_nbr by 1]         double  vector of steady state values
%   params        [M_.param_nbr by 1]        double  vector of parameter values in declaration order
%   it_           scalar                     double  time period for exogenous variables for which
%                                                    to evaluate the model
%
% Output:
%   T           [#temp variables by 1]       double  vector of temporary terms
%

assert(length(T) >= 42);

T(1) = y(27)^params(2);
T(2) = T(1)/y(58)*params(4);
T(3) = 1/params(7);
T(4) = (y(58)/y(56))^T(3);
T(5) = (y(58)/y(59))^T(3);
T(6) = (y(58)/y(60))^T(3);
T(7) = y(31)^params(3);
T(8) = y(27)^(-params(2));
T(9) = 1-params(12)*(y(26)/y(2)-1)^2;
T(10) = T(9)-y(26)*(y(26)/y(2)-1)*2*params(12)/y(2);
T(11) = (params(7)-1)/params(7);
T(12) = y(60)^T(11);
T(13) = params(5)*y(56)^T(11)+params(6)*y(59)^T(11)+x(it_, 2)*T(12);
T(14) = y(31)^(1-params(8));
T(15) = (y(19)*y(4))^params(8);
T(16) = params(11)*(1+y(34))^params(10);
T(17) = (y(56)+y(29))*params(23)/(1-params(23));
T(18) = y(50)^params(17)+y(49)^params(17);
T(19) = T(18)^(1/params(17));
T(20) = 1+(y(49)/y(50))^params(17);
T(21) = (y(56)+y(29))/(y(29)*(1-params(23)));
T(22) = params(14)*(1-y(69))*y(71)+y(73)*(y(69)-y(70))*T(21);
T(23) = params(14)*(1-y(69))+(y(69)-y(70))*T(21);
T(24) = params(16)^2;
T(25) = exp(params(15)+T(24)/2);
T(26) = (T(24)+params(15)-log(y(52)))/params(16);
T(27) = .5-.5*erf((log(y(52))-params(15))/(params(16)*1.414213562373095));
T(28) = (T(24)+params(15)-log(y(53)))/params(16);
T(29) = .5-.5*erf((log(y(53))-params(15))/(params(16)*1.414213562373095));
T(30) = y(69)-y(70)+(abs(y(52)-y(53))<1e-7);
T(31) = (log(y(52))-params(15))/params(16);
T(32) = (log(y(53))-params(15))/params(16);
T(33) = T(14)*1/y(38);
T(34) = (1-params(8))*y(33)/y(30);
T(35) = T(34)^((1-params(8))/params(8));
T(36) = T(34)^(1/params(8));
T(37) = y(19)*T(36);
T(38) = params(1)*(y(27)+y(85))^(-params(2));
T(39) = 2*params(12)*y(32)*(y(25)+y(86))*((y(26)+y(87))/y(26)-1);
T(40) = ((y(26)+y(87))/y(26))^2;
T(41) = (1+y(34)+y(84))^params(10);
T(42) = (1+y(34)+y(84))^(params(10)-1);

end
