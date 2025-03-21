function T = dynamic_g1_tt(T, y, x, params, steady_state, it_)
% function T = dynamic_g1_tt(T, y, x, params, steady_state, it_)
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

assert(length(T) >= 66);

T = transition.dynamic_resid_tt(T, y, x, params, steady_state, it_);

T(43) = (-(params(12)*(-y(26))/(y(2)*y(2))*2*(y(26)/y(2)-1)));
T(44) = (-(params(12)*2*(y(26)/y(2)-1)*1/y(2)));
T(45) = params(4)*getPowerDeriv(y(27),params(2),1)/y(58);
T(46) = getPowerDeriv(y(27),(-params(2)),1);
T(47) = params(1)*getPowerDeriv(y(27)+y(85),(-params(2)),1);
T(48) = getPowerDeriv(y(19)*y(4),params(8),1);
T(49) = (y(69)-y(70))*(y(29)*(1-params(23))-(y(56)+y(29))*(1-params(23)))/(y(29)*(1-params(23))*y(29)*(1-params(23)));
T(50) = getPowerDeriv(T(34),(1-params(8))/params(8),1);
T(51) = getPowerDeriv(T(34),1/params(8),1);
T(52) = (-((y(36)+y(88))*y(32)*params(11)*getPowerDeriv(1+y(34)+y(84),params(10),1)));
T(53) = (-((y(37)+y(89))*y(32)*params(11)*getPowerDeriv(1+y(34)+y(84),params(10)-1,1)));
T(54) = (-((-(y(42)*(y(39)+y(90))))/((1+y(34)+y(84))*(y(39)+y(90))*(1+y(34)+y(84))*(y(39)+y(90)))));
T(55) = (-((-(y(41)*(y(39)+y(90))))/((1+y(34)+y(84))*(y(39)+y(90))*(1+y(34)+y(84))*(y(39)+y(90)))));
T(56) = (-((-(y(42)*(1+y(34)+y(84))))/((1+y(34)+y(84))*(y(39)+y(90))*(1+y(34)+y(84))*(y(39)+y(90)))));
T(57) = (-((-(y(41)*(1+y(34)+y(84))))/((1+y(34)+y(84))*(y(39)+y(90))*(1+y(34)+y(84))*(y(39)+y(90)))));
T(58) = getPowerDeriv(T(18),1/params(17),1);
T(59) = getPowerDeriv(y(49),params(17),1)*T(58);
T(60) = getPowerDeriv(y(49)/y(50),params(17),1);
T(61) = T(58)*getPowerDeriv(y(50),params(17),1);
T(62) = getPowerDeriv(y(58)/y(56),T(3),1);
T(63) = getPowerDeriv(T(13),params(7)/(params(7)-1),1);
T(64) = params(4)*(-T(1))/(y(58)*y(58));
T(65) = getPowerDeriv(y(58)/y(59),T(3),1);
T(66) = getPowerDeriv(y(58)/y(60),T(3),1);

end
