function residual = dynamic_resid(T, y, x, params, steady_state, it_, T_flag)
% function residual = dynamic_resid(T, y, x, params, steady_state, it_, T_flag)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T             [#temp variables by 1]     double   vector of temporary terms to be filled by function
%   y             [#dynamic variables by 1]  double   vector of endogenous variables in the order stored
%                                                     in M_.lead_lag_incidence; see the Manual
%   x             [nperiods by M_.exo_nbr]   double   matrix of exogenous variables (in declaration order)
%                                                     for all simulation periods
%   steady_state  [M_.endo_nbr by 1]         double   vector of steady state values
%   params        [M_.param_nbr by 1]        double   vector of parameter values in declaration order
%   it_           scalar                     double   time period for exogenous variables for which
%                                                     to evaluate the model
%   T_flag        boolean                    boolean  flag saying whether or not to calculate temporary terms
%
% Output:
%   residual
%

if T_flag
    T = transition.dynamic_resid_tt(T, y, x, params, steady_state, it_);
end
residual = zeros(60, 1);
    residual(1) = (y(32)*y(46)/(1+y(34)+y(84))) - (1-T(2)*params(5)*T(4));
    residual(2) = (y(32)/(1+y(34)+y(84))) - (1-T(2)*params(6)*T(5));
    residual(3) = (y(32)*params(26)/(1+y(34)+y(84))) - (1-T(2)*x(it_, 2)*T(6));
    residual(4) = (y(30)) - (T(7)/T(8));
    residual(5) = (y(32)) - (T(38)/T(8));
    residual(6) = (1) - (y(25)*T(10)+T(39)*T(40));
    residual(7) = (y(28)) - (y(26)*T(9)+(1-params(9))*y(19)*y(4));
    residual(8) = (y(58)) - (T(13)^(params(7)/(params(7)-1)));
    residual(9) = (y(24)) - (T(15)*T(33));
    residual(10) = (1) - (params(11)*(1+y(34))^(params(10)-1)+(1-params(11))*y(35)^(1-params(10)));
    residual(11) = (y(35)) - (y(36)/y(37));
    residual(12) = (y(36)) - (y(24)*params(10)/(params(10)-1)*y(33)+y(32)*params(11)*T(41)*(y(36)+y(88)));
    residual(13) = (y(37)) - (y(24)+y(32)*params(11)*T(42)*(y(37)+y(89)));
    residual(14) = (y(38)) - ((1-params(11))*y(35)^(-params(10))+T(16)*y(8));
    residual(15) = (y(39)) - ((y(40)+y(25)*(1-params(9)))/y(1));
    residual(16) = (y(40)) - (params(8)*y(33)*T(35));
    residual(17) = (y(31)) - (y(4)*T(37));
    residual(18) = (y(25)*y(28)) - (y(29)*params(14)*(1-y(69))+(y(56)+y(29))/(1-params(23))*(y(69)-y(70)));
    residual(19) = (y(61)) - ((1-y(69))*y(29)*params(14)*params(23)+(y(69)-y(70))*T(17));
    residual(20) = (y(29)) - (params(13)*(y(4)*y(19)*y(39)*y(1)-y(11)*y(14)/(1+y(34))-y(21)*y(20)/(1+y(34))+y(10)*y(15)/(1+y(34))+y(45)*y(17)/(1+y(34))-y(13)*y(18)/(1+y(34))));
    residual(21) = (y(52)) - (y(42)/((1+y(34)+y(84))*(y(39)+y(90))));
    residual(22) = (y(53)) - (y(41)/((1+y(34)+y(84))*(y(39)+y(90))));
    residual(23) = (1/y(74)) - (y(41));
    residual(24) = (y(46)) - ((1-y(69))*y(42)+y(70)*y(41)+y(73)*(1+y(34)+y(84))*(y(69)-y(70))*(y(39)+y(90)));
    residual(25) = (y(49)) - ((1-y(69))*(y(29)*(params(14)*(1-params(23))-1)-y(56)));
    residual(26) = (y(50)) - ((y(56)+y(29))*y(70)-y(55));
    residual(27) = (y(47)) - (y(50)/T(19));
    residual(28) = (y(48)) - (y(49)/T(19));
    residual(29) = (y(42)) - (y(47)*y(51)*y(43)+(1-y(47)*y(51))*y(44));
    residual(30) = (y(41)) - (y(44)*y(48)*(1-y(51))+y(43)*(1-y(48)*(1-y(51))));
    residual(31) = (y(51)) - (1/T(20));
    residual(32) = (y(43)) - (params(19)*(y(12)+(1-y(16))*params(21))+(1-params(19))*(x(it_, 1)+params(18)*(y(34)-params(20)))-(1-y(51))*params(21));
    residual(33) = (y(44)) - (y(43)+params(21));
    residual(34) = (y(62)) - (y(43)-params(22));
    residual(35) = (y(61)+y(54)+y(49)*(1-y(47))) - (y(60)+y(59)+y(50)*(1-y(48)));
    residual(36) = (y(54)) - ((steady_state(1))*params(24)*params(25));
    residual(37) = ((steady_state(1))*params(24)) - (y(55)+y(54));
    residual(38) = (y(45)) - (1/y(23));
    residual(39) = (y(57)) - (T(22)/T(23));
    residual(40) = (y(27)) - (y(24)-y(26));
    residual(41) = (y(71)) - (T(25)*normcdf(T(26),0,1)/T(27));
    residual(42) = (y(72)) - (T(25)*normcdf(T(28),0,1)/T(29));
    residual(43) = (y(73)) - ((y(72)*(1-y(70))-(1-y(69))*y(71))/T(30));
    residual(44) = (y(69)) - (normcdf(T(31),0,1));
    residual(45) = (y(70)) - (normcdf(T(32),0,1));
    residual(46) = (y(68)) - (1);
    residual(47) = (y(75)) - (y(46)*(1+y(34)+y(84))-1);
    residual(48) = (y(63)) - (log(y(27)));
    residual(49) = (y(64)) - (params(4)*log(y(58)));
    residual(50) = (y(65)) - (y(31)^(1+params(3))/(1+params(3)));
    residual(51) = (y(66)) - (y(65)+y(63)+y(64)+params(1)*(y(66)+y(91)));
    residual(52) = (y(67)) - (y(65)+y(63)+y(64));
    residual(53) = (y(76)) - (y(34)-y(5));
    residual(54) = (y(77)) - (y(27)-y(3));
    residual(55) = (y(78)) - (y(25)-y(1));
    residual(56) = (y(79)) - (y(26)-y(2));
    residual(57) = (y(80)) - (y(36)-y(6));
    residual(58) = (y(81)) - (y(37)-y(7));
    residual(59) = (y(82)) - (y(39)-y(9));
    residual(60) = (y(83)) - (y(66)-y(22));

end
