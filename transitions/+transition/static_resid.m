function residual = static_resid(T, y, x, params, T_flag)
% function residual = static_resid(T, y, x, params, T_flag)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T         [#temp variables by 1]  double   vector of temporary terms to be filled by function
%   y         [M_.endo_nbr by 1]      double   vector of endogenous variables in declaration order
%   x         [M_.exo_nbr by 1]       double   vector of exogenous variables in declaration order
%   params    [M_.param_nbr by 1]     double   vector of parameter values in declaration order
%                                              to evaluate the model
%   T_flag    boolean                 boolean  flag saying whether or not to calculate temporary terms
%
% Output:
%   residual
%

if T_flag
    T = transition.static_resid_tt(T, y, x, params);
end
residual = zeros(60, 1);
    residual(1) = (y(9)*y(23)/(1+y(11)+y(53))) - (1-T(2)*params(5)*T(4));
    residual(2) = (y(9)/(1+y(11)+y(53))) - (1-T(2)*params(6)*T(5));
    residual(3) = (y(9)*params(26)/(1+y(11)+y(53))) - (1-T(2)*x(2)*T(6));
    residual(4) = (y(7)) - (T(7)/T(8));
    residual(5) = (y(9)) - (T(9)/T(8));
    residual(6) = (1) - (y(2)+T(10)*T(11));
    residual(7) = (y(5)) - (y(3)+y(5)*(1-params(9))*y(34));
    residual(8) = (y(35)) - (T(13)^(params(7)/(params(7)-1)));
    residual(9) = (y(1)) - (T(14)*T(16));
    residual(10) = (1) - (params(11)*(1+y(11))^(params(10)-1)+(1-params(11))*y(12)^(1-params(10)));
    residual(11) = (y(12)) - (y(13)/y(14));
    residual(12) = (y(13)) - (y(1)*params(10)/(params(10)-1)*y(10)+y(9)*params(11)*T(17)*(y(13)+y(57)));
    residual(13) = (y(14)) - (y(1)+y(9)*params(11)*T(18)*(y(14)+y(58)));
    residual(14) = (y(15)) - ((1-params(11))*y(12)^(-params(10))+y(15)*T(19));
    residual(15) = (y(16)) - ((y(17)+y(2)*(1-params(9)))/y(2));
    residual(16) = (y(17)) - (params(8)*y(10)*T(21));
    residual(17) = (y(8)) - (y(5)*T(23));
    residual(18) = (y(2)*y(5)) - (y(6)*params(14)*(1-y(46))+(y(33)+y(6))/(1-params(23))*(y(46)-y(47)));
    residual(19) = (y(38)) - ((1-y(46))*y(6)*params(14)*params(23)+(y(46)-y(47))*T(24));
    residual(20) = (y(6)) - (params(13)*(y(5)*y(34)*y(2)*y(16)-y(19)*y(26)/(1+y(11))-y(38)*y(39)/(1+y(11))+y(18)*y(27)/(1+y(11))+y(22)*y(32)/(1+y(11))-y(23)*y(33)/(1+y(11))));
    residual(21) = (y(29)) - (y(19)/((1+y(11)+y(53))*(y(16)+y(59))));
    residual(22) = (y(30)) - (y(18)/((1+y(11)+y(53))*(y(16)+y(59))));
    residual(23) = (1/y(51)) - (y(18));
    residual(24) = (y(23)) - ((1-y(46))*y(19)+y(47)*y(18)+y(50)*(1+y(11)+y(53))*(y(46)-y(47))*(y(16)+y(59)));
    residual(25) = (y(26)) - ((1-y(46))*(y(6)*(params(14)*(1-params(23))-1)-y(33)));
    residual(26) = (y(27)) - ((y(33)+y(6))*y(47)-y(32));
    residual(27) = (y(24)) - (y(27)/T(26));
    residual(28) = (y(25)) - (y(26)/T(26));
    residual(29) = (y(19)) - (y(24)*y(28)*y(20)+(1-y(24)*y(28))*y(21));
    residual(30) = (y(18)) - (y(21)*y(25)*(1-y(28))+y(20)*(1-y(25)*(1-y(28))));
    residual(31) = (y(28)) - (1/T(27));
    residual(32) = (y(20)) - (params(19)*(y(20)+(1-y(28))*params(21))+(1-params(19))*(x(1)+params(18)*(y(11)-params(20)))-(1-y(28))*params(21));
    residual(33) = (y(21)) - (y(20)+params(21));
    residual(34) = (y(39)) - (y(20)-params(22));
    residual(35) = (y(38)+y(31)+y(26)*(1-y(24))) - (y(37)+y(36)+y(27)*(1-y(25)));
    residual(36) = (y(31)) - ((y(1))*params(24)*params(25));
    residual(37) = ((y(1))*params(24)) - (y(32)+y(31));
    residual(38) = (y(22)) - (1/y(51));
    residual(39) = (y(34)) - (T(29)/T(30));
    residual(40) = (y(4)) - (y(1)-y(3));
    residual(41) = (y(48)) - (T(32)*normcdf(T(33),0,1)/T(34));
    residual(42) = (y(49)) - (T(32)*normcdf(T(35),0,1)/T(36));
    residual(43) = (y(50)) - ((y(49)*(1-y(47))-(1-y(46))*y(48))/T(37));
    residual(44) = (y(46)) - (normcdf(T(38),0,1));
    residual(45) = (y(47)) - (normcdf(T(39),0,1));
    residual(46) = (y(45)) - (1);
    residual(47) = (y(52)) - (y(23)*(1+y(11)+y(53))-1);
    residual(48) = (y(40)) - (log(y(4)));
    residual(49) = (y(41)) - (params(4)*log(y(35)));
    residual(50) = (y(42)) - (y(8)^(1+params(3))/(1+params(3)));
    residual(51) = (y(43)) - (y(42)+y(40)+y(41)+params(1)*(y(43)+y(60)));
    residual(52) = (y(44)) - (y(42)+y(40)+y(41));
residual(53) = y(53);
residual(54) = y(54);
residual(55) = y(55);
residual(56) = y(56);
residual(57) = y(57);
residual(58) = y(58);
residual(59) = y(59);
residual(60) = y(60);

end
