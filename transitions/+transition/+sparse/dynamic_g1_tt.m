function [T_order, T] = dynamic_g1_tt(y, x, params, steady_state, T_order, T)
if T_order >= 1
    return
end
[T_order, T] = transition.sparse.dynamic_resid_tt(y, x, params, steady_state, T_order, T);
T_order = 1;
if size(T, 1) < 66
    T = [T; NaN(66 - size(T, 1), 1)];
end
T(43) = (-(params(12)*(-y(63))/(y(3)*y(3))*2*(y(63)/y(3)-1)));
T(44) = (-(params(12)*2*(y(63)/y(3)-1)*1/y(3)));
T(45) = params(4)*getPowerDeriv(y(64),params(2),1)/y(95);
T(46) = getPowerDeriv(y(64),(-params(2)),1);
T(47) = params(1)*getPowerDeriv(y(64)+y(174),(-params(2)),1);
T(48) = getPowerDeriv(y(34)*y(5),params(8),1);
T(49) = (y(106)-y(107))*(y(66)*(1-params(23))-(y(93)+y(66))*(1-params(23)))/(y(66)*(1-params(23))*y(66)*(1-params(23)));
T(50) = getPowerDeriv(T(34),(1-params(8))/params(8),1);
T(51) = getPowerDeriv(T(34),1/params(8),1);
T(52) = (-((y(73)+y(177))*y(69)*params(11)*getPowerDeriv(1+y(71)+y(173),params(10),1)));
T(53) = (-((y(74)+y(178))*y(69)*params(11)*getPowerDeriv(1+y(71)+y(173),params(10)-1,1)));
T(54) = (-((-(y(79)*(y(76)+y(179))))/((1+y(71)+y(173))*(y(76)+y(179))*(1+y(71)+y(173))*(y(76)+y(179)))));
T(55) = (-((-(y(78)*(y(76)+y(179))))/((1+y(71)+y(173))*(y(76)+y(179))*(1+y(71)+y(173))*(y(76)+y(179)))));
T(56) = (-((-(y(79)*(1+y(71)+y(173))))/((1+y(71)+y(173))*(y(76)+y(179))*(1+y(71)+y(173))*(y(76)+y(179)))));
T(57) = (-((-(y(78)*(1+y(71)+y(173))))/((1+y(71)+y(173))*(y(76)+y(179))*(1+y(71)+y(173))*(y(76)+y(179)))));
T(58) = getPowerDeriv(T(18),1/params(17),1);
T(59) = getPowerDeriv(y(86),params(17),1)*T(58);
T(60) = getPowerDeriv(y(86)/y(87),params(17),1);
T(61) = T(58)*getPowerDeriv(y(87),params(17),1);
T(62) = getPowerDeriv(y(95)/y(93),T(3),1);
T(63) = getPowerDeriv(T(13),params(7)/(params(7)-1),1);
T(64) = params(4)*(-T(1))/(y(95)*y(95));
T(65) = getPowerDeriv(y(95)/y(96),T(3),1);
T(66) = getPowerDeriv(y(95)/y(97),T(3),1);
end
