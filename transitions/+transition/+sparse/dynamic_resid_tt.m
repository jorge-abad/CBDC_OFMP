function [T_order, T] = dynamic_resid_tt(y, x, params, steady_state, T_order, T)
if T_order >= 0
    return
end
T_order = 0;
if size(T, 1) < 42
    T = [T; NaN(42 - size(T, 1), 1)];
end
T(1) = y(64)^params(2);
T(2) = T(1)/y(95)*params(4);
T(3) = 1/params(7);
T(4) = (y(95)/y(93))^T(3);
T(5) = (y(95)/y(96))^T(3);
T(6) = (y(95)/y(97))^T(3);
T(7) = y(68)^params(3);
T(8) = y(64)^(-params(2));
T(9) = 1-params(12)*(y(63)/y(3)-1)^2;
T(10) = T(9)-y(63)*(y(63)/y(3)-1)*2*params(12)/y(3);
T(11) = (params(7)-1)/params(7);
T(12) = y(97)^T(11);
T(13) = params(5)*y(93)^T(11)+params(6)*y(96)^T(11)+x(2)*T(12);
T(14) = y(68)^(1-params(8));
T(15) = (y(34)*y(5))^params(8);
T(16) = params(11)*(1+y(71))^params(10);
T(17) = (y(93)+y(66))*params(23)/(1-params(23));
T(18) = y(87)^params(17)+y(86)^params(17);
T(19) = T(18)^(1/params(17));
T(20) = 1+(y(86)/y(87))^params(17);
T(21) = (y(93)+y(66))/(y(66)*(1-params(23)));
T(22) = params(14)*(1-y(106))*y(108)+y(110)*(y(106)-y(107))*T(21);
T(23) = params(14)*(1-y(106))+(y(106)-y(107))*T(21);
T(24) = params(16)^2;
T(25) = exp(params(15)+T(24)/2);
T(26) = (T(24)+params(15)-log(y(89)))/params(16);
T(27) = .5-.5*erf((log(y(89))-params(15))/(params(16)*1.414213562373095));
T(28) = (T(24)+params(15)-log(y(90)))/params(16);
T(29) = .5-.5*erf((log(y(90))-params(15))/(params(16)*1.414213562373095));
T(30) = y(106)-y(107)+(abs(y(89)-y(90))<1e-7);
T(31) = (log(y(89))-params(15))/params(16);
T(32) = (log(y(90))-params(15))/params(16);
T(33) = T(14)*1/y(75);
T(34) = (1-params(8))*y(70)/y(67);
T(35) = T(34)^((1-params(8))/params(8));
T(36) = T(34)^(1/params(8));
T(37) = y(34)*T(36);
T(38) = params(1)*(y(64)+y(174))^(-params(2));
T(39) = 2*params(12)*y(69)*(y(62)+y(175))*((y(63)+y(176))/y(63)-1);
T(40) = ((y(63)+y(176))/y(63))^2;
T(41) = (1+y(71)+y(173))^params(10);
T(42) = (1+y(71)+y(173))^(params(10)-1);
end
