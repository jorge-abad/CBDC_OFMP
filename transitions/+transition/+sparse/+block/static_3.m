function [y, T] = static_3(y, x, params, sparse_rowval, sparse_colval, sparse_colptr, T)
  y(42)=y(8)^(1+params(3))/(1+params(3));
  y(41)=params(4)*log(y(35));
  y(40)=log(y(4));
  y(44)=y(42)+y(40)+y(41);
end
