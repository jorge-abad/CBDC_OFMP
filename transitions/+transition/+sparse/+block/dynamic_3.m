function [y, T] = dynamic_3(y, x, params, steady_state, sparse_rowval, sparse_colval, sparse_colptr, T)
  y(101)=params(4)*log(y(95));
  y(100)=log(y(64));
  y(102)=y(68)^(1+params(3))/(1+params(3));
end
