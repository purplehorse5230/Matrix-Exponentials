function err=trueerr(y,t,H_hat)

% t is a scalar time.
% y is the numerical solution at time t.
%evaluate the true err at t if you know the eact solution at t,
%otherwise let err=10 (a dummy value) is ok.

err=max(max(abs(y-expm(-H_hat*t))));


return;