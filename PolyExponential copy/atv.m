function ax = atv(x,t,dt,A,dfdy,n,m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This is the matrix vector product subroutine.
%
% Input:
%   x: the input vector.
%   t: the points where we want the solution.
%   dt: the step size.
%   A: the integration matrix.
%   dfdy: the Jacobian matrix.
%   n: the number of Gaussian points.
%
% Question: (Jingfang)
%   Are the Jacobian matrix the same?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global iformulation;

switch iformulation
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Note: for the Picard formulation, one first multiple the Jacobian
  %       matrix, then the integration matrix.
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  case 0           % The Picard formulation.
    x=reshape(x,n+1,m);
    for i=1:n+1,
        dum0(i,:)=x(i,:)*dfdy(:,:,i);
    end
    
    dum=dt*A*dum0;

    rhs=x-dum;
    delta=updated(rhs,t,dfdy,n,m);  % Apply the preconditioner.

    ax=delta(:); % Output the results.
    
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Note: for the yp formulation, one first multipole the integration matrix,
  %       then the Jacobian matrix.
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  case 1          % the yp formulation.
    x=reshape(x,n+1,m);
    dum=dt*A*x;

    for i=1:n+1
      dum(i,:)=dum(i,:)*dfdy(:,:,i); %Note: the Jacobian matrix is written so we have
      % a right matrix vector product.
    end

    rhs=x-dum;
    delta=updated(rhs,t,dfdy,n,m);  % Apply the preconditioner.

    ax=delta(:); % Output the results.
  otherwise
    stop
end

return