function [y,reskeep,errs1,errs2,iters,ierrmsg]=...
  onestep(m,t0,y0,dt,n,T_hat,V_hat,H_hat,tc,kmax,gtol,etol,A,B)
global num_mult
global sdcint
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  one time step for sdc
%
%  Input:
%    m: the dimension of the problem.
%    t0: starting time.
%    y0: initial value.
%    dt: time step.
%    n: number of quadrature points.
%    tc: The quadrature points.
%    gtol: Stop tolerance for GMRES.
%    etol: stop tolerance for the residual.
%    A: the integration matrix.
%
%  Output:
%    y: the final solution.
%    reskeep: the residual from GMRES.
%    errs1: right hand side error after each iteration.
%    errs2: true error after each iteration.
%    iters: number of iterations for each GMRES.
%
%  Note: one may change the following:
%    predictor.m: to use better predictor.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ierrmsg=0;          % future error message.

t=dt*tc+t0;  % set up the Gaussian points in [t0,t1].
t1=dt+t0;    % t1 is the final time.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% call predictor. 
% ynew is the approximate solution.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ynew=predictor(t0,y0,t,m,n,T_hat,V_hat,H_hat);

% solution at the final time point.
y=zeros([m,m]);
y(:,:)=ynew(n+1,:,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute the residue.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rhsint=tensorprod(dt*A, ynew, 2,1);
eps=zeros([n+1,m,m]);
for k=1:n+1
  temp0=zeros([m,m]); temp0(:,:)=rhsint(k,:,:);
  temp1=zeros([m,m]);temp1(:,:)=ynew(k,:,:);
  eps(k,:,:)=y0-H_hat*temp0-temp1;
  num_mult = num_mult + 1;
end
% epstemp=y0all+tensorprod(rhsint,-H_hat,2,2)-ynew;
% this is the residual defined in the paper.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize monitor vectors.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
errest=max(max(max(abs(eps)))); % estimate error (residue).
lerrest=errest; % used to check if residual decays or not.
lcount=0;       % used to check if the number of nondcay error iterations.
errs=trueerr(y,t1,H_hat); % the true error.
reskeep=[];         % gmres error.
iters=[];           % gmres iterations.
errs1=errest;       % right hand side error, the residual.
errs2=errs;         % the true error.
count1=0;           % the number of GMRES corrections.
error=[];           % each gmres error.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SDC corrections.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while count1<sdcint && errest>etol  %kmax: max number of iterations allowed.
  count1=count1+1;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Find a low order approximation of the error
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  delta=zeros([n+1,m,m]); % the initial guess
  % dt=t(1)-t0;
  % epsmat=zeros([m,m]); epsmat(:,:)=eps(1,:,:);
  % delta(1,:,:)=epsmat+(dt/2)*(-H_hat)*epsmat;
  % for k=2:n+1
  %     dt=t(k)-t(k-1);
  %     epsmat0=zeros([m,m]); 
  %     epsmat0(:,:)=delta(k-1,:,:)-eps(k-1,:,:)+dt/2*eps(k-1,:,:);
  %     epsmat1=zeros([m,m]);
  %     epsmat(:,:)=eps(k,:,:);
  %     epsmat1=epsmat+dt/2*(-H_hat)*epsmat;
  %     delta(k,:,:)=expm(-T_hat*(t(k)-t(k-1)))*expm(-V_hat*(t(k)-t(k-1)))*epsmat0+epsmat1;
  % end
  for k=2:n+1
      dtk=t(k)-t(k-1);
      dmat0=zeros([m,m]); dmat0(:,:)=delta(k-1,:,:);
      epsmat0=zeros([m,m]); epsmat0(:,:)=eps(k-1,:,:);
      epsmat1=zeros([m,m]);epsmat1(:,:)=eps(k,:,:);
      delta(k,:,:)=expm(-H_hat*dtk)*dmat0+epsmat1-expm(-H_hat*dtk)*epsmat0...
          -dtk/2*H_hat*(expm(-H_hat*dtk)*epsmat0+epsmat1);
      num_mult = num_mult + 4;
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % update approximate solution.
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ynew=ynew+delta;                               % the new solution.

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % compute the residue.
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  rhsint=tensorprod(dt*A, ynew, 2,1);
  eps=zeros([n+1,m,m]);
  for k=1:n+1
      temp0=zeros([m,m]); temp0(:,:)=rhsint(k,:,:);
      temp1=zeros([m,m]);temp1(:,:)=ynew(k,:,:);
      eps(k,:,:)=y0-H_hat*temp0-temp1;
      num_mult = num_mult + 1;
  end
  % this is the residual defined in the paper.

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % monitor the errors.
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  errest=max(max(max(abs(eps)))); % right hand side error.
  reskeep=[reskeep error]; % error from gmres.
  errs1=[errs1,errest];  % right hand side error.
  % iters(count1)=total_iters; %
  % errs2(count1)=trueerr(y,t1);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %  find the best solution at the end point.
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % rhsint=tensorprod(dt*B, ynew, [1],[1]);
  % y=y0+tensorprod(rhsint,-H_hat,[1],[2]);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % check how many non-decay iterations have happened.
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if errest>=lerrest*0.5
    if lcount==400
      % ierrmsg=1
      break;
    else
      lcount=lcount+1;
    end
  else
    lerrest=errest; %lcount=0;
  end

end

% after convergence, construct the best solution.

rhsint=tensorprod(dt*B, ynew, 2,1);
temp0=zeros([m,m]); temp0(:,:)=rhsint(1,:,:);
y=y0-H_hat*temp0;

return
