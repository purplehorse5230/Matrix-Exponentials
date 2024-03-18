%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initiallization. Note: You need to add your table path.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; format short e; clc; close all;
addpath 'GaussNode'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construct the T_hat and V_hat matrices.
%   Code from Rachel.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% T_hat, V_hat construction
MATRIX_NX = 2;
MATRIX_M = 1;
MATRIX_G = 1;
MATRIX_LAMBDA = 0.2;

T_hat = hat_T_3d_two_body_constructor(MATRIX_NX, MATRIX_M);
V_hat = hat_V_3d_two_body_constructor(MATRIX_NX, MATRIX_G * MATRIX_LAMBDA);
H_hat=T_hat+V_hat;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Parameter settings. 
%  One can change the following parameters:
%    etol, n, h0, k0.
%
%  One may also change: 
%     predictor.m: which predictor to use.
%     updated.m: which preconditioner to use.
%     chebnodes.m: which node points to use.
%     lcount==4 (see onestep.m)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kmax=20;            %  maximum number of SDC corrections.
gtol=1e-2;           %  Tolerance for GMRES call(original 1e-11)
                     %  for linear problems, this is
                     %  the total stopping procedure
                     %  For nonlinear problems, it is the
                     %  tolerance for each Newton iteration
etol=1e-14;          %  Total error tolerance for SDC step(original 1e-13)
                     %  the residual in Picard form must be 
                     %  smaller than etol
t0=0;                % to is the starting time.
tfinal=10;            %  tfinal is the final time.
global num_mult
global sdcint
num_mult = 0;
n_h=zeros([10,25]);
n_fun=zeros([10,25]);
n_err=zeros([10,25]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for nk = 2:4
  % n =10+(nk-1)*14;   % loop over different number of nodes for each big time step.
  n=30;
  
  for k =1:25        % loop over different big time step sizes h0.

    h0=tfinal*0.5^(k+1);  % different big time step size
    
    % figure(3*nk -2); % output figures.
    % subplot(2,2,k);
    % figure(3*nk-1);
    % subplot(2,2,k);
    % figure(3*nk);
    % subplot(2,2,k);
    
%%%%%%%%%%%%%%%%%%%%%%Main Subroutine%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    mm=size(H_hat); 
    m=mm(1); 
    y0=eye(m); % initial value is the identity matrix, not a vector in classical ODE IVP.
    [ysolfin,res,indres,errrhs,errtrue,inderr,iter,ierrmsg]=...
          sdckdc(m,t0,tfinal,y0,h0,n,T_hat,V_hat,H_hat,kmax,gtol,etol);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The required parameters of sdcgmres are explained below.
% Output
%    ysolfin: the final solution.
%    res: residue output by GMRES.
%    indres (row vector) index for res : res(indres(i-1)+1:indres(i)), 
%        residues outputs from gmres' during the i_th time step
%    errrhs: the error after each GMRES acceleration.    
%    errtrue: the true error after each gmres acceleration.
%    inderr: (row vector) index for err: err(inderr(i-1)+1:inderr(i)),
%        errors outputs from gmres' during the i_th time step
%    iter (row vector) iteration numbers outputed by each gmres.
%    ierrmsg: error message index.
% 
% Input:
%    m: size of the problemdex.
%    tfinal:  fintal time
%    y0: (row vector) initial value for y
%    h0:  step size
%    n: number of grid points used for each time step
%    gtol: tolerance for gmres
%    etol: tolerance for err/res. Note that this is for right hand side.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    h0
    num_mult
    n_h(nk,k)=h0;
    n_fun(nk,k)=num_mult;
    num_mult=0;
      % check how many digits we have.
    maxerr=  max(max(abs(ysolfin-analsolu(tfinal,H_hat))))
    n_err(nk,k) =maxerr;
    ndigits=log10(maxerr);
      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Output:
%
%  figure 1: gmres error. Number of digits, and number of function
%            evaluations.
%  figure 2: right hand error.
%  figure 3: true error.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % if (ierrmsg==0),
      %   figure(3*nk-2);
      %   semilogy(res,scolor);
      %   title(strcat('j: ',int2str(length(iter)*(n+1)),...
      %     '  f: ',int2str(length(res)*n),'  e: ',int2str(ndigits)));
      %   hold on;
      % 
      %   figure(3*nk-1);
      %   semilogy(errrhs,scolor);
      %   title(strcat('rhs error  n: ',int2str(n)));
      %   hold on;
      % 
      %   figure(3*nk);
      %   semilogy(errtrue,scolor);
      %   hold on;
      %   title(strcat('true error  n: ',int2str(n)));
      % 
      % end

  end
end

%%%%%%%%%%%%%%%%%%%%The End***********************************************

