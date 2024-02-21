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
%    kmax, gtol, etol, n, h0, k0.
%
%  One may also change: 
%     predictor.m: which predictor to use.
%     updated.m: which preconditioner to use.
%     chebnodes.m: which node points to use.
%     lcount==4 (see onestep.m)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kmax=100;            %  maximum number of outer GMRES iterations
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for nk = 1:1,
  % n =10+(nk-1)*14;   % loop over different number of nodes for each big time step.
  n=20;
  
  for k =1:1,        % loop over different big time step sizes h0.

    h0=tfinal*0.5^(k+1);  % different big time step size
    h0=0.1;
    
    % figure(3*nk -2); % output figures.
    % subplot(2,2,k);
    % figure(3*nk-1);
    % subplot(2,2,k);
    % figure(3*nk);
    % subplot(2,2,k);
    
    for j = 7:7        % loop over different terms in Krylov subspace method k0
      k0=max(n+4,81);  % Number of terms in Krylov subspace for the 
                       % restarted GMRES(k0).
                       % Note: instead of GMRES(k0), BiCGStab or
                       %  TFQMR may be applied, which should be more 
                       %  efficient and require less storage.

%%%%%%%%%%%%%%%%%%%%%%Main Subroutine%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      mm=size(H_hat); 
      m=mm(1); 
      y0=eye(m); % initial value is the identity matrix, not a vector in classical ODE IVP.
      [ysolfin,res,indres,errrhs,errtrue,inderr,iter,ierrmsg]=...
          sdckdc(m,t0,tfinal,y0,h0,n,T_hat,V_hat,H_hat,kmax,gtol,etol,k0);
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
%    kmax: maximal number that the gmres is done
%    gtol: tolerance for gmres
%    etol: tolerance for err/res. Note that this is for right hand side.
%    k0: gmresk selecting parameter
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ierrmsg 
      iter;
      % check how many digits we have.
      maxerr=  max(max(abs(ysolfin-analsolu(tfinal,H_hat))),1e-60);
      ndigits=log10(maxerr)
      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Output:
%
%  figure 1: gmres error. Number of digits, and number of function
%            evaluations.
%  figure 2: right hand error.
%  figure 3: true error.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      if (ierrmsg==0),
        figure(3*nk-2);
        semilogy(res,scolor);
        title(strcat('j: ',int2str(length(iter)*(n+1)),...
          '  f: ',int2str(length(res)*n),'  e: ',int2str(ndigits)));
        hold on;
        
        figure(3*nk-1);
        semilogy(errrhs,scolor);
        title(strcat('rhs error  n: ',int2str(n)));
        hold on;
        
        figure(3*nk);
        semilogy(errtrue,scolor);
        hold on;
        title(strcat('true error  n: ',int2str(n)));
        
      end

    end
  end
end

%%%%%%%%%%%%%%%%%%%%The End***********************************************

