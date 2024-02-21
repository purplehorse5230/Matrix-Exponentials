function delta=updated(rhs,t,dfdy,n,m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Given right hand side, this subroutine returns
% (1-h*tilde{A}*alpha)^(-1)*rhs.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global iformulation;

switch iformulation;
  case 0
    delta=zeros(n+1,m);
    M=eye(m);

    for i=1:n,
      delta(i+1,:)=(delta(i,:)+rhs(i+1,:)-rhs(i,:))/(M-(t(i+1)-t(i))*dfdy(:,:,i+1));
    end

  case 1
    E=eye(m);
    delta=zeros(n+1,m);
    sdelta=zeros(1,m);

    for i=1:n,
      dum=dfdy(:,:,i+1);
      delta(i+1,:)=(rhs(i+1,:)+sdelta*dum)/(E-(t(i+1)-t(i))*dum);
      sdelta=sdelta+(t(i+1)-t(i))*delta(i+1,:);
    end
  otherwise
end
return
