function ynew=predictor(t0,y0,t,m,n,T_hat,V_hat,H_hat)

% first use constant approximation. This can be improved by a low order
% method later using Rachel's methods.

ichoice=2;
ytemp=zeros([m,m]);
ynew=zeros([n+1,m,m]);
switch ichoice
    case 1
        %constant approximation
        for k=1:n+1
            ynew(k,:,:)=y0;
        end
    case 2
        %low order approximation
        ynew(1,:,:)=expm(-T_hat*(t(1)-t0))*expm(-V_hat*(t(1)-t0))*y0;
        for k=1:n
            ytemp(:,:)=ynew(k,:,:);
            ynew(k+1,:,:)=expm(-T_hat*(t(k+1)-t(k)))*expm(-V_hat*(t(k+1)-t(k)))*ytemp;
        end
    case 3
        % accurate approximation for debugging purposes
        ynew(1,:,:)=expm(-H_hat*(t(1)-t0))*y0;
        for k=1:n
            ytemp(:,:)=ynew(k,:,:);
            ynew(k+1,:,:)=expm(-H_hat*(t(k+1)-t(k)))*ytemp;
        end
end

return