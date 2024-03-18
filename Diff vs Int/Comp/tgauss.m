function [A,tLegendre]=tgauss(n);
%A=zeros(n+1,n+1);
eval(strcat('tlgr1',num2str(n)));
return