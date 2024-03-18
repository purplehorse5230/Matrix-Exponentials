function [A,tLegendre]=tgaussint(n);
%A=zeros(n+1,n+1);
eval(strcat('tlgrInt1',num2str(n)));
return