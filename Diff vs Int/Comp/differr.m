for n=3:80
    [A,tc]=tgauss(n);

    fval=sin(tc);
    fdiffanal=(cos(tc))';
    fdiffnum=A*fval';
    %error(n)=norm(fdiffanal(2:end)-fdiffnum(2:end));
    errordiff(n)=norm(fdiffanal(2:end)-fdiffnum(2:end),"inf");
end

semilogy(errordiff)

for n=3:79
    [A,tc]=tgaussint(n);

    fval=sin(tc);
    fintanal=(1-cos(tc))';
    fintnum=A*fval';
    %error(n)=norm(fdiffanal(2:end)-fdiffnum(2:end));
    errint(n)=norm(fintanal(2:end)-fintnum(2:end),"inf");
end

nn=1:80;
semilogy(nn(3:80),errordiff(3:80),'r',nn(3:79),errint(3:79),'b')
xlabel("# of Gaussian Nodes")
ylabel("Approximation Error")
legend('Spectral Differentiation','Spectral Integration')
