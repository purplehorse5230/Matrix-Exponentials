% save("errorn10ex6.mat","n_err","n_h","n_fun")


% loglog(n_h(2,:),n_err(2,:),'r',n_h(3,:),n_err(3,:),'b',n_h(4,:),n_err(4,:),'g')
% title('How errors decay for different time-step sizes (Standard Matrix n=10)')
% xlabel('Step-size')
% ylabel('Approximation Error')
% legend('2 corrections','3 corrections','4 corrections')
% saveas(gcf,'errordecay.png')
% 
% figure(3)
% loglog(n_fun(2,:),n_err(2,:),'r',n_fun(3,:),n_err(3,:),'b',n_fun(4,:),n_err(4,:),'g')
% title('How errors decay for different function evaluations (Standard Matrix n=10)')
% xlabel('Function Evaluations')
% ylabel('Approximation Error')
% legend('2 corrections','3 corrections','4 corrections')
% saveas(gcf,'errordecaydifffunction.png')
% 
% figure(4)
% nn=8
% loglog(n_h(2,1:nn),n_err(2,1:nn),'r',n_h(3,1:nn),n_err(3,1:nn),'b',n_h(4,1:nn),n_err(4,1:nn),'g')
% title('How errors decay for different time-step sizes (High-Rank matrix n=4)')
% xlabel('Step-size')
% ylabel('Approximation Error')
% legend('2 corrections','3 corrections','4 corrections')
% saveas(gcf,'errordecaytimestep3n10.png')
% 
% figure(3)
% nn=8
% loglog(n_fun(2,1:nn),n_err(2,1:nn),'r',n_fun(3,1:nn),n_err(3,1:nn),'b',n_fun(4,1:nn),n_err(4,1:nn),'g')
% title('How errors decay for different function evaluations (High Rank matrix n=4)')
% xlabel('Function Evaluations')
% ylabel('Approximation Error')
% legend('2 corrections','3 corrections','4 corrections')
% saveas(gcf,'errordecaydifffunction5n4.png')
% 
