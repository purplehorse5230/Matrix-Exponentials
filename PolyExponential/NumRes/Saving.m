save("errorn10.mat","n_err","n_h","n_fun")


loglog(n_h(2,:),n_err(2,:),'r',n_h(3,:),n_err(3,:),'b',n_h(4,:),n_err(4,:),'g')