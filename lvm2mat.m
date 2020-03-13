function lvm2mat(N,c,an,cp,mu,cm,am,ind)

[A,r,Tp,tl] = LV_M(c,N,cp,1,mu,1,cm,1,am,1,@d_c,an,1,0.2,0,1,1,1,0.2,0.2,0,1,1,0,500,false);

fname=sprintf('Lvm_N(%d)c(%d)an(%d)cp(%d)mu(%d)cm(%d)an(%d)ind(%d)',N,c,an,cp,mu,cm,an,ind);
save(strcat(fname,'.mat'),'A','r','Tp','tl','N','an','cp','mu','cm','am','ind');

end

