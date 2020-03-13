function [A,r,Tp,tl] = LV_M(c,N,cp,fcp,mu,fmu,cm,fcm,am,fam,dnt,an,fan,rit,mit,fas,mnti,nb,sg,sgd,tsp,msp,max_r,min_mort,it,bc)
[Tr,b]=Tr_matrix(c,N);
C=NT_Comunity(Tr,cp,fcp,mu,fmu,cm,fcm,am,fam,dnt,an,fan,rit,mit,fas,mnti,nb,sg,sgd,tsp,msp,max_r,min_mort,true,100,99.99999999999,bc);
for l=1:it
    
    [Tr,b]=Tr_matrix(c,N);
    C=NT_Comunity(Tr,cp,fcp,mu,fmu,cm,fcm,am,fam,dnt,an,fan,rit,mit,fas,mnti,nb,sg,sgd,tsp,msp,max_r,min_mort,true,100,99.99999999999,bc);
    if(C.est>0)
        break;
    end
    l
end
% [row,col] = find(C.adj);
A=C.com;
r=C.R;
Tp=C.tp;
tl=C.Nv;
% E=zeros(N,N,N);
% 
% for i=1:length(row)
%     for j=length(col)
%         if(rand(1) < hp)
%             E(row(i),col(j),randi([1 N],1))=hif*(rand(1)*2-1);
%         end
%     end
% end

end