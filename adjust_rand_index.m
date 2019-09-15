function ari = adjust_rand_index(Y, Y0)

K=max(Y); 
K0=max(Y0);
nk=zeros(K,K0);
for i=1:K
    for j=1:K0
        nk(i,j)=sum((Y==i)&(Y0==j));
    end;
end;
sums=0;
for i=1:K
    for j=1:K0
        sums=sums+nk(i,j)*(nk(i,j)-1)/2;
    end;
end;
nk1=sum(nk,1);
sumd1=0;
for j=1:K0
    sumd1=sumd1+nk1(j)*(nk1(j)-1)/2;
end;  
nk2=sum(nk,2);
sumd2=0;
for i=1:K
    sumd2=sumd2+nk2(i)*(nk2(i)-1)/2;
end; 
N=numel(Y); sumN=N*(N-1)/2;
ari=(sums-sumd1*sumd2/sumN)/(0.5*(sumd1+sumd2)-sumd1*sumd2/sumN);

