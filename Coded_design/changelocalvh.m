function gdmd=changelocalvh(gdmd,n1,n2,MH,tam_mask)
[N1,N2,S]=size(gdmd);
N=N1;
J=randperm(N1-2*n1)+n1;
aux=randperm(N1);
J=[J aux(1:MH)];
aux=randperm(N2);
H=randperm(N2-2*n2)+n2;
H=[H aux(1:MH)];
for j=1:length(J)
    for h=1:length(H)       
        for r=1:S
            a2=zeros(1,S);
            a2(r)=1;
            gdmd(J(j),H(h),:)=a2(:);
            mask=fspecial('gaussian',tam_mask);
            suma=neigvorintegrado(gdmd,n1,n2,J(j),H(h), mask);
            ya=max(suma);
            ss(r)=ya;
        end        
        [a]=find(ss==min(ss));
        u=randperm(length(a));
        gdmd(J(j),H(h),:)=0;
        gdmd(J(j),H(h),a(u(1)))=1;        
    end
end