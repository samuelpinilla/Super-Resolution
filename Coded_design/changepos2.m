function     gdmd=changepos2(gdmd,gdmd2,shot,Sh1,Sh2)
[N,N,S]=size(gdmd);
[a,b]=find(gdmd2);
pos=[];
for j=1:S
    if j~=shot
        pos=[pos j];
    end
end
for j=1:length(a);
    p=randperm(length(pos));
    for k=1:length(pos)
        if sum(gdmd(a(j),max(b(j)-Sh1,1):min(b(j)+Sh2,N),pos(p(k))))<=1
           gdmd(a(j),b(j),pos(p(k)))=1;
           gdmd(a(j),b(j),shot)=0;
           break;
        end
    end

    
end

    
    



