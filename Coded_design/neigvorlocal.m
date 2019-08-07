function [suma,aux]=neigvorlocal(gdmd,n1,n2,j,h)
[N,N,S]=size(gdmd);
suma=0;
n1=n1-1;
n2=n2-1;
for r=1:S
    p1=max(1,j-n1/2);
    p2=min(N,j+n1/2);
    p3=max(1,h-n2/2);
    p4=min(N,h+n2/2);
    u=gdmd(p1:p2,p3:p4,r);
    du=bwmorph(u,'clean'); %Remove isolated pixels in the window
    suma(r)=sum((u(:)));
    aux(r)=sum(du(:));
end

