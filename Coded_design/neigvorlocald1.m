function suma=neigvorlocald1(gdmd,n1,n2,j,h)
[N,N,S]=size(gdmd);
suma=0;
u=min(n1,n2);
U=eye(u);
%U2=U(:,end:-1:1);
%U=(U+U2)>0;
U2=zeros(n1,n2);
Mask=U;
if n1>n2
    b=round((n1-n2)/2);
    Z1=zeros(b,n2);
    Z2=zeros(n1,b);
    Mask=[Z1;U;Z1];
   % Mask=[Z2 Mask Z2];
else
    if n2>n1
            b=round((n2-n1)/2);
            Z1=zeros(b,n2);
            Z2=zeros(n1,b);
            Mask=[Z2 U Z2];            
        %    Mask=[Z1;Mask;Z1];
    end
end



n1=n1-1;
n2=n2-1;
Mask=Mask(max(1,-j+n1/2+2):min(n1+1,N-j+n1/2+1),max(1,-h+n2/2+2):min(n2+1,N-h+n2/2+1));

for r=1:S
    p1=max(1,j-n1/2);
    p2=min(N,j+n1/2);
    p3=max(1,h-n2/2);
    p4=min(N,h+n2/2);
    suma(r)=sum(sum(gdmd(p1:p2,p3:p4,r).*Mask));
end

