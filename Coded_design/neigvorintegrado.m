function suma=neigvorintegrado(gdmd,n1,n2,j,h,mask)
[M,N,S]=size(gdmd);
suma=zeros(1,S);
for k=1:S
   suma(1,k)=sum(sum(gdmd(j-n1:j+n1,h-n2:h+n2,k).*mask)); 
end

end