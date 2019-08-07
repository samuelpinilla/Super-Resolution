function gdmd=gdmddesign(N,shots)
gdmd = DMD_Boolean(N,N,shots);    % Coded aperture Random entries
iter=100;
sh=1;
Sh1=2;
Sh2=1;
iper=0;
for h=1:10
    if iper==1 % Used only when iper=1 -- Not used this time
        for k=1:iter %Optimization per rows
            p=randperm(shots);
            gdmd2=blueoptimal(gdmd,N,shots,1);
            gdmd=changepos(gdmd,gdmd2(:,:,p(1)),p(1));
        end
        for k=1:iter %Optimization per columns
            p=randperm(shots);
            gdmd2=blueoptimalcolumns(gdmd,N,shots,1);
            gdmd=changepos(gdmd,gdmd2(:,:,p(1)),p(1));
        end
    end
    p=randperm(shots);
    gdmd2=blueoptimalcolumns(gdmd,N,shots,1);
    gdmd=changepos2(gdmd,gdmd2(:,:,p(1)),p(1),Sh1,Sh2);
end