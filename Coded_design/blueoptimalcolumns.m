function gdmd2=blueoptimalcolumns(gdmd,N,shots,sh)
    for j=1:shots
        gdmd2(:,:,j)=[[gdmd(:,:,j); zeros(sh,N)].*[zeros(sh,N); gdmd(:,:,j)]];
    end
    gdmd2=gdmd2(2:2+N-1,:,:);
%     gdmd2=zeros(N,N,shots);
%     gdmd2(:,1:N+sh+1,:)=gdmd(:,-sh:N,:);
%     for j=1:shots
%         gdmd2(:,1:N,j)=gdmd(:,:,j).*gdmd2(:,1:N,j);
%     end
%    


h=1;

