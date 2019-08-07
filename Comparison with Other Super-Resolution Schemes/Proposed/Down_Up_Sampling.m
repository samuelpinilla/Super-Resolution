function [z_upsampled,z_downsampled] = Down_Up_Sampling(z,N,M)
[yN,xN,L]=size(z);

for ll = 1:L
   for sy=1:yN/N;%% %%N %% yN*qq/N
       for sx=1:xN/M;
 temp0=z(( sy-1)*N+(1:N),(sx-1)*M+(1:M),ll);
 z_upsampled(( sy-1)*N+(1:N),(sx-1)*M+(1:M),ll)=ones(N,M)*mean(temp0(:));
 z_downsampled(sy,sx,ll)=sum(temp0(:));
 end
    end
end
end

