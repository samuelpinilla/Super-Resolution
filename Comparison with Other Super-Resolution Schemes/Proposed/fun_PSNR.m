function psnr = fun_PSNR(img,res)

[M,N,L]=size(img);
mn=0;
for i=1:L
temp=1./(M*N)*sum(sum((img(:,:,i)-res(:,:,i)).^2));
psnr=10*log10(max(max(img(:,:,i).^2./temp)));
mn=mn+psnr;
end
psnr=mn/i;