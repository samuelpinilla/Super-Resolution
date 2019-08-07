function [z,z1,Relerrs] = PRSF_Super_V2(x,y, Params, A, At,z)   
 

x1 = zeros(size(z));
Relerrs=[];
intervalo= Params.itervalo;
  filt = zeros(size(z));
   filt(intervalo,intervalo) =1;

for j = 1 : 20
    j
   u = Params.u0;
    for t = 1: 20,
        grad = compute_grad_super_v2(z,u,y,Params,A,At,x1);
        z = z-Params.mu*grad;
         z = z.*filt;
        
        Relerrs = [Relerrs, norm(x - exp(-1i*angle(trace(x'*z))) * z, 'fro')/norm(x,'fro')]; % Relative error
        if norm(compute_grad_super_v2(z,u,y,Params,A,At,x1),'fro')< Params.y*u
            u = Params.y1*u;
        end
    end
%     figure,imagesc(angle(z))
%     figure,imagesc(wrap(angle(z)))
    %% otro
    tau1       = Params.tau1;
    tau2       = Params.tau1;
    lmb_max1   = norm(abs(z), inf);
    lmb_max2   = norm(angle(z), inf);
    
    N1 = Params.n1;
      % magnitude
      x1m   = zeros(N1);
      tem   = abs(z(intervalo,intervalo));
     [temp2]    = admm_tv(tem(:), tau1*lmb_max1, Params.rho, [128 128],10);
      x1m(intervalo,intervalo)   =  reshape(temp2,[128 128]);
     %phase 
     x1f        =zeros(N1);
     temp       =puma_ho(angle(z(intervalo,intervalo)),.5,'verbose','no');    
     [temp2]    = admm_tv(temp(:), tau2*lmb_max2, Params.rho, [128 128],10); 
     x1f(intervalo,intervalo) = reshape(temp2,[128 128]);
    x1=x1m.*exp(1i.*x1f);
    
    Relerrs = [Relerrs, norm(x - exp(-1i*angle(trace(x'*z))) * z, 'fro')/norm(x,'fro')]; 
    z1 = reshape(x1,size(z));
    
    Params.rho  = Params.rho./1.12; 
    
%         figure,imagesc(angle(z1))
  Relerrs(end)
  colormap gray
    subplot(3,1,1),plot(Relerrs), title(j);
    z_rc = puma_ho(angle(z(intervalo,intervalo)),.5,'verbose','no');
    subplot(3,1,2),imagesc(z_rc); title(num2str(Relerrs(end)))
    subplot(3,1,3),imagesc(puma_ho(angle(z1(intervalo,intervalo)),.5,'verbose','no'));
    pause(0.01)
end

    function [ FDh, FDinv ] = totl_vart(N)
% horizontal difference operator
FDh         = zeros(N);
FDh(1,1)    = -1;
FDh(1,end)  = 1;
FDh         = fft2(FDh);
FDhH        = conj(FDh);

% vertical difference operator
FDv         = zeros(N);
FDv(1,1)    = -1;
FDv(end,1)  = 1;
FDv         = fft2(FDv);
FDvH        = conj(FDv);

FDinv       = 1./( FDhH.* FDh + FDvH.* FDv + 1); 

function [ Y ] = soft(X, Kappa)
Y          = max( 0, X - Kappa ) - max( 0, -X - Kappa );