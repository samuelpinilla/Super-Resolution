clc;
clear all;
close all;
clear clases;
% 
addpath('DATA');
addpath('PUMA');
varphi=double(imread('lena512.bmp'))/255*pi/2;
varphi = imresize(varphi,0.5);
Bo=ones(size(varphi))*1;                           %% Uniform amplitude

%% Parametros
rm = 4;  % super_resolution factor
Np = 1024; %
No = 128;
itervalo = (Np/2-No/2)+1:(Np/2+No/2);


x = zeros(Np);
x(itervalo,itervalo)=imresize(Bo.*exp(1j*varphi),No/256);
[n1,n2] = size(x);

%% Inicialization parameters
if exist('Params')                == 0,  Params.n1          = n1;       end
if isfield(Params, 'n2')          == 0,  Params.n2          = n2;       end
if isfield(Params, 'Rn')          == 0,  Params.Rn         = 1;        end  % escala de la imagen
if isfield(Params, 'L')           == 0,  Params.L           = 2;        end 
if isfield(Params, 'p')           == 0,  Params.p           = 7;      end
if isfield(Params, 'npower_iter') == 0,  Params.npower_iter = 200;      end
if isfield(Params, 'u0')          == 0,  Params.u0          = 30;       end

L           = Params.L; % 
Rn          = Params.Rn; % 
d1          = n1/Rn; % 
d2          = n2/Rn;  % 
Params.m    = d1*d2*L;
Params.itervalo = itervalo;
%% Mask
rand('seed',1001);
Masks = zeros(n1,n2,L);
for tt=1:L
    for sy=1:n1/rm;%%
        for sx=1:n2/rm;
%             Masks(( sy-1)*rm+(1:rm),(sx-1)*rm+(1:rm),tt)=exp(1j*randn(1,1)*pi/2/2)*(1+0*(rand(1,1)));
            Masks(( sy-1)*rm+(1:rm),(sx-1)*rm+(1:rm),tt)=randsrc(1,1,[1,1i,-1,1i],round(rand(1)*1000));
        end
    end
end

%% para
lambda = 632.8e-9;                       %% wavelength
delta_x_z = 5.2e-6/rm;  delta_y_z = 5.2e-6/rm; % the sensor pitch sizes
delta_x_0 = 5.2e-6/rm; delta_y_0 = 5.2e-6/rm;   % the object (SLM) pitch sizes
f=delta_x_0*delta_x_z/lambda*Np;  
const = delta_x_0^2/lambda/f;

Params.const = const;
Params.u0          = 10*const;
%% Operators

A = @(I) fft2(conj(Masks) .* reshape(repmat(I,[1 L]), size(I,1), size(I,2), L)).*const;
At = @(I) sum(Masks .*ifft2((I)), 3) * d1 * d2./const;


%% Initialization
y = abs(A(x));

[z_upsampled,z_downsampled] = Down_Up_Sampling(y.^2,rm,rm);
y = sqrt(z_upsampled);
tic
[z0,Relerrs] = Inicialization(x,y,Params, A, At);
toc
figure,imshow(puma_ho(angle(z0(itervalo,itervalo)),.5,'verbose','no'),[])
pause(0.01)
figure,
%% Methods
if isfield(Params, 'mu')          == 0,  Params.mu          = 0.6;      end
if isfield(Params, 'y1')          == 0,  Params.y1          = 0.5;      end 
if isfield(Params, 'y')           == 0,  Params.y           = 0.9;      end

% v_2
if isfield(Params, 'rho')         == 0,  Params.rho         = 0.8;      end
if isfield(Params, 'tau1')         == 0,  Params.tau1      = 0.01;    end   

%% --------------------------------------------
% Main algorithm
%% ----------------------------------------------
tic
[z,z1,RelerrsT] = PRSF_Super_V2(x,y, Params, A, At,z0)  ;

z = exp(-1i*angle(trace(x'*z))) * z;
z1 = exp(-1i*angle(trace(x'*z1))) * z1;
z_r1 = puma_ho(angle(z1(itervalo,itervalo)),.5,'verbose','no');
z_rm1    = abs(z1(itervalo,itervalo));
toc
z_rc = puma_ho(angle(z(itervalo,itervalo)),.5,'verbose','no');
z_rm    = abs(z(itervalo,itervalo));
x_r   = x(itervalo,itervalo);

%% Metricas
norm_pha = norm(angle(x_r)-z_rc,'fro')./norm(angle(x_r));
norm_pha2 = norm(angle(x_r)-z_r1,'fro')./norm(angle(x_r));

norm_mag = norm(abs(x_r)-z_rm,'fro')./norm(abs(x_r));
norm_mag2 = norm(abs(x_r)-z_rm1,'fro')./norm(abs(x_r));

PSNR_mag = fun_PSNR(abs(x_r),z_rm)
PSNR_mag_2 = fun_PSNR(abs(x_r),z_rm1)

PSNR_pha = fun_PSNR(angle(x_r),z_rc)
PSNR_pha_2 = fun_PSNR(angle(x_r),z_r1)

% close all
figure(10), subplot(2,2,1), imshow(z_rc,[]), title(['Relative error_{phase} =', num2str(norm_pha)],'FontSize',14),...
    subplot(2,2,2), imshow(z_rm,[.5 1.5]), title(['Relative error_{Magnitude} =', num2str(norm_mag)],'FontSize',14)
    subplot(2,2,3), plot(angle(x_r(No/2,:))),  hold on, grid on, plot(z_rc(No/2,:),'-r'),title(['PSNR =', num2str(PSNR_pha),' dB'],'FontSize',14),...,
    subplot(2,2,4),  plot(abs(x_r(No/2,:))) , hold on, grid on,plot(z_rm(No/2,:),'-r'),title(['PSNR =', num2str(PSNR_mag),' dB'],'FontSize',14) %% title(['IterNum=', num2str(0),'; RMSE-ampl=', num2str(RMSE_ampl(1),3)])

figure(11), subplot(2,2,1), imshow(z_r1,[]), title(['Relative error_{phase} =', num2str(norm_pha2)],'FontSize',14),...
    subplot(2,2,2), imshow(z_rm1,[.5 1.5]), title(['Relative error_{Magnitude} =', num2str(norm_mag2)],'FontSize',14)
    subplot(2,2,3), plot(angle(x_r(No/2,:))),  hold on, grid on, plot(z_r1(No/2,:),'-r'),title(['PSNR =', num2str(PSNR_pha_2),' dB'],'FontSize',14),...,
    subplot(2,2,4),  plot(abs(x_r(No/2,:))) , hold on, grid on,plot(z_rm1(No/2,:),'-r'),title(['PSNR =', num2str(PSNR_mag),' dB'],'FontSize',14) %% title(['IterNum=', num2str(0),'; RMSE-ampl=', num2str(RMSE_ampl(1),3)])

