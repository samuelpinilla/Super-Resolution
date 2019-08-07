clc;
clear all;
close all;
clear clases;
% 
addpath('DATA');
addpath('PUMA');
varphi=double(imread('lena512.bmp'))/255*pi/2;
Bo=ones(size(varphi))*10;                           %% Uniform amplitude

%% Parametros
rm = 1;
Np = 1024; %512
No = 256;
itervalo = (Np/2-No/2)+1:(Np/2+No/2);


x = zeros(Np);
x(itervalo,itervalo)=imresize(Bo.*exp(1j*varphi),No/512);
[n1,n2] = size(x);

%% Inicialization parameters
if exist('Params')                == 0,  Params.n1          = n1;       end
if isfield(Params, 'n2')          == 0,  Params.n2          = n2;       end
if isfield(Params, 'Rn')          == 0,  Params.Rn          = 1;        end  % escala de la imagen
if isfield(Params, 'L')           == 0,  Params.L           = 1;        end 
if isfield(Params, 'p')           == 0,  Params.p           = 7;      end
if isfield(Params, 'npower_iter') == 0,  Params.npower_iter = 200;      end
if isfield(Params, 'u0')          == 0,  Params.u0          = 30;       end

L           = Params.L; 
Rn          = Params.Rn;
d1          = n1/Rn; 
d2          = n2/Rn;  

Params.m    = d1*d2*L;
Params.itervalo = itervalo;
%% Mask
rand('seed',1001);
for tt=1:L
Masks(:,:,tt)   = randsrc(n1,n2,[1,1i,-1,1i],round(rand(1)*100));
end
%% para
const = 1;
Params.const = const;

lambda = 0.633*10^-6; % wavelength, unit: m
delta  = 10*lambda;   % sampling period, unit: m
M      = n1;          % object size
z      = 0.07;        % propagation distance; m
c      = 1:M;
r      = 1:M;
[C, R] = meshgrid(c, r);
THOR   = ((R-M/2-1).^2+(C-M/2-1).^2).^0.5;
A      = THOR.*delta;
Q1 = exp(-1i*pi/lambda/z.*(A.^2));

A  = @(I) ifft2(reshape(repmat(I.*Q1,[1 L]), size(I,1), size(I,2), L).*conj(Masks))...
    .* size(I,1) * size(I,2);
At = @(I) sum(fft2(I).*reshape(repmat(conj(Q1),[1 L]), size(I,1), size(I,2), L).*Masks , 3);

%% Initialization
y = abs(A(x));

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
toc

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
% figure,imshow(z_rc,[]);
PSNR_mag = fun_PSNR(abs(x_r),z_rm)
PSNR_mag_2 = fun_PSNR(abs(x_r),z_rm1)

PSNR_pha = fun_PSNR(angle(x_r),z_rc)
PSNR_pha_2 = fun_PSNR(angle(x_r),z_r1)

% close all
figure(10), subplot(2,2,1), imshow(z_rc,[]), title(['RMSE_{phase} =', num2str(norm_pha)],'FontSize',14),...
    subplot(2,2,2), imshow(z_rm,[9 11]), title(['RMSE_{ampl} =', num2str(norm_mag)],'FontSize',14)
    subplot(2,2,3), plot(angle(x_r(No/2,:))),  hold on, grid on, plot(z_rc(No/2,:),'-r'),title(['PSNR =', num2str(PSNR_pha),' dB'],'FontSize',14),...,
    subplot(2,2,4),  plot(abs(x_r(No/2,:))) , hold on, grid on,plot(z_rm(No/2,:),'-r'),title(['PSNR =', num2str(PSNR_mag),' dB'],'FontSize',14) %% title(['IterNum=', num2str(0),'; RMSE-ampl=', num2str(RMSE_ampl(1),3)])

figure(11), subplot(2,2,1), imshow(z_r1,[]), title(['RMSE_{phase} =', num2str(norm_pha2)],'FontSize',14),...
    subplot(2,2,2), imshow(z_rm1,[9 11]), title(['RMSE_{ampl} =', num2str(norm_mag2)],'FontSize',14)
    subplot(2,2,3), plot(angle(x_r(No/2,:))),  hold on, grid on, plot(z_r1(No/2,:),'-r'),title(['PSNR =', num2str(PSNR_pha_2),' dB'],'FontSize',14),...,
    subplot(2,2,4),  plot(abs(x_r(No/2,:))) , hold on, grid on,plot(z_rm1(No/2,:),'-r'),title(['PSNR =', num2str(PSNR_mag),' dB'],'FontSize',14) %% title(['IterNum=', num2str(0),'; RMSE-ampl=', num2str(RMSE_ampl(1),3)])

