close all
clear all
addpath .\DATA
addpath .\PUMA


%% PARAMETERS of BM3D-FRAME FILTER

 Nstep = 3;  N11=8;  N22=8; 
threshType='h';                           %% hard 'h', soft 's' - thresholding
%%%%%%%%% Thresholding parameters 

threshold_ampl=1.4; % para el bm3D
threshold_phi=1.4;

N=1; M=1;  
ampl_for_masks=0;
if_superresoltion=1;     %% Switch for super-resolution
support_part=1;        %% Percentage of active pixels is equal to support_part^2
                           %% support_part=.5 gives 25% support_part=4/5  is equal to 64%
%% support_part=1/2; 
%% support_part=1; 

    %% OPTICAL SETUP PARAMETERS %%%%%%%%%%%%%
    
N_sensor=1024;           %% Square sensor size in pixels
Np=N_sensor/support_part;  %% FFT number of pixels
N_object=512/4;            %% Object size
           
%% for superresolution
IterNum=50;               %% ITERATION NUMBER
L = 2;                   %% Number of experiments

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% OPTICAL WAVEFIELD PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%
lambda = 632.8e-9;                       %% wavelength
delta_x_z = 5.2e-6;  delta_y_z = 5.2e-6; % the sensor pitch sizes
delta_x_0 = 5.2e-6; delta_y_0 = 5.2e-6;   % the object (SLM) pitch sizes

if if_superresoltion==1
N=4; M=N;       %% N is a ratio of the sensor-pixel size to the object-pixel size
                    %% R_s=N;
delta_x_z= delta_x_z/M;  delta_y_z= delta_y_z/N; % the computational sensor pitch sizes
delta_x_0 = delta_x_0/M; delta_y_0=delta_y_0/N;   % the computational object (SLM) pitch sizes

end


f=delta_x_0*delta_x_z/lambda*Np      %% focal length of lens


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CONTROL PARAMETERS of ALGORITHMS
filtering=1;               %%  1 for BM3D filtering
unwrapping=1;              %%  1 -means unwrapping on before BM3D phase filtering, 0- means unwrapping off
Fienup=1;                  %%  1 means that the Fienup rules is used instead of the optimal solution for data denoising


%% Parameter of Poisson distribution

KAPPA_set=[100000 ];
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                      
 %% Complex-valued Object Design %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Bo_Interval=[Np/2-N_object/2+1:Np/2+N_object/2];
Sensor_Interval=[Np/2-N_sensor/2+1:Np/2+N_sensor/2];
varphi=zeros(Np); 
%%%%%%%%%%%% Object phase %%%%%%%%%%%%%%%%%%%%%%%%%%%%

test_phase = 6;     %% Phase image selection
switch(test_phase)
 
   case 1
load DEMO_1                      %% Gaussian abs phase
   case 3
     load DEMO_3                      %% Truncated Gaussian abs phase 
   case 4
   case 5
   case 6
      varphi0=double(imread('image_Lena256.png'))/255*2*pi/4;
      varphi0=imresize(varphi0,0.5);
      %%   varphi0=double(imread('image_Lena512.png'))/255*2*pi/4;
          [yN,xN]=size(varphi0);
          Bo_Interval=[Np/2-N_object/2+1:Np/2+N_object/2];
 varphi=zeros(Np);

 varphi(Bo_Interval,Bo_Interval)=varphi0; 
        
 figure,imagesc(varphi)
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Pbject Amplitude 
test_ampl = 1;
switch(test_ampl)
   case 1
 
      Bo=zeros(Np); 
 
 Bo(Bo_Interval,Bo_Interval)=1.;

end
  
x=Bo.*exp(1j*varphi);                           %% TRUE COMPLEX-VALUED OBJECT
%%%%%%%%%%%%%
[yN,xN]=size(varphi);

Object_Size=delta_x_0*N_object ; %% in mm

Sensor_Size=delta_x_z*Np*support_part; %% in mm

FresnelNumber=(Object_Size/2)^2/f/lambda
resolution_in_wavelength=delta_x_0/lambda     %% R_lambda

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% OBSERVATION MODEL, MASKS and PROPAGATION %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Masks   = zeros(yN,xN,L);        %% Storage for L modulation masks, each of dim yN x xN
m       = yN*xN*L;                %% Data length
% % 
% % % rand('seed',10001), randn('seed',1001);
rand('seed',1001);


if (if_superresoltion==1)&(N~=1)&(M~=1)
    Mask_scalled=zeros(yN,xN,L);
for ll = 1:L
   for sy=1:yN/N;%% 
   for sx=1:xN/M;
%           Mask_scalled(( sy-1)*N+(1:N),(sx-1)*M+(1:M),ll)=exp(1j*randn(1,1)*pi/2/2)*(1+0*(rand(1,1)));
            Mask_scalled(( sy-1)*N+(1:N),(sx-1)*M+(1:M),ll)=randsrc(1,1,[1,1i,-1,1i],round(rand(1)*1000));
    end
    end
%   
end
    Masks=Mask_scalled;
end
% load('Mask2.mat')
clear  Mask_scalled

%% If there is no phase modulation 
%% Masks=ones(yN,xN,L);

%% OBSERVATION MODELING
%%%%%%%%%%%%%% SENSOR SUPPORT%%%%%%%%%%%
    Masks_Support=zeros(Np,Np,L);                               %% 

support_y=[Np/2-floor(Np/2*support_part)+1:Np/2+floor(Np/2*support_part)];
support_x=support_y;


for ll = 1:L, Masks_Support(support_y,support_x,ll) = 1; 
   
end

figure,imagesc(angle(Masks(:,:,1)))

figure,imagesc(abs(Masks_Support(:,:,1)))

for ww=1:L                              %% True intensities at the sensor plane
    

 temp=Masks_Support(:,:,ll).*fftshift((fft2(x.*Masks(:,:,ww)))*delta_x_0^2/lambda/f).*Masks_Support(:,:,ww);

 y(:,:,ww)=abs(temp).^2;                     %% True intensity array
end
    figure,imagesc(y(:,:,1))
%%%%%%%%%%%%%% MAIN PROGRAM %%%%%%%%%%%%%%
clear temp
s_KAPPA=0;
for KAPPA=KAPPA_set                          %% Loop for KAPPA values
                                             %% Noisy observations 
s_KAPPA=s_KAPPA+1
gamma_1=1/KAPPA;                             %% Parameter of the iterative algorithm
rand('seed',100001), randn('seed',1001);
%%%%%%%%%%%%%%%%%%%%%%%
z=poissrnd(KAPPA*y);
zz=z;
z_downsampled=z;y_downsampled=y;
if if_superresoltion==1;
% z=zeros(yN,xN,L);
% z_down_sampl=zeros(yN/N,xN/M,L);              %% Low resolution sensor
 
[z_upsampled,z_downsampled] = Down_Up_Sampling(z,N,M);
[y_upsampled,y_downsampled] = Down_Up_Sampling(y,N,M);
z=z_upsampled;
zz=z;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end 

 clear  z_upsampled
 
%%z  = poissrnd(y*KAPPA);                      %% Noisy Poissonain Intensity Observations

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
Number_Photons_per_Pixel(s_KAPPA)=sum(z_downsampled(:))/N_sensor/N_sensor/L;                 %% Mean umber of Photons per Pixel
SNR(s_KAPPA)=10*log10(sum(y_downsampled(:).^2)/sum((y_downsampled(:)-z_downsampled(:)/KAPPA).^2)) %% SNR for Poisson Data
% SNR(s_KAPPA)=10*log10(sum(y(:).^2)/sum((y(:)-z(:)/KAPPA).^2)) %% SNR for Poisson Data



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% SPAR algorithm %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic                 %% for SPAR
%% INITIALIZATION FOR SPAR

if 1
V=ones(size(x)).* exp(1j*pi*randn(size(x))*.0)*1.3;  %% RANDOM PHASE INITIALIZATION

else                                                  %% Initialization by sensor backpropagation
temp=zeros(yN,xN,L); 

for ww=1:L          %% backward propagation to the
    
    temp(:,:,ww)=(ifft2(fftshift(sqrt(z(:,:,ww))/sqrt(KAPPA))))./Masks(:,:,ww)*lambda*f/delta_x_0^2;
  
 end

 
V=sum(temp,3)/L;

end
%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

V_init=V; 


%% Accuracy of initialization
temp=abs(x)-abs(V_init);
RMSE_ampl(s_KAPPA,1)=norm(temp(Bo_Interval,Bo_Interval),'fro')/(length(Bo_Interval));
clear temp

%% Phase unwrapping
if unwrapping==1
potential.quantized = 'no';
potential.threshold = pi;
verbose='no';
phi_Uo_init=angle(V_init);
phi_Uo_init0=puma_ho(angle(V_init(Bo_Interval,Bo_Interval)),.5,'verbose',verbose); %% for demo unwrapping from raw data
phi_Uo_init(Bo_Interval,Bo_Interval)=phi_Uo_init0;

bb_init=mean(mean(varphi(Bo_Interval,Bo_Interval) - phi_Uo_init(Bo_Interval,Bo_Interval)));

temp=varphi - phi_Uo_init-bb_init;
RMSE_phi(s_KAPPA,1)=norm(temp(Bo_Interval,Bo_Interval),'fro')/(length(Bo_Interval));
clear temp bb_init

else
bb=mean(mean(wrap(varphi(Bo_Interval,Bo_Interval)) - angle(V_init(Bo_Interval,Bo_Interval)))); 
temp=wrap(angle(x) - angle(V)-bb);

RMSE_phi(s_KAPPA,1)=norm(temp(Bo_Interval,Bo_Interval),'fro')/(length(Bo_Interval));
clear temp
end

%% RMSE of the phase and amplitude estimates for use in BM3D filtered estimates
sigma_phi(1)=function_stdEst2D(angle(V(Bo_Interval,Bo_Interval)),2);
sigma_ampl(1)=function_stdEst2D(abs(V(Bo_Interval,Bo_Interval)),2);

%% MAIN ITERATIONS of SPAR PHASE RETRIEVAL
gamma_10=gamma_1;
 
%% delta_diff_2=1; delta_diff_1=0;
for ss=1:IterNum
  ss;  
  if (ss>3)&(mod(ss,2)==1)
ss0=ss;       %% ss0 is iteration when BM3D grouping is updated
  else
      ss0=ss;
  end
%%V_=V;
gamma_1=gamma_10/(1+.2*(ss-1));
%% STEP 1:     %% Forward propagation 

for ww=1:L          %% forward propagation to the sensor plane
        
    temp=fft2(V.*Masks(:,:,ww))*delta_x_0^2/lambda/f;
  
    Us(:,:,ww)=fftshift(temp);                     %% true wavefront at the sesnor plane
end
clear temp

%% STEP 2-0 Upsampling
% lambda=.2/(1+(ss-1)*0.1);
% for ll=1:L
%    
%     
% temp=imresize(abs(Us(:,:,ll)).^2,1/N,'nearest')-z_down_sampl(:,:,ll)/KAPPA;
%    
% y_var(:,:,ll)=max(y_var(:,:,ll)-lambda*(imresize(temp,N,'bicubic')),0);
% 
% 
% end
% keyboard


%% STEP 2:     %% Poissonian Filtering
if Fienup==1   %% No Filtering at Sensor Plane
  
 Us=sqrt(zz).*Masks_Support/sqrt(KAPPA).*exp(1j*angle(Us))+Us.*(1-Masks_Support);
    
else
    
    [Us]=NoiseFiltering_3D_propag(Us,zz,KAPPA,gamma_1,L).*Masks_Support+Us.*(1-Masks_Support); %% Filtering at Sensor Plane
end

%% STEP 3: Backward propagation
temp=zeros(yN,xN,L); 

for ww=1:L          %% backward propagation to the
    
   temp(:,:,ww)=(ifft2(fftshift(Us(:,:,ww))))./Masks(:,:,ww)*lambda*f/delta_x_0^2;
   
   
 end

 
Uo=mean(temp,3);

clear temp

%% STEP 4: Phase unwrapping ;
if unwrapping==1
potential.quantized = 'no';
potential.threshold = pi;
verbose='no';

phi_Uo=angle(Uo);
phi_Uo_0=puma_ho(angle(Uo(Bo_Interval,Bo_Interval)),.5,'verbose',verbose); %% for demo unwrapping from raw data
phi_Uo(Bo_Interval,Bo_Interval)=phi_Uo_0;
% 
else
    phi_Uo=double(angle(Uo));
    
end   
    
%% STEP 5: BM3D phase and amplitude filtering
%% BM3D filtering of phase
   
 sigma_phi(ss+1)=function_stdEst2D(phi_Uo(Bo_Interval,Bo_Interval),2);  %% STD of phase
 sigma_ampl(ss+1)=function_stdEst2D(abs(Uo(Bo_Interval,Bo_Interval)),2);
 

if filtering==1
   abs_Uo=abs(Uo);
phi_u0_SPAR=BM3D_SPAR_UNWRAP_PHIp(phi_Uo(Bo_Interval,Bo_Interval),threshType, sigma_phi(ss+1)*threshold_phi, N11, N22, Nstep,filtering,ss,ss0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

B_u0_SPAR=BM3D_SPAR_ABSp(abs_Uo(Bo_Interval,Bo_Interval),threshType, sigma_ampl(ss+1)*threshold_ampl, N11, N22, Nstep,filtering,ss,ss0);
%%%%%%%%%%%
%%%%%%%% new
phi_Uo(Bo_Interval,Bo_Interval)=phi_u0_SPAR;
phi_u0_SPAR=phi_Uo;
abs_Uo(Bo_Interval,Bo_Interval)=B_u0_SPAR;
B_u0_SPAR=abs_Uo;

else
phi_u0_SPAR=phi_Uo;
B_u0_SPAR=abs(Uo);


end

 %% STEP 6: UPDATE of x
V=B_u0_SPAR.*exp(1j*phi_u0_SPAR).*Bo;

% 
% if isOdd(ss)==0
%   delta_diff_1=mean(mean(abs(V(Bo_Interval,Bo_Interval)-V_(Bo_Interval,Bo_Interval))));
%   
% else 
%     delta_diff_2=mean(mean(abs(V(Bo_Interval,Bo_Interval)-V_(Bo_Interval,Bo_Interval))));
%     
% end

if unwrapping==1
    
bb=mean(mean(varphi(Bo_Interval,Bo_Interval) - phi_u0_SPAR(Bo_Interval,Bo_Interval)));                                       %% absolute phase shift


RMSE_phi(s_KAPPA,ss+1) = norm(varphi(Bo_Interval,Bo_Interval) - phi_u0_SPAR(Bo_Interval,Bo_Interval)-bb, 'fro')/sqrt(yN*xN) %% RMSE absolute phase error

V_cor=V.*exp(1j*bb);
else

bb=mean(mean(wrap(varphi(Bo_Interval,Bo_Interval)) - angle(V(Bo_Interval,Bo_Interval)))); 

temp=wrap(varphi - angle(V)-bb);
V_cor=V.*exp(1j*bb);
RMSE_phi(s_KAPPA,ss+1)=norm(temp(Bo_Interval,Bo_Interval),'fro')/(length(Bo_Interval))
clear temp
end

temp=abs(x)-abs(V);
RMSE_ampl(s_KAPPA,ss+1)=norm(temp(Bo_Interval,Bo_Interval),'fro')/(length(Bo_Interval));
clear temp



if 0   %% ITERATION SHOW
if unwrapping==1
    pause(.1)
figure(11), subplot(2,2,1), imshow(phi_Uo(Bo_Interval,Bo_Interval)+bb,[]), title(['IterNum=', num2str(ss+1),'; RMSE-abs.phase=', num2str(RMSE_phi(s_KAPPA,ss+1),3)]),...
    subplot(2,2,2), imshow(abs(V(Bo_Interval,Bo_Interval)),[]), title(['IterNum=', num2str(ss+1),'; RMSE-ampl=', num2str(RMSE_ampl(s_KAPPA,ss+1),3)])
    subplot(2,2,3), imshow(phi_Uo_init(Bo_Interval,Bo_Interval)+bb_init,[]), title(['IterNum=', num2str(0),'; RMSE-abs.phase=', num2str(RMSE_phi(s_KAPPA,1),3)]),...
    subplot(2,2,4), imshow(abs(V_init(Bo_Interval,Bo_Interval)),[]), title(['IterNum=', num2str(0),'; RMSE-ampl=', num2str(RMSE_ampl(s_KAPPA,1),3)])
else
        pause(2)
 figure(44), subplot(2,2,1), imshow(angle(V_cor(Bo_Interval,Bo_Interval)),[]), title(['IterNum=', num2str(ss+1),'; RMSE-abs.phase=', num2str(RMSE_phi(s_KAPPA,ss+1),3)]),...
    subplot(2,2,2), imshow(abs(V_cor(Bo_Interval,Bo_Interval)),[]), title(['IterNum=', num2str(ss+1),'; RMSE-ampl=', num2str(RMSE_ampl(s_KAPPA,ss+1),3)])
    


subplot(2,2,3), imshow(angle(V_cor_init(Bo_Interval,Bo_Interval)),[]), title(['IterNum=', num2str(0),'; RMSE-abs.phase=', num2str(RMSE_phi(s_KAPPA,1),3)]),...
    subplot(2,2,4), imshow(abs(V_cor_init(Bo_Interval,Bo_Interval)),[]), title(['IterNum=', num2str(0),'; RMSE-ampl=', num2str(RMSE_ampl(s_KAPPA,1),3)])
    end
end
%%[ss RMSE_phi(s_KAPPA,ss) RMSE_ampl(s_KAPPA,ss)]

figure(22), subplot(1,2,1), plot(RMSE_phi(:,1:ss+1)','-o'), grid on, title('RMSE phase '), xlabel('ITERATIONS'),...,
    subplot(1,2,2), plot(RMSE_ampl(:,1:ss+1)','-o'), grid on,title('RMSE ampl'),xlabel('ITERATIONS')
pause(.1)

%%%%%%%%%%%% STOPPING RULE %%%%%%%%%%%
% if 0
% if (ss>10)&(delta_diff_1>1.5*delta_diff_2)
%    
%     'break'
%     break
%     
% end
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
SPAR_time(s_KAPPA)=toc    %% for SPAR

if 0
if unwrapping==1
figure(3), subplot(2,2,1), imshow(phi_Uo+bb,[]), title(['IterNum=', num2str(ss+1),'; RMSE-abs.phase=', num2str(RMSE_phi(ss+1),3)]),...
    subplot(2,2,2), imshow(abs(V),[]), title(['IterNum=', num2str(ss+1),'; RMSE-ampl=', num2str(RMSE_ampl(ss+1),3)])
    subplot(2,2,3), imshow(phi_Uo_init+bb_init,[]), title(['IterNum=', num2str(0),'; RMSE-abs.phase=', num2str(RMSE_phi(1),3)]),...
    subplot(2,2,4), imshow(abs(V_init),[]), title(['IterNum=', num2str(0),'; RMSE-ampl=', num2str(RMSE_ampl(1),3)])
    else
 figure(3), subplot(2,2,1), imshow(angle(V_cor(Bo_Interval)),[]), title(['IterNum=', num2str(ss+1),'; RMSE-abs.phase=', num2str(RMSE_phi(ss+1),3)]),...
    subplot(2,2,2), imshow(abs(V_cor(Bo_Interval)),[]), title(['IterNum=', num2str(ss+1),'; RMSE-ampl=', num2str(RMSE_ampl(ss+1),3)])
    subplot(2,2,3), imshow(angle(V_cor_init(Bo_Interval)),[]), title(['IterNum=', num2str(0),'; RMSE-abs.phase=', num2str(RMSE_phi(1),3)]),...
    subplot(2,2,4), imshow(abs(V_cor_init(Bo_Interval)),[]), title(['IterNum=', num2str(0),'; RMSE-ampl=', num2str(RMSE_ampl(1),3)])
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% Images for results %%%%%%

%% Metricas

x_r = x(Bo_Interval,Bo_Interval);
V_cor = exp(-1i*angle(trace(x'*V_cor))) * V_cor;
z_rc = angle(V_cor(Bo_Interval,Bo_Interval));
z_rm = abs(V_cor(Bo_Interval,Bo_Interval));


norm_pha = norm(angle(x_r)-z_rc,'fro')./norm(angle(x_r));
norm_mag = norm(abs(x_r)-z_rm,'fro')./norm(abs(x_r));

PSNR_pha = fun_PSNR(angle(x_r),z_rc)

%%PLOTS_Combined
figure(s_KAPPA), subplot(2,2,1), imshow(angle(V_cor(Bo_Interval,Bo_Interval)),[]), title(['RMSE_{phase} =', num2str(norm_pha)],'FontSize',14),...
    subplot(2,2,2), imshow(abs(V_cor(Bo_Interval,Bo_Interval)),[0.5 1.5]), title(['RMSE_{ampl} =', num2str(norm_mag)],'FontSize',14)
    subplot(2,2,3), plot(angle(x(Np/2,Bo_Interval))),  hold on, grid on, plot(angle(V_cor(Np/2,Bo_Interval)),'-r'),title(['PSNR =', num2str(PSNR_pha),', R_{\lambda} = ', num2str(resolution_in_wavelength,3),' \lambda'],'FontSize',14),...,
    subplot(2,2,4),  plot(abs(x(Np/2,Bo_Interval))) , hold on, grid on,plot(abs(V_cor(Np/2,Bo_Interval)),'-r'),title(['\chi =', num2str(KAPPA,'%10.1e\n'),', SNR =', num2str(SNR(s_KAPPA),3),' dB'],'FontSize',14) %% title(['IterNum=', num2str(0),'; RMSE-ampl=', num2str(RMSE_ampl(1),3)])
pause(.1)

%%%%%%%%%%%%%%%%%% mesh for phase
% if unwrapping==0
% % potential.quantized = 'no';
% % potential.threshold = pi;
% % verbose='no';
% % %% phi_Uo=double(puma_ho(angle(Uo),.5,'verbose',verbose)); %% for demo unwrapping from raw data
% % 
% % phi_Uo=phi_u0_SPAR;
% % phi_Uo_0=puma_ho(phi_Uo(Bo_Interval,Bo_Interval),.5,'verbose',verbose); %% for demo unwrapping from raw data
% % phi_Uo(Bo_Interval,Bo_Interval)=phi_Uo_0;
% % % 
% % % 
% % bb=mean(mean(varphi(Bo_Interval,Bo_Interval) - phi_Uo_0));
% % % 
% % RMSE_phi11(s_KAPPA) = norm(varphi(Bo_Interval,Bo_Interval) - phi_Uo_0-bb, 'fro')/sqrt(yN*xN); %% RMSE Initial rel. phase error
% % figure(s_KAPPA+30), surfl(phi_Uo_0-bb);shading interp;colormap(gray),...
% %     title(['ABS.PHASE, SNR = ', num2str(SNR(s_KAPPA),3),' dB',', RMSE_{phase} =', num2str(RMSE_phi11(s_KAPPA),3)],'FontSize',14)
% % %% mesh(phi_u0_SPAR(Bo_Interval,Bo_Interval)), 
% else
 
%     figure(s_KAPPA+30), surfl(phi_u0_SPAR(Bo_Interval,Bo_Interval)-bb);shading interp;colormap(gray);...
%     title(['ABS.PHASE, SNR = ', num2str(SNR(s_KAPPA),3),' dB',', RMSE_{phase} =', num2str(RMSE_phi(s_KAPPA,end),3)],'FontSize',14)


% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end  %% KAPPA_set

Time_total=toc









