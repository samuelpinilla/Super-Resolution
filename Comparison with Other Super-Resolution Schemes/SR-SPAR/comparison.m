close all
clear all
clc;
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
    
N_sensor=1024*4;           %% Square sensor size in pixels
Np=N_sensor/support_part;  %% FFT number of pixels
N_object=512/2;            %% Object size
           
%% for superresolution
IterNum=50;               %% ITERATION NUMBER
L = 4;                   %% Number of experiments

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
filtering=0;               %%  1 for BM3D filtering
unwrapping=0;              %%  1 -means unwrapping on before BM3D phase filtering, 0- means unwrapping off
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
 const =delta_x_0^2/lambda/f;
 
 
Masks   = zeros(yN,xN,L);        %% Storage for L modulation masks, each of dim yN x xN
m       = yN*xN*L;                %% Data length

rand('seed',10001), randn('seed',1001);



if (if_superresoltion==1)&(N~=1)&(M~=1)
    Mask_scalled=zeros(yN,xN,L);
for ll = 1:L
   for sy=1:yN/N;%% 
   for sx=1:xN/M;
          Mask_scalled(( sy-1)*N+(1:N),(sx-1)*M+(1:M),ll)=exp(1j*randn(1,1)*pi/2/2)*(1+ampl_for_masks*(rand(1,1)));
    end
    end
%   
end
    Masks=Mask_scalled;
end

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
%         cop(:,:,ww) = fftshift((fft2(x.*Masks(:,:,ww))).*delta_x_0^2/lambda/f);
 y(:,:,ww)=abs(temp).^2;                     %% True intensity array
end
    figure,imagesc(y(:,:,1))
%%%%%%%%%%%%%% MAIN PROGRAM %%%%%%%%%%%%%%
clear temp
                                             %% Noisy observations 

rand('seed',100001), randn('seed',1001);
%%%%%%%%%%%%%%%%%%%%%%%
z=y;
zz=z;
z_downsampled=z;y_downsampled=y;
if if_superresoltion==1;
[z_upsampled,z_downsampled] = Down_Up_Sampling(z,N,M);
[y_upsampled,y_downsampled] = Down_Up_Sampling(y,N,M);
z=z_upsampled;
zz=z;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end 

 clear  z_upsampled
 
 %% Inicialization parameters
if exist('Params')                == 0,  Params.n1          = Np;       end
if isfield(Params, 'n2')          == 0,  Params.n2          = Np;       end
if isfield(Params, 'Rn')          == 0,  Params.Rn          = 1;        end  % escala de la imagen
if isfield(Params, 'L')           == 0,  Params.L           = L;        end 
if isfield(Params, 'p')           == 0,  Params.p           = 2;        end
if isfield(Params, 'npower_iter') == 0,  Params.npower_iter = 200;      end
if isfield(Params, 'u0')          == 0,  Params.u0          = 10;       end
Params.m  = Np*Np*L;
 
 %% main algorithm 
    
Params.Masks = Masks;
Params.const = const;

B = @(I) fft2(Masks .* reshape(repmat(I,[1 L]), size(I,1), size(I,2), L))*const;
Bt = @(I) sum(conj(Masks) .*ifft2(I), 3) * size(I,1)* size(I,2).*(1/const);

y_2 = abs(B(x));
  
 tic
[z0,Relerrs] = Inicialization(x,y_2,Params,B,Bt);
toc
 
figure,imagesc(angle(x))
figure,imagesc(angle(z0))
