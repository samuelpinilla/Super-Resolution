clc;
close all;
clear;

%% algorithm to design uniformly distributed coded apertures
%% variables

d     = [1i,1,-1,-1i]; % random variable
N     = 64;            % size of the coded aperture
Cd    = length(d);     % number of coding elements
L     = 4;             % number of shots
masks = zeros(N,N,L);  % coded apertures

for l=1:L
    aux = coded_design(length(d),N);
    for t=1:length(d)
        masks(:,:,l) = masks(:,:,l) + aux(:,:,t)*d(t);
        figure, imagesc(aux(:,:,t))
    end
    figure,imagesc(angle(masks(:,:,l))), title(strcat('experiment: ',num2str(l)));
end

%% save results
save('masks.mat','masks');