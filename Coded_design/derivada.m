function [imhor, imver]=derivada(gdmd)
%  load('data4shot2014sep18.mat')
%  A=gdmd(:,:,1);
 [N,N,S]=size(gdmd);
 for j=1:S
     A=gdmd(:,:,j);
     f1=imfilter(A,[1 -1]);
     f2=imfilter(A,[-1 1 0]);
     f3=imfilter(A,[1;-1]);
     f4=imfilter(A,[-1;1;0]);
     imhor(:,:,j)=(A.*(f1==0)+A.*(f2==0))>0;
     imver(:,:,j)=(A.*(f3==0)+A.*(f4==0))>0;
 end
% close all;imagesc(A(1:8,1:15));figure;imagesc(f);colorbar