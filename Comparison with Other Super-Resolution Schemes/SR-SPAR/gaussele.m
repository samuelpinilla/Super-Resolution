function  z = gaussele(dimx,dimy,top,sx,sy)
%         z = gaussele(dimx,dimy,top,sx,sy);
%         Genetares a gaussian elevation of size 
%         MATLAB 5.x
%        
%   [dimx,dimy]= dim(z)
%   top   -> max(max(z))
%   sx,sy -> along x and y standard deviations
%
%   Author: J.M. Bioucas Dias, 2000
%   Topic -  Interferometry

[x y] = meshgrid(-dimy/2+1:dimy/2,-dimx/2+1:dimx/2);
z=top*(exp(-x.^2/(2*sx^2)-y.^2/(2*sy^2)));
