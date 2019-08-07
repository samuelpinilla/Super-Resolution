function [DMD] = DMD_Boolean(M,N, shots)
%Generates a MxN Boolean coded aperture set for a specific number of shots.
%This set has just one one-valued entry on each spatial position per shot.
%Returns the 3-dimensional structure of codes

DMD = zeros(M*N, shots);
for j=1:N^2
    DMD(j, ceil(shots*rand)) = 1;
end
DMD = reshape(DMD, [M N shots]);
