function [nlat, nlon, thet, phi, weight] = GridOnSphere(N)
%% Create grid for representation
nlat = N+1;
nlon = 2*N+1;
[thet, weight] = gaqdm(nlat);
phi = 0:2*pi/nlon:(1-1/nlon)*(2*pi);