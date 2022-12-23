function [nlat, nlon, thet, phi] = EquallySpacedGridOnSphere(N)
%% Create equally-spaced grid for the nearest points
nlat = N+1;
nlon = 2*N+1;
thet = linspace(0,pi,nlat)';
phi = 0:2*pi/nlon:(1-1/nlon)*(2*pi);