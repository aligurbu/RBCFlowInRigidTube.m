function [V, A, Un] = getVolumeArea(axi, bxi, u)
%% Compute volume, surface area and
%% rate of volume change: \int u.n d\Gamma
%%
xi = shsgcm(axi,bxi);
[gxi_thet, gxi_phi] = gradgcm(axi, bxi);
normalVector = cross(gxi_thet, gxi_phi, 3); % cross product of vectors
Jacobian = sqrt(sum(normalVector.^2,3));
xi_dot_normalVector = sum(xi.*normalVector,3);
nlat = size(xi,1);
nlon = size(xi,2);
[~,wg] = gaqdm(nlat);
V = sum(sum(xi_dot_normalVector.*wg))*(2*pi/nlon)/3;
A = sum(sum(Jacobian.*wg))*(2*pi/nlon);

if nargin == 3
    u_dot_gamm = sum(u.*normalVector,3);
    Un = sum(sum(u_dot_gamm.*wg))*(2*pi/nlon);
end