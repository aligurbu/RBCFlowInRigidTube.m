function [F, M, IntnormF, IntnormM] = integrateForce(axi, bxi, f)
%%
% f = the membrane forces per unit area of sphere
% (NOT per unit area of deformed configuration = f./Jxibrev;)
% Therefore, this integration scheme does not involve the Jacobian
% determinant, Jxibrev, of deformed configuration
%%
xi = shsgcm(axi,bxi);
m = cross(xi, f, 3); % cross product of vectors in DIM 3
nlat = size(xi,1);
nlon = size(xi,2);
[~,wg] = gaqdm(nlat);
F = reshape(sum(sum(f.*wg))*(2*pi/nlon),3,1);
M = reshape(sum(sum(m.*wg))*(2*pi/nlon),3,1);

normF = sqrt(f(:,:,1).^2 + f(:,:,2).^2 + f(:,:,3).^2);
normM = sqrt(m(:,:,1).^2 + m(:,:,2).^2 + m(:,:,3).^2);

IntnormF = sum(sum(normF.*wg))*(2*pi/nlon);
IntnormM = sum(sum(normM.*wg))*(2*pi/nlon);