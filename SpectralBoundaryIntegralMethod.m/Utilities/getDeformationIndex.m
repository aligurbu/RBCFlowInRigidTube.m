function [dA, dB, dC, D, OriVec, Vs] = getDeformationIndex(axi,bxi)
%% Computing the inertia tensor of the RBC shape
%% to compute the deformation index and the orientation

xi = shsgcm(axi,bxi);
[gxi_thet, gxi_phi] = gradgcm(axi, bxi);
normalVector = cross(gxi_thet, gxi_phi, 3);
xi_dot_normalVector = sum(xi.*normalVector,3);
xi_xi_dot_normalVector = xi.*xi_dot_normalVector;
nlat = size(xi,1);
nlon = size(xi,2);
[~,wg] = gaqdm(nlat);
V = sum(sum(xi_dot_normalVector.*wg))*(2*pi/nlon)/3;
Int_xi_xi_dot_normalVector_dA = ...
                        sum(sum(xi_xi_dot_normalVector.*wg))*(2*pi/nlon)/4;
CenterMass = reshape(Int_xi_xi_dot_normalVector_dA./V,3,1);
r(:,:,1) = xi(:,:,1) - CenterMass(1);
r(:,:,2) = xi(:,:,2) - CenterMass(2);
r(:,:,3) = xi(:,:,3) - CenterMass(3);
TT(:,:,1) = sum(r.*r,3) - r(:,:,1).*r(:,:,1);
TT(:,:,2) = sum(r.*r,3) - r(:,:,2).*r(:,:,2);
TT(:,:,3) = sum(r.*r,3) - r(:,:,3).*r(:,:,3);
TT(:,:,4) = - r(:,:,1).*r(:,:,2);
TT(:,:,5) = - r(:,:,2).*r(:,:,3);
TT(:,:,6) = - r(:,:,1).*r(:,:,3);
TT = TT.*sum(r.*normalVector,3);
Inertia = (1/5)*sum(sum(TT.*wg))*(2*pi/nlon);
InertiaTensor = Inertia([1 4 6; 4 2 5; 6 5 3]);
[EigVec,D] = eig(InertiaTensor);
[d,ind] = sort(diag(D));
Ds = D(ind,ind);
Vs = EigVec(:,ind);
T = diag(Ds);
%% Three axes of ellipsoid with mass and volume of RBC
dA = 2*sqrt((5*(T(2)+T(3)-T(1)))/(2*V));
dB = 2*sqrt((5*(T(3)+T(1)-T(2)))/(2*V));
dC = 2*sqrt((5*(T(1)+T(2)-T(3)))/(2*V));

%% Deformation parameter
D = (dA-dB)/(dA+dB);

%% Orientation vector
OriVec = Vs(:,3);
