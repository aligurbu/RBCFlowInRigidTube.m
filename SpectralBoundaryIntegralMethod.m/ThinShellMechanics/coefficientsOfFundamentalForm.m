function [gXi_thet, gXi_phi, EE, FF, GG, WW, JXibrev, ...
          unitNormalXi, gUnitNormalXi_thet, gUnitNormalXi_phi, ...
          LL, MM, NN] = ...
                    coefficientsOfFundamentalForm(aXi, bXi, UpSampleFactor)
%% Compute the first and second fundamental form coefficients
%%
N = size(aXi,1) - 1; % degree of spherical harmonic expansion
Nup = UpSampleFactor*N;

%% Up-sample
aXiup = UpSampling(aXi,N,Nup);
bXiup = UpSampling(bXi,N,Nup);

%%
[gXi_thet,gXi_phi] = gradgcm(aXiup,bXiup);

%% The first fundamental form coefficients
EE = sum(gXi_thet.*gXi_thet,3);
FF = sum(gXi_thet.*gXi_phi,3);
GG = sum(gXi_phi.*gXi_phi,3);
WW = EE.*GG - FF.^2;
JXibrev = sqrt(WW); % Jacobian determinant of undeformed surface

%% The unit normal to the undeformed surface
unitNormalXi(:,:,1) = (gXi_thet(:,:,2).*gXi_phi(:,:,3) - ...
                       gXi_thet(:,:,3).*gXi_phi(:,:,2))./JXibrev;
unitNormalXi(:,:,2) = (gXi_thet(:,:,3).*gXi_phi(:,:,1) - ...
                       gXi_thet(:,:,1).*gXi_phi(:,:,3))./JXibrev;
unitNormalXi(:,:,3) = (gXi_thet(:,:,1).*gXi_phi(:,:,2) - ...
                       gXi_thet(:,:,2).*gXi_phi(:,:,1))./JXibrev;
                   
%% Gradient of unit normal vectors 
% Note that there are many paper talks about low accuracy of 
% taking derivative of an unit normal vector in this way, the proposed 
% solution is to use the recurrence relation for the second-derivative 
% instead of taking twice a first-derivative. 
[aUnitNormalXi,bUnitNormalXi] = shagcm(unitNormalXi);
[gUnitNormalXi_thet,gUnitNormalXi_phi] = gradgcm(aUnitNormalXi,bUnitNormalXi);

%% The second fundamental form coefficients
LL = sum(gXi_thet.* gUnitNormalXi_thet, 3);
MM = sum(gXi_thet.* gUnitNormalXi_phi, 3);
NN = sum(gXi_phi .* gUnitNormalXi_phi, 3);