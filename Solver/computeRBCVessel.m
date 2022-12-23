function X = computeRBCVessel(X, cX_Vessel, RBCVessel_LHS, ...
                              mask_a, mask_b, N, nlat, nlon)
%%
aX = zeros(size(mask_a));   bX = zeros(size(mask_b));
aX(mask_a) = X(1:3*(N+1)*(N+2)/2);
bX(mask_b) = X(3*(N+1)*(N+2)/2+1:3*(N+1)^2);
X_ = shsgcm(aX, bX);
X__ = reshape(X_ , nlat*nlon,3)';
X = X__(:);

X = X + RBCVessel_LHS * cX_Vessel;

[aX, bX] = shagcm(reshape(reshape(X,3,nlat*nlon)',nlat,nlon,3));
X = [aX(mask_a); bX(mask_b)];