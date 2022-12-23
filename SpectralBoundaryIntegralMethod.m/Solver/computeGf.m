function cGf = computeGf(cf, cxi, CHI, RotationMatrix, ...
                         eta, wg, mask_a, mask_b, mu, ...
                         nlat, nlon, N, NGSphere, NGphi)
%% Computation of integral with Gf
%%
BI = zeros(nlat*nlon,3);
% parfor m = 1:nlat*nlon
% parfor works but just little bit faster with a lot of cpu usage
for m = 1:nlat*nlon
    chi = reshape(CHI(:,m),1,1,3);
    Rot = RotationMatrix{m};
    %% Sample points
    xihat = InvTransformUpSampleRotateSHCoeff(cxi, Rot, ...
                                              mask_a, mask_b, ...
                                              N, NGSphere);
    %% Density
    fhat = InvTransformUpSampleRotateSHCoeff(cf, Rot, ...
                                             mask_a, mask_b, ...
                                             N, NGSphere);
    %%
    rr = chi - xihat;
    normrr = sqrt(sum(rr.^2,3));
    rrhat = rr./normrr;
    BI_ = sum(sum((fhat + rrhat.*sum((rrhat.*fhat),3)).* ...
                  (eta./normrr).*wg*(2*pi/NGphi)));
    BI(m,:) = reshape(BI_,1,3);
end
%%
BI = reshape(BI,nlat,nlon,3);
[aBI,bBI] = shagcm(BI);
cGf = [aBI(mask_a); bBI(mask_b)]/(8*pi*mu);
end