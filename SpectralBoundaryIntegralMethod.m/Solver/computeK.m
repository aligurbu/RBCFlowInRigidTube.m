function Kernel = computeK(cxi, CHI, RotationMatrix, ...
                           eta, wg, mask_a, mask_b, ...
                           nlat, nlon, N, NGSphere, NGthet, NGphi)
%% Pre-compute the integrand with K kernel
Kernel = cell(nlat*nlon,1);

% parfor m = 1:nlat*nlon 
% parfor works and it is faster with 100% cpu usage
for m = 1:nlat*nlon
    KK = zeros(NGthet, NGphi,6);
    Normal = zeros(NGthet, NGphi,3);
    chi = reshape(CHI(:,m),1,1,3);
    Rot = RotationMatrix{m};

    %% Sample points 
    [Xihat, ahatXi_, bhatXi_] = ...
                      InvTransformUpSampleRotateSHCoeff(cxi, Rot, ...
                                                        mask_a, mask_b, ...
                                                        N, NGSphere);
    %%
    [gXihat_thet,gXihat_phi] = gradgcm(ahatXi_,bhatXi_);
    rr = chi - Xihat;
    normrr = sqrt(sum(rr.^2,3));
    rrhat = rr./normrr;

    Normal(:,:,1) = gXihat_thet(:,:,2).*gXihat_phi(:,:,3) - ...
                    gXihat_thet(:,:,3).*gXihat_phi(:,:,2);
    Normal(:,:,2) = gXihat_thet(:,:,3).*gXihat_phi(:,:,1) - ...
                    gXihat_thet(:,:,1).*gXihat_phi(:,:,3);
    Normal(:,:,3) = gXihat_thet(:,:,1).*gXihat_phi(:,:,2) - ...
                    gXihat_thet(:,:,2).*gXihat_phi(:,:,1);

    KK(:,:,1) = rrhat(:,:,1).*rrhat(:,:,1);
    KK(:,:,2) = rrhat(:,:,2).*rrhat(:,:,2);
    KK(:,:,3) = rrhat(:,:,3).*rrhat(:,:,3);
    KK(:,:,4) = rrhat(:,:,1).*rrhat(:,:,2);
    KK(:,:,5) = rrhat(:,:,2).*rrhat(:,:,3);
    KK(:,:,6) = rrhat(:,:,1).*rrhat(:,:,3);
    Kernel{m} = 6*KK.*(sum(rrhat.*Normal,3)./normrr) ...
                    .*(eta./normrr).*wg*(2*pi/NGphi);
end
end