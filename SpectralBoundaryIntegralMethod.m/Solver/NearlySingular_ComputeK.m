function Kernel = NearlySingular_ComputeK(cxi, xiEqSpaced, ...
                                          EvaluationPts, ...
                                          RotationMatrix_EqSpaced, ...
                                          wg, mask_a, mask_b, ...
                                          N, NGSphere, NGthet, NGphi)
%% Pre-compute the integrand with K kernel for the nearly-singular integration
Kernel = cell(size(EvaluationPts,2),1);
for m = 1:size(EvaluationPts,2)
    KK = zeros(NGthet, NGphi,6);
    Normal = zeros(NGthet, NGphi,3);
    chi = reshape(EvaluationPts(:,m),1,1,3);
    [~,IImin] = min(sqrt(sum((xiEqSpaced-EvaluationPts(:,m)).^2,1)));
    Rot = RotationMatrix_EqSpaced{IImin};

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
                    .*(1./normrr).*wg*(2*pi/NGphi);
end
end