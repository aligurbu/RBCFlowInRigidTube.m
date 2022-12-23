function Gf = NearlySingular_ComputeGf(cf, cxi, xiEqSpaced, ...
                                       EvaluationPts, EvaluationPtsLocation, ...
                                       RotationMatrix_EqSpaced, ...
                                       wg, mask_a, mask_b, mu, lam, ...
                                       N, NGSphere, NGphi)
%% Computation of integral with Gf when the target point 
%% is at exterior or interior to the cell boundary
%%
BI = zeros(size(EvaluationPts));
for n = 1 : size(EvaluationPts,2)
    chi = reshape(EvaluationPts(:,n),1,1,3);
    [~,IImin] = min(sqrt(sum((xiEqSpaced-EvaluationPts(:,n)).^2,1)));
    Rot = RotationMatrix_EqSpaced{IImin};
        
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
                  (1./normrr).*wg*(2*pi/NGphi)));
    
    if EvaluationPtsLocation(n) % The target point inside
        BI(:,n) = BI_/lam;
    else % The target point outside
        BI(:,n) = BI_;
    end
end
%%
Gf = BI(:)/(8*pi*mu);
end