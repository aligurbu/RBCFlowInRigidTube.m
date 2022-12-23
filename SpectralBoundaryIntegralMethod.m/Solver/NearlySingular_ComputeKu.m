function Ku = NearlySingular_ComputeKu(Kernel, cu, xiEqSpaced, uEqSpaced, ...
                                       EvaluationPts, EvaluationPtsLocation, ...
                                       RotationMatrix_EqSpaced, ...
                                       mask_a, mask_b, lam, ...
                                       N, NGSphere)
%% Computation of integral with K(u-u_\perp) when the target point 
%% is at exterior or interior to the cell boundary
%%
BI = zeros(size(EvaluationPts));
for n = 1 : size(EvaluationPts,2)
    chi = EvaluationPts(:,n);
    [~,IImin] = min(sqrt(sum((xiEqSpaced-chi).^2,1)));
    Rot = RotationMatrix_EqSpaced{IImin};

    uchiperp = reshape(uEqSpaced(:,IImin),1,1,3);
    %% Density
    uhat = InvTransformUpSampleRotateSHCoeff(cu, Rot, ...
                                             mask_a, mask_b, ...
                                             N, NGSphere);
    K = Kernel{n};

    delta_u = uhat - uchiperp;
    BI_ = [sum(sum(delta_u(:,:,1).*K(:,:,1) + ...
                   delta_u(:,:,2).*K(:,:,4) + ...
                   delta_u(:,:,3).*K(:,:,6)));
           sum(sum(delta_u(:,:,1).*K(:,:,4) + ...
                   delta_u(:,:,2).*K(:,:,2) + ...
                   delta_u(:,:,3).*K(:,:,5)));
           sum(sum(delta_u(:,:,1).*K(:,:,6) + ...
                   delta_u(:,:,2).*K(:,:,5) + ...
                   delta_u(:,:,3).*K(:,:,3)))];
               
    if EvaluationPtsLocation(n) % The target point inside
       BI(:,n) = ((lam-1)/lam)*(-uEqSpaced(:,IImin) + BI_/(8*pi));
    else % The target point outside
       BI(:,n) = ((lam-1)/(8*pi))*BI_;
    end
end
Ku = BI(:);
end