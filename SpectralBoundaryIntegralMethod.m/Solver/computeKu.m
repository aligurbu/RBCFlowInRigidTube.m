function cKu = computeKu(Kernel, cu, u, RotationMatrix, ...
                         mask_a, mask_b, lam, N, nlat, nlon, NGSphere)
%% Computation of integral with K(u-u_chi) integrand
%%
BI = zeros(nlat*nlon,3);
u_chi = reshape(u,nlat*nlon,3)';
%%
for m = 1:nlat*nlon
    uchi = reshape(u_chi(:,m),1,1,3);
    Rot = RotationMatrix{m};
    %% Density
    uhat = InvTransformUpSampleRotateSHCoeff(cu, Rot, ...
                                             mask_a, mask_b, ...
                                             N, NGSphere);
    K = Kernel{m};

    delta_u = uhat - uchi;
    BI_ = [sum(sum(delta_u(:,:,1).*K(:,:,1) + ...
                   delta_u(:,:,2).*K(:,:,4) + ...
                   delta_u(:,:,3).*K(:,:,6)));
           sum(sum(delta_u(:,:,1).*K(:,:,4) + ...
                   delta_u(:,:,2).*K(:,:,2) + ...
                   delta_u(:,:,3).*K(:,:,5)));
           sum(sum(delta_u(:,:,1).*K(:,:,6) + ...
                   delta_u(:,:,2).*K(:,:,5) + ...
                   delta_u(:,:,3).*K(:,:,3)))]*((lam-1)/8/pi);
    BI(m,:)= u_chi(:,m) + BI_;
end
%%
BI = reshape(BI,nlat,nlon,3);
[aBI,bBI] = shagcm(BI);
cKu = [aBI(mask_a); bBI(mask_b)];
end