function [hatXup, ahatXup, bhatXup] = ...
             InvTransformUpSampleRotateSHCoeff(cX, Rot, mask_a, mask_b, ...
                                               N, NGSphere)
%% Rotate spherical harmonic coefficients 
chatX = Rot*cX;
%% Upsample 
ahatX = zeros(size(mask_a));
bhatX = zeros(size(mask_b));
ahatX(mask_a) = chatX(1:3*(N+1)*(N+2)/2);
bhatX(mask_b) = chatX(3*(N+1)*(N+2)/2+1:3*(N+1)^2);
ahatXup = UpSampling(ahatX,N,NGSphere);
bhatXup = UpSampling(bhatX,N,NGSphere);
%% Inverse transform 
hatXup = shsgcm(ahatXup,bhatXup);