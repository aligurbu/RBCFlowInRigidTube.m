function [KN, K] = RegularIntegrals_KN_K_M(chi,xie,nhat,wJ, ...
                                           BasisFn,numDofPerNode, ...
                                           numDofPerElem)
%% Computes the regular (non-singular) element integrals
%%
r = chi - xie;
normr = sqrt(r(1,:).^2 + r(2,:).^2 + r(3,:).^2);
rhat = r./normr;

rrhat = [rhat(1,:).*rhat(1,:); rhat(1,:).*rhat(2,:); ...
         rhat(1,:).*rhat(3,:); rhat(2,:).*rhat(2,:); ...
         rhat(2,:).*rhat(3,:); rhat(3,:).*rhat(3,:)];

KK_ = rrhat.*sum(rhat.*nhat,1).*wJ.*(3./(4*pi*normr.*normr));

KN_ = KK_*BasisFn';

KN = reshape(KN_([1 2 3 2 4 5 3 5 6],:), numDofPerNode, numDofPerElem);

K__ = sum(KN_,2);
K_ = K__([1 2 3 2 4 5 3 5 6]);
K = reshape(K_,numDofPerNode,numDofPerNode);