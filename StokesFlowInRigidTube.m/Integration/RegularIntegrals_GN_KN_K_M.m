function [GN, KN, K] = RegularIntegrals_GN_KN_K_M(chi,xie,nhat,wJ, ...
                                                  BasisFn, mu,...
                                                  numDofPerNode, ...
                                                  numDofPerElem)
%% Computes the regular (non-singular) element integrals
%%
r = chi - xie;
normr = sqrt(r(1,:).^2 + r(2,:).^2 + r(3,:).^2);
rhat = r./normr;

rrhat = [rhat(1,:).*rhat(1,:); rhat(1,:).*rhat(2,:); ...
         rhat(1,:).*rhat(3,:); rhat(2,:).*rhat(2,:); ...
         rhat(2,:).*rhat(3,:); rhat(3,:).*rhat(3,:)];

GG_ = ([1;0;0;1;0;1] + rrhat).*wJ./(8*pi*mu*normr);
KK_ = rrhat.*sum(rhat.*nhat,1).*wJ.*(3./(4*pi*normr.*normr));

GN_ = GG_*BasisFn';
KN_ = KK_*BasisFn';

GN = reshape(GN_([1 2 3 2 4 5 3 5 6],:), numDofPerNode, numDofPerElem);
KN = reshape(KN_([1 2 3 2 4 5 3 5 6],:), numDofPerNode, numDofPerElem);

K__ = sum(KN_,2);
K_ = K__([1 2 3 2 4 5 3 5 6]);
K = reshape(K_,numDofPerNode,numDofPerNode);