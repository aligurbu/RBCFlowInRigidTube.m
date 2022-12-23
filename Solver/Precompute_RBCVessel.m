function [RBCVessel_LHS, RBCVessel_RHS] = ...
            Precompute_RBCVessel(RBCxi, coord, connect, ...
                                 inletelem, elemDofNum, ...
                                 NeumannDofs, DirichletElem, NeumannElem, ...
                                 FieldPts, NormalV, Weights, BasisFn, ...
                                 Telem, ...
                                 grx, grw, gtx, gtw, mu, numGaussPoints, ...
                                 numNodes, numDofPerNode, numDofPerElem)
%% Similar to post-processing function for the computation of velocity 
%% field inside the vessel domain for given evaluation point, but here 
%% the evaluation point is the coordinates of sample points on RBC.
%% The nearly singular integrals on the surface of vessel are considered 
%% if the sample point on RBC is in the near vicinity. 
%%
numNodesRBC = size(RBCxi,2); % Total number of nodes on RBC surface
numInletElem = length(inletelem);
numNeumannElem = length(NeumannElem);
numDirichletElem = length(DirichletElem);
%% Build matrices
bcell = cell(numNodesRBC,1);
Acell = cell(numNodesRBC,1);
%%
parfor n = 1:numNodesRBC
    chi = RBCxi(:,n);
    tempb = zeros(numDofPerNode,1);
    tempA = zeros(numDofPerNode, numNodes*numDofPerNode);
    %%
    for mm = 1:numInletElem
        m = inletelem(mm);
        ind =  numGaussPoints^2*(m-1) + (1:numGaussPoints^2);
        xi = FieldPts(:,ind);
        wJ = Weights(ind);
        
        %%
        xi_e = coord(:,connect(:,m));
        l1 = norm(xi_e(:,6) - xi_e(:,9)) + norm(xi_e(:,9) - xi_e(:,8));
        l2 = norm(xi_e(:,7) - xi_e(:,9)) + norm(xi_e(:,9) - xi_e(:,5));
        LengE = max(l1,l2);

        [dmin, dind] = min(sqrt(sum((xi_e - chi).^2,1)));
        if dmin/LengE < 1
            zetap = closest_zeta_9nodequad(chi, xi_e, dmin, dind);
            GN = NearlySingular_ElementIntegrals_GxpN(xi_e, chi, mu, ...
                                                      grx, grw, gtx, gtw, ...
                                                      zetap(1), zetap(2));
        else
%             GN = RegularIntegrals_GN_M(chi,xi,wJ,BasisFn,mu,...
%                                        numDofPerNode,numDofPerElem);
            GN = RegularIntegrals_GN(chi,xi,wJ,BasisFn,mu);
        end
        tempb = tempb + GN * Telem(:,m);
    end
    %%
    for mm = 1:numNeumannElem
        m = NeumannElem(mm);
        ind =  numGaussPoints^2*(m-1) + (1:numGaussPoints^2);
        xi = FieldPts(:,ind);
        nhat = NormalV(:,ind);
        wJ = Weights(ind);
        
        %%
        xi_e = coord(:,connect(:,m));
        l1 = norm(xi_e(:,6) - xi_e(:,9)) + norm(xi_e(:,9) - xi_e(:,8));
        l2 = norm(xi_e(:,7) - xi_e(:,9)) + norm(xi_e(:,9) - xi_e(:,5));
        LengE = max(l1,l2);
    
        [dmin, dind] = min(sqrt(sum((xi_e - chi).^2,1)));
        if dmin/LengE < 1
            zetap = closest_zeta_9nodequad(chi, xi_e, dmin, dind);
            KNi = NearlySingular_ElementIntegrals_KxpN_M(xi_e, chi, ...
                                                         grx, grw, gtx, gtw, ...
                                                         zetap(1), zetap(2));
            % Can we improve accuracy of K when it is nearly singular.
            [~,K] = RegularIntegrals_KN_K(chi,xi,nhat,wJ,BasisFn);
            Nzetap = interpolate_9nodequad(zetap(1), zetap(2))';            
            KN = reshape(K(:)*Nzetap, numDofPerNode, numDofPerElem) + KNi;
        else
%             [KN, ~] = RegularIntegrals_KN_K_M(chi,xi,nhat,wJ, ...
%                                               BasisFn,numDofPerNode, ...
%                                               numDofPerElem);
            [KN, ~] = RegularIntegrals_KN_K(chi,xi,nhat,wJ,BasisFn);
        end
        elemDofind = ismember(elemDofNum(:,m),NeumannDofs);
        elemDofNumind = elemDofNum(elemDofind,m);
        tempA(:,elemDofNumind) = tempA(:,elemDofNumind) + KN(:,elemDofind);
    end
    %%
    for mm = 1:numDirichletElem
        m = DirichletElem(mm);
        ind =  numGaussPoints^2*(m-1) + (1:numGaussPoints^2);
        xi = FieldPts(:,ind);
        wJ = Weights(ind);
        
        %%
        xi_e = coord(:,connect(:,m));
        l1 = norm(xi_e(:,6) - xi_e(:,9)) + norm(xi_e(:,9) - xi_e(:,8));
        l2 = norm(xi_e(:,7) - xi_e(:,9)) + norm(xi_e(:,9) - xi_e(:,5));
        LengE = max(l1,l2);
        
        [dmin, dind] = min(sqrt(sum((xi_e - chi).^2,1)));
        if dmin/LengE < 1
            zetap = closest_zeta_9nodequad(chi, xi_e, dmin, dind);
            GN = NearlySingular_ElementIntegrals_GxpN(xi_e, chi, mu, ...
                                                      grx, grw, gtx, gtw, ...
                                                      zetap(1), zetap(2));
        else
%             GN = RegularIntegrals_GN_M(chi,xi,wJ,BasisFn,mu,...
%                                        numDofPerNode,numDofPerElem);
            GN = RegularIntegrals_GN(chi,xi,wJ,BasisFn,mu);
        end
        tempA(:,elemDofNum(:,m)) = tempA(:,elemDofNum(:,m)) - GN;
    end
    %%
    bcell{n} = tempb;
    Acell{n} = tempA;
end
%%
RBCVessel_RHS = cell2mat(bcell);
RBCVessel_LHS = cell2mat(Acell);