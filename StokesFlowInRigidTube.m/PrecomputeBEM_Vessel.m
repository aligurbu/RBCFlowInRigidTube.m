function [A, b, Minv] = ...
             PrecomputeBEM_Vessel(coord, connect, inletelem, outletelem, ...
                                  elemDofNum, ...
                                  NeumannDofs, NeumannNode, DirichletElem, ...
                                  FieldPts, NormalV, Weights, BasisFn, ...
                                  Telem, ...
                                  grx, grw, gtx, gtw, mu, numGaussPoints, ...
                                  numNodes, numDofPerNode, numDofPerElem)
%% Build matrices
bcell = cell(numNodes,1);
Acell = cell(numNodes,1);
Mcell = cell(numNodes,1);

numDirichletElem = length(DirichletElem);
numInletElem = length(inletelem);
numOutletElem = length(outletelem);
%%
found_mex_RegularIntegrals_GN_K = ...
                      ~isempty(which('Integration/RegularIntegrals_GN_K'));
found_mex_RegularIntegrals_GN_KN_K = ...
                   ~isempty(which('Integration/RegularIntegrals_GN_KN_K'));
found_mex_RegularIntegrals_KN_K = ...
                      ~isempty(which('Integration/RegularIntegrals_KN_K'));
found_mex_WeaklySingularIntegrals_GxN_KxN = ...
            ~isempty(which('Integration/WeaklySingularIntegrals_GxN_KxN'));
found_mex_WeaklySingularIntegrals_KxN = ...
                ~isempty(which('Integration/WeaklySingularIntegrals_KxN'));
found_mex_WeaklySingularIntegrals_GxN = ...
                ~isempty(which('Integration/WeaklySingularIntegrals_GxN'));
parfor n = 1:numNodes
    nodeDofNum = (n-1)*numDofPerNode + (1:numDofPerNode);
    NeumannNodeind = ismember(n,NeumannNode);
    chi = coord(:,n);
    tempb = zeros(numDofPerNode,1);
    tempA = zeros(numDofPerNode, numNodes*numDofPerNode);
    tempM = zeros(numDofPerNode, numNodes*numDofPerNode);
    %%
    for mm = 1:numInletElem
        m = inletelem(mm);
        ind =  numGaussPoints^2*(m-1) + (1:numGaussPoints^2);
        xi = FieldPts(:,ind);
        nhat = NormalV(:,ind);
        wJ = Weights(ind);

        %% Check if the node (n) belongs to the element (m)
        %% find the position of the node in the element: xnodenum
        xnodenum = find(connect(:,m)==n);
        if (isempty(xnodenum))
            % if the node does not belong to the element (m)
            xnodenum = 0;
        end
        %%
        if xnodenum == 0
            if found_mex_RegularIntegrals_GN_KN_K
                [GN, KN, K] = RegularIntegrals_GN_KN_K(chi,xi,nhat,wJ, ...
                                                       BasisFn, mu);
            else
                [GN, KN, K] = RegularIntegrals_GN_KN_K_M(chi,xi,nhat,wJ, ...
                                                         BasisFn, mu,...
                                                         numDofPerNode, ...
                                                         numDofPerElem);
            end
            tempb = tempb + GN * Telem(:,m);
            elemDofind = ismember(elemDofNum(:,m),NeumannDofs);
            elemDofNumind = elemDofNum(elemDofind,m);
            tempA(:,elemDofNumind) = ...
                                tempA(:,elemDofNumind) + KN(:,elemDofind);
            if NeumannNodeind
                tempA(:,nodeDofNum) = tempA(:,nodeDofNum) - K;
                tempM(:,nodeDofNum) = tempM(:,nodeDofNum) - K;
            end
        else
            xi_e = coord(:,connect(:,m));
            if found_mex_WeaklySingularIntegrals_GxN_KxN
                [GxN, KxN] = WeaklySingular_ElementIntegrals_GxN_KxN ...
                            (chi, xi_e, mu, xnodenum, grx, grw, gtx, gtw);
            else
                [GxN, KxN] = WeaklySingular_ElementIntegrals_GxN_KxN_M ...
                            (chi, xi_e, mu, xnodenum, grx, grw, gtx, gtw);
            end
            tempb = tempb + GxN * Telem(:,m);
            elemDofind = ismember(elemDofNum(:,m),NeumannDofs);
            elemDofNumind = elemDofNum(elemDofind,m);
            tempA(:,elemDofNumind) = ...
                                tempA(:,elemDofNumind) + KxN(:,elemDofind);
            tempM(:,elemDofNumind) = ...
                                tempM(:,elemDofNumind) + KxN(:,elemDofind);
        end
    end
    %%
    for mm = 1:numOutletElem
        m = outletelem(mm);
        ind =  numGaussPoints^2*(m-1) + (1:numGaussPoints^2);
        xi = FieldPts(:,ind);
        nhat = NormalV(:,ind);
        wJ = Weights(ind);

        %% Check if the node (n) belongs to the element (m)
        %% find the position of the node in the element: xnodenum
        xnodenum = find(connect(:,m)==n);
        if (isempty(xnodenum))
            % if the node does not belong to the element (m)
            xnodenum = 0;
        end
        %%
        if xnodenum == 0
            if found_mex_RegularIntegrals_KN_K
                [KN, K] = RegularIntegrals_KN_K(chi,xi,nhat,wJ,BasisFn);
            else
                [KN, K] = RegularIntegrals_KN_K_M(chi,xi,nhat,wJ, ...
                                                  BasisFn,numDofPerNode, ...
                                                  numDofPerElem);
            end
            elemDofind = ismember(elemDofNum(:,m),NeumannDofs);
            elemDofNumind = elemDofNum(elemDofind,m);
            tempA(:,elemDofNumind) = ...
                                tempA(:,elemDofNumind) + KN(:,elemDofind);
            if NeumannNodeind
                tempA(:,nodeDofNum) = tempA(:,nodeDofNum) - K;
                tempM(:,nodeDofNum) = tempM(:,nodeDofNum) - K;
            end
        else
            xi_e = coord(:,connect(:,m));
            if found_mex_WeaklySingularIntegrals_KxN
                KxN = WeaklySingular_ElementIntegrals_KxN...
                                                      (chi, xi_e, ...
                                                       xnodenum, ...
                                                       grx, grw, gtx, gtw);
            else
                KxN = WeaklySingular_ElementIntegrals_KxN_M...
                                                      (chi, xi_e, ...
                                                       xnodenum, ...
                                                       grx, grw, gtx, gtw);
            end
            elemDofind = ismember(elemDofNum(:,m),NeumannDofs);
            elemDofNumind = elemDofNum(elemDofind,m);
            tempA(:,elemDofNumind) = ...
                                tempA(:,elemDofNumind) + KxN(:,elemDofind);
            tempM(:,elemDofNumind) = ...
                                tempM(:,elemDofNumind) + KxN(:,elemDofind);
        end
    end
    %%
    for mm = 1:numDirichletElem
        m = DirichletElem(mm);
        ind =  numGaussPoints^2*(m-1) + (1:numGaussPoints^2);
        xi = FieldPts(:,ind);
        nhat = NormalV(:,ind);
        wJ = Weights(ind);

        %% Check if the node (n) belongs to the element (m)
        %% find the position of the node in the element: xnodenum
        xnodenum = find(connect(:,m)==n);
        if (isempty(xnodenum))
            % if the node does not belong to the element (m)
            xnodenum = 0;
        end
        %%
        if xnodenum == 0
            if found_mex_RegularIntegrals_GN_K
                [GN, K] = RegularIntegrals_GN_K(chi,xi,nhat,wJ,BasisFn,mu);
            else
                [GN, K] = RegularIntegrals_GN_K_M(chi,xi,nhat,wJ,BasisFn,mu,...
                                                  numDofPerNode,numDofPerElem);
            end
            tempA(:,elemDofNum(:,m)) = tempA(:,elemDofNum(:,m)) - GN;
            if NeumannNodeind
                tempA(:,nodeDofNum) = tempA(:,nodeDofNum) - K;
                tempM(:,nodeDofNum) = tempM(:,nodeDofNum) - K;
            end
        else
            xi_e = coord(:,connect(:,m));
            if found_mex_WeaklySingularIntegrals_GxN
                GxN = WeaklySingular_ElementIntegrals_GxN...
                                                      (chi, xi_e, mu, ...
                                                       xnodenum, ...
                                                       grx, grw, gtx, gtw);
            else
                GxN = WeaklySingular_ElementIntegrals_GxN_M...
                                                      (chi, xi_e, mu, ...
                                                       xnodenum, ...
                                                       grx, grw, gtx, gtw);
            end
            tempA(:,elemDofNum(:,m)) = tempA(:,elemDofNum(:,m)) - GxN;
            tempM(:,elemDofNum(:,m)) = tempM(:,elemDofNum(:,m)) - GxN;
        end
    end
    %%
    bcell{n} = tempb;
    Acell{n} = tempA;
    Mcell{n} = tempM;
end
%%
b = cell2mat(bcell);
A = cell2mat(Acell);
MTube = cell2mat(Mcell);

Minv = inv(MTube);