function UU = PostProcessing(EvaluationPts, coord, connect, ...
                             inletelem, elemDofNum, ...
                             DirichletElem, NeumannElem, ...
                             FieldPts, NormalV, Weights, BasisFn, ...
                             unodal, Telem, ...
                             grx, grw, gtx, gtw, mu, numGaussPoints, ...
                             numDofPerNode, numDofPerElem)
%% Post-Processing
%%
numPts = size(EvaluationPts,2); % Total number of points
numDirichletElem = length(DirichletElem);
numNeumannElem = length(NeumannElem);
%%
UU = zeros(numPts*numDofPerNode,1);

%%
found_mex_RegularIntegrals_GN = ...
                        ~isempty(which('Integration/RegularIntegrals_GN'));
found_mex_RegularIntegrals_KN_K = ...
                      ~isempty(which('Integration/RegularIntegrals_KN_K'));
found_mex_NearlySingularIntegrals_GxpN = ...
               ~isempty(which('Integration/NearlySingularIntegrals_GxpN'));
found_mex_NearlySingularIntegrals_KxpN = ...
               ~isempty(which('Integration/NearlySingularIntegrals_KxpN'));
parfor mm = 1:length(inletelem)
    m = inletelem(mm);
    ind =  numGaussPoints^2*(m-1) + (1:numGaussPoints^2);
    xi = FieldPts(:,ind);
    wJ = Weights(ind);

    GNInlet = zeros(numPts*numDofPerNode,numDofPerElem);

    xi_e = coord(:,connect(:,m));
    l1 = norm(xi_e(:,6) - xi_e(:,9)) + norm(xi_e(:,9) - xi_e(:,8));
    l2 = norm(xi_e(:,7) - xi_e(:,9)) + norm(xi_e(:,9) - xi_e(:,5));
    LengE = max(l1,l2);

    for c = 1:numPts
        nodeDofNum = (c-1)*numDofPerNode + (1:numDofPerNode);
        chi = EvaluationPts(:,c);

        [dmin, dind] = min(sqrt(sum((xi_e - chi).^2,1)));
        if dmin/LengE < 1
            zetap = closest_zeta_9nodequad(chi, xi_e, dmin, dind);
            if found_mex_NearlySingularIntegrals_GxpN
                GN = NearlySingular_ElementIntegrals_GxpN...
                                                   (xi_e, chi, mu, ...
                                                    grx, grw, gtx, gtw, ...
                                                    zetap(1), zetap(2));
            else
                GN = NearlySingular_ElementIntegrals_GxpN_M...
                                                   (xi_e, chi, mu, ...
                                                    grx, grw, gtx, gtw, ...
                                                    zetap(1), zetap(2));
            end
        else
            if found_mex_RegularIntegrals_GN
                GN = RegularIntegrals_GN(chi,xi,wJ,BasisFn,mu);
            else
                GN = RegularIntegrals_GN_M(chi,xi,wJ,BasisFn,mu,...
                                           numDofPerNode,numDofPerElem);
            end
        end
        GNInlet(nodeDofNum, :) = GN;
    end
    UU = UU + GNInlet * Telem(:,m);
end
%%
parfor mm = 1:numNeumannElem
    m = NeumannElem(mm);
    ind =  numGaussPoints^2*(m-1) + (1:numGaussPoints^2);
    xi = FieldPts(:,ind);
    nhat = NormalV(:,ind);
    wJ = Weights(ind);

    KNNeuman = zeros(numPts*numDofPerNode,numDofPerElem);

    xi_e = coord(:,connect(:,m));
    l1 = norm(xi_e(:,6) - xi_e(:,9)) + norm(xi_e(:,9) - xi_e(:,8));
    l2 = norm(xi_e(:,7) - xi_e(:,9)) + norm(xi_e(:,9) - xi_e(:,5));
    LengE = max(l1,l2);

    for c = 1:numPts
        nodeDofNum = (c-1)*numDofPerNode + (1:numDofPerNode);
        chi = EvaluationPts(:,c);

        [dmin, dind] = min(sqrt(sum((xi_e - chi).^2,1)));
        if dmin/LengE < 1
            zetap = closest_zeta_9nodequad(chi, xi_e, dmin, dind);
            if found_mex_NearlySingularIntegrals_KxpN
                KNi = NearlySingular_ElementIntegrals_KxpN...
                                                   (xi_e, chi, ...
                                                    grx, grw, gtx, gtw, ...
                                                    zetap(1), zetap(2));
            else
                KNi = NearlySingular_ElementIntegrals_KxpN_M...
                                                   (xi_e, chi, ...
                                                    grx, grw, gtx, gtw, ...
                                                    zetap(1), zetap(2));
            end
            % Can we improve accuracy of K when it is nearly singular.
            if found_mex_RegularIntegrals_KN_K
                [~,K] = RegularIntegrals_KN_K(chi,xi,nhat,wJ,BasisFn);
            else
                [~,K] = RegularIntegrals_KN_K_M(chi,xi,nhat,wJ, ...
                                                BasisFn,numDofPerNode, ...
                                                numDofPerElem)
            end
            Nzetap = interpolate_9nodequad(zetap(1), zetap(2))';
            KN = reshape(K(:)*Nzetap, numDofPerNode, numDofPerElem) + KNi;
        else
            if found_mex_RegularIntegrals_KN_K
                [KN, ~] = RegularIntegrals_KN_K(chi,xi,nhat,wJ,BasisFn);
            else
                [KN, ~] = RegularIntegrals_KN_K_M(chi,xi,nhat,wJ, ...
                                                  BasisFn,numDofPerNode, ...
                                                  numDofPerElem)    
            end
        end
        KNNeuman(nodeDofNum, :) = KN;
    end
    UU = UU - KNNeuman * unodal(elemDofNum(:,m));
end
%%
parfor mm = 1:numDirichletElem
    m = DirichletElem(mm);
    ind =  numGaussPoints^2*(m-1) + (1:numGaussPoints^2);
    xi = FieldPts(:,ind);
    wJ = Weights(ind);

    GNDirichlet = zeros(numPts*numDofPerNode,numDofPerElem);

    xi_e = coord(:,connect(:,m));
    l1 = norm(xi_e(:,6) - xi_e(:,9)) + norm(xi_e(:,9) - xi_e(:,8));
    l2 = norm(xi_e(:,7) - xi_e(:,9)) + norm(xi_e(:,9) - xi_e(:,5));
    LengE = max(l1,l2);

    for c = 1:numPts
        nodeDofNum = (c-1)*numDofPerNode + (1:numDofPerNode);
        chi = EvaluationPts(:,c);

        [dmin, dind] = min(sqrt(sum((xi_e - chi).^2,1)));
        if dmin/LengE < 1
            zetap = closest_zeta_9nodequad(chi, xi_e, dmin, dind);
            if found_mex_NearlySingularIntegrals_GxpN
                GN = NearlySingular_ElementIntegrals_GxpN...
                                                   (xi_e, chi, mu, ...
                                                    grx, grw, gtx, gtw, ...
                                                    zetap(1), zetap(2));
            else
                GN = NearlySingular_ElementIntegrals_GxpN_M...
                                                   (xi_e, chi, mu, ...
                                                    grx, grw, gtx, gtw, ...
                                                    zetap(1), zetap(2));
            end
        else
            if found_mex_RegularIntegrals_GN
                GN = RegularIntegrals_GN(chi,xi,wJ,BasisFn,mu);
            else
                GN = RegularIntegrals_GN_M(chi,xi,wJ,BasisFn,mu,...
                                           numDofPerNode,numDofPerElem);
            end
        end
        GNDirichlet(nodeDofNum, :) = GN;
    end
    UU = UU + GNDirichlet * Telem(:,m);
end
UU = reshape(UU,size(EvaluationPts,1),size(EvaluationPts,2));