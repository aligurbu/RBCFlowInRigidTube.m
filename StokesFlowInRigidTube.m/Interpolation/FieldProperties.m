function [FieldPts, BasisFn, NormalV, Weights] = ...
                    FieldProperties(coord, connect, numElem, ...
                                    numDofPerNode, numNodesPerElem, gx, gw)

NG = length(gx);

%% Field points (Gauss Quadrature nodes)
BasisFn = zeros(numNodesPerElem,NG^2);
FieldPts = zeros(numDofPerNode,numElem*NG^2);
NormalV = zeros(numDofPerNode,numElem*NG^2);
Weights = zeros(1, numElem*NG^2);
for s1 = 1:NG
    for s2 = 1:NG
        Nind = NG*(s1-1) + s2;
        BasisFn(:,Nind) = interpolate_9nodequad(gx(s1),gx(s2));
    end
end
for m = 1:numElem
    xi_e = coord(:,connect(:,m));
    for s1 = 1:NG
        for s2 = 1:NG
            ind = NG^2*(m-1) + NG*(s1-1) + s2;
            Nind = NG*(s1-1) + s2;
            FieldPts(:,ind) = xi_e * BasisFn(:,Nind);
            [nhat, JJ] = jacobian_9nodequad(gx(s1), gx(s2), xi_e);
            NormalV(:,ind) = nhat;
            Weights(ind) = JJ*gw(s1)*gw(s2);
        end
    end
end