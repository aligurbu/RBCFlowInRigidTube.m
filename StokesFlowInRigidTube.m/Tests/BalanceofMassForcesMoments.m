function [inflow, outflow, inletForce, SumofForces, ...
          SumofMoments, SumofNormofMoments] = ...
          BalanceofMassForcesMoments(FieldPts, BasisFn, NormalV, Weights, ...
                                     inletelem, outletelem, wallelem, ...
                                     elemDofNum, ...                            
                                     unodal, Telem, ...
                                     numGaussPoints)
%% Check conservation of mass and forces

%% Inlet elements
% Inflow 
ind =  numGaussPoints^2*(inletelem-1) + (1:numGaussPoints^2)';
inflowFieldx = (unodal(elemDofNum(1:3:27,inletelem))'*BasisFn)';
inflowFieldy = (unodal(elemDofNum(2:3:27,inletelem))'*BasisFn)';
inflowFieldz = (unodal(elemDofNum(3:3:27,inletelem))'*BasisFn)';
inflowField = [(inflowFieldx(:))'; ...
               (inflowFieldy(:))'; ...
               (inflowFieldz(:))'];
inflow = sum(sum(inflowField.*(-NormalV(:,ind(:)))).*Weights(1,ind(:)),2);

% inlet force and moment
inletForceFieldx = (Telem(1:3:27,inletelem)'*BasisFn)';
inletForceFieldy = (Telem(2:3:27,inletelem)'*BasisFn)';
inletForceFieldz = (Telem(3:3:27,inletelem)'*BasisFn)';
inletForceField = [(inletForceFieldx(:))'; ...
                   (inletForceFieldy(:))'; ...
                   (inletForceFieldz(:))'];
inletMomentField = cross(FieldPts(:,ind(:)), inletForceField,1);
normInletMomentField = sqrt(inletMomentField(1,:).^2 + ...
                            inletMomentField(2,:).^2 + ...
                            inletMomentField(3,:).^2);

inletForce = sum(inletForceField.*Weights(1,ind(:)),2);
inletMoment = sum(inletMomentField.*Weights(1,ind(:)),2);
normInletMoment = sum(normInletMomentField.*Weights(1,ind(:)),2);

%% Outlet elements
% Outflow
ind =  numGaussPoints^2*(outletelem-1) + (1:numGaussPoints^2)';
outflowFieldx = (unodal(elemDofNum(1:3:27,outletelem))'*BasisFn)';
outflowFieldy = (unodal(elemDofNum(2:3:27,outletelem))'*BasisFn)';
outflowFieldz = (unodal(elemDofNum(3:3:27,outletelem))'*BasisFn)';
outflowField = [(outflowFieldx(:))'; ...
                (outflowFieldy(:))'; ...
                (outflowFieldz(:))'];
outflow = sum(sum(outflowField.*(NormalV(:,ind(:)))).*Weights(1,ind(:)),2);

% outlet force and moment
outletForceFieldx = (Telem(1:3:27,outletelem)'*BasisFn)';
outletForceFieldy = (Telem(2:3:27,outletelem)'*BasisFn)';
outletForceFieldz = (Telem(3:3:27,outletelem)'*BasisFn)';
outletForceField = [(outletForceFieldx(:))'; ...
                    (outletForceFieldy(:))'; ...
                    (outletForceFieldz(:))'];
outletMomentField = cross(FieldPts(:,ind(:)), outletForceField,1);
normOutletMomentField = sqrt(outletMomentField(1,:).^2 + ...
                             outletMomentField(2,:).^2 + ...
                             outletMomentField(3,:).^2);

outletForce = sum(outletForceField.*Weights(1,ind(:)),2);
outletMoment = sum(outletMomentField.*Weights(1,ind(:)),2);
normOutletMoment = sum(normOutletMomentField.*Weights(1,ind(:)),2);

%% Wall elements
ind =  numGaussPoints^2*(wallelem-1) + (1:numGaussPoints^2)';
wallForceFieldx = (Telem(1:3:27,wallelem)'*BasisFn)';
wallForceFieldy = (Telem(2:3:27,wallelem)'*BasisFn)';
wallForceFieldz = (Telem(3:3:27,wallelem)'*BasisFn)';
wallForceField = [(wallForceFieldx(:))'; ...
                  (wallForceFieldy(:))'; ...
                  (wallForceFieldz(:))'];
wallMomentField = cross(FieldPts(:,ind(:)), wallForceField,1);
normWallMomentField = sqrt(wallMomentField(1,:).^2 + ...
                           wallMomentField(2,:).^2 + ...
                           wallMomentField(3,:).^2);

wallForce = sum(wallForceField.*Weights(1,ind(:)),2);
wallMoment = sum(wallMomentField.*Weights(1,ind(:)),2);
normWallMoment = sum(normWallMomentField.*Weights(1,ind(:)),2);

%% Sum of forces and moments
SumofForces = inletForce + outletForce + wallForce;
SumofMoments = inletMoment + outletMoment + wallMoment;
SumofNormofMoments = normInletMoment + normOutletMoment + normWallMoment;