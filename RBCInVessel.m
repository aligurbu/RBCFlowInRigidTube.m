%% Red blood cell flow inside a vessel
% Spectral Galerkin BIE for an RBC motion and deformation and
% direct BEM for the rigid vessel implementations are put together
% to analyze the stresses develop on the cell membrane while flowing
% in a vessel with a size comparable to red cell or while squeezing
% through a constriction.
%%
clear all; close all; clc;
addpath(genpath('../RBCFlowInRigidTube.m'))
verbose_Plot = false;

%% Input the model and parameters for the analysis from Models folder
LoadElasRBC_Short_Pr4_2_Time0_75s
% LoadElasRBC_RefCons_6mic_Pr8
% LoadElasRBC_LongConVes_Pr8

% LoadMemVisRBC_Short_muMem10_Pr4_2_Time0_75
% LoadMemVisRBC_RefCons_6mic_muMem_3_18_Pr8
% LoadMemVisRBC_LongConVes_muMem_3_18_P40

%% Vessel set-up
ParametersForTheAnalysis_Vessel

%% RBC set-up
ParametersForTheAnalysis_RBC

%% Check the geometry and position of RBC
if verbose_Plot
    figure('Color','white')
    hold on
    VisualizeGeometry(nlat, nlon, aXi, bXi, 'r', true)
    Patch_Mesh(coord, connect, 0.1)
    axis on
    xlabel('X'); ylabel('Y'); zlabel('Z');
    view([0 90])
end

%% Field points (Gauss quadrature nodes)
[FieldPts, BasisFn, NormalV, Weights] = ...
                   FieldProperties(coord, connect, numElem, ...
                                   numDofPerNode, numNodesPerElem, gx, gw);

%% Precomputations Vessel
[VesselVessel_LHS, VesselVessel_RHS, VesselVessel_Preconditioner] = ...
             PrecomputeBEM_Vessel(coord, connect, inletelem, outletelem, ...
                                  elemDofNum, ...
                                  NeumannDofs, NeumannNode, DirichletElem, ...
                                  FieldPts, NormalV, Weights, BasisFn, ...
                                  Telem, ...
                                  grx, grw, gtx, gtw, mu, numGaussPoints, ...
                                  numNodes, numDofPerNode);

%% Masks to go between Spherepack and vector representations of SH coeff
mask_a = repmat(triu(true(N+1),0),1,1,3);
mask_b = mask_a;
mask_b(1,:,:) = false;

%% Initialization
xi = Xi; axi = aXi; bxi = bXi;

viscousStress_prev = zeros(UpSampleFactor*N+1, 2*UpSampleFactor*N+1, 4);
epsilbrev_prev = zeros(UpSampleFactor*N+1, 2*UpSampleFactor*N+1, 4);

%% Set up time integration
Solution = zeros(ModelSize,1);

%% Set up output file
fidTime = fopen(['Time_',name,'.dat'],'w');
fidCoord = fopen(['Coord_',name,'.dat'],'w');
fidMemFor = fopen(['MemFor_',name,'.dat'],'w');
fidSol = fopen(['Sol_',name,'.dat'],'w');

%% Time-stepping
for nstep = 0:NSTEPS
    if (nstep==0 || mod(nstep,SaveAtIncrementalSteps)==0)
        %% Write the time to file
        fwrite(fidTime, Time, 'double');
        %% Write the position of cell to file
        cxi = [axi(mask_a); bxi(mask_b)];
        fwrite(fidCoord, cxi, 'double');
    end

    [u, au, bu, Solution, f, viscousStress_prev, epsilbrev_prev, ITER] = ...
         SolveSpectralBIE_BEM(VesselVessel_LHS, VesselVessel_RHS, ...
                              VesselVessel_Preconditioner, ...
                              coord, connect, coordLocation, ...
                              numDofPerElem, elemDofNum, NeumannDofs, ...
                              inletelem, DirichletElem, NeumannElem, ...
                              FieldPts, NormalV, Weights, BasisFn, ...
                              Telem, ...
                              EE, FF, GG, WW, JXibrev, ...
                              LL, MM, NN, ...
                              UpSampleFactor, ...
                              ES, ED, EB, ...
                              StVenantKirchhoff, ...
                              NeoHookean, Skalak, ...
                              MembraneViscoelasticity, ...
                              mu_Mem, Tau, DT, ...
                              viscousStress_prev, ...
                              epsilbrev_prev, ...
                              axi, bxi, xi, ...
                              RotationMatrix, ...
                              nlatEqSpaced, nlonEqSpaced, N_EqSpaced, ...
                              RotationMatrix_EqSpaced, ...
                              NGSphere, NGthet, NGphi, ...
                              eta, wg, ...
                              nlat, nlon, mask_a, mask_b, ...
                              mu, lam, N, ...
                              grx, grw, gtx, gtw, numGaussPoints, ...
                              numNodes, numDofPerNode, ModelSize, ...
                              ToleranceGMRES, Solution);

    if (nstep==0 || mod(nstep,SaveAtIncrementalSteps)==0)
        %% Write output to file
        writef = f(:)';
        fwrite(fidMemFor,writef,'double');
        writeSolution = Solution';
        fwrite(fidSol,writeSolution,'double');
    end

    %% Display
    f_norm = sqrt(reshape(f(:,:,1).^2+f(:,:,2).^2+f(:,:,3).^2,nlat,nlon));
    u_norm = sqrt(reshape(u(:,:,1).^2+u(:,:,2).^2+u(:,:,3).^2,nlat,nlon));

    fprintf('# of step %d;\t remaining # of steps %d; \n', ...
             nstep, NSTEPS-nstep)
    fprintf('DT: %g;\t Time %d \n', DT, Time);
    fprintf('normf = %g;\t normu = %g\n',max(max(f_norm)),max(max(u_norm)))
    fprintf('Position of RBC min: %g max: %g \n', ...
             min(min(xi(:,:,1))), max(max(xi(:,:,1))))
    fprintf('Number of GMRES iteration: %d\n', ITER(2))
    fprintf('# of stableCounter %d\n', stableCounter)
    fprintf('# of counterTime %d\n', counterTime)
    fprintf('\n')

    %% Update state (Forward Euler time scheme)
    axi = axi + DT*au;
    bxi = bxi + DT*bu;
    xi = shsgcm(axi,bxi);

    if constrainCenterOfMassToInitialPosition
        %% Set the center of the mass of RBC to initial position, InitXi
        [gxi_thet, gxi_phi] = gradgcm(axi, bxi);
        normalVector = cross(gxi_thet, gxi_phi, 3);
        xi_dot_normalVector = sum(xi.*normalVector,3);
        xi_xi_dot_normalVector = xi.*xi_dot_normalVector;
        V = sum(sum(xi_dot_normalVector.*weight))*(2*pi/nlon)/3;
        Int_xi_xi_dot_normalVector_dA = ...
                    sum(sum(xi_xi_dot_normalVector.*weight))*(2*pi/nlon)/4;
        CenterMass = reshape(Int_xi_xi_dot_normalVector_dA./V,3,1);
        xi = xi - reshape((CenterMass - InitXi),1,1,3);
    end
    [axi,bxi] = shagcm(xi);

end