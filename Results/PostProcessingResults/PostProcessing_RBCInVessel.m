%% Vessel set-up
ParametersForTheAnalysis_Vessel

%% Field points (Gauss quadrature nodes)
[FieldPts, BasisFn, NormalV, Weights] = ...
                   FieldProperties(coord, connect, numElem, ...
                                   numDofPerNode, numNodesPerElem, gx, gw);

%% RBC set-up
ParametersForTheAnalysis_RBC
%% Note that eta is NOT needed for the nearly-singular integration!

%% Masks to go between Spherepack and vector representations of SH coeff
mask_a = repmat(triu(true(N+1),0),1,1,3);
mask_b = mask_a;
mask_b(1,:,:) = false;

%% Patch faces (connectivity matrix)
faces = PatchFaces(nlat, nlon);

%%
Numframe = 0;
for nstep = 1:NSTEPS
    if (nstep~=1 && nstep~=NSTEPS && mod(nstep,timeStepIncrement)~=0)
        continue
    end
    Numframe = Numframe + 1;
end

%%
xi_GaussX = zeros(nlat,nlon,Numframe);
xi_GaussY = zeros(nlat,nlon,Numframe);
xi_GaussZ = zeros(nlat,nlon,Numframe);
CellX = zeros(size(faces,1),size(faces,2),Numframe);
CellY = zeros(size(faces,1),size(faces,2),Numframe);
CellZ = zeros(size(faces,1),size(faces,2),Numframe);

u_GaussX = zeros(nlat,nlon,Numframe);
u_GaussY = zeros(nlat,nlon,Numframe);
u_GaussZ = zeros(nlat,nlon,Numframe);

%%
VelEvalPtsX = zeros(Numframe,size(EvaluationPts,2));
VelEvalPtsY = zeros(Numframe,size(EvaluationPts,2));
VelEvalPtsZ = zeros(Numframe,size(EvaluationPts,2));

%%
nframe = 0;
T_step = zeros(Numframe,1);
for nstep = 1:NSTEPS
    %% Read from file
    cxi = fread(fidCoord,3*(N+1)^2,'double');
    axi = zeros(size(mask_a));  bxi = zeros(size(mask_b));
    axi(mask_a) = cxi(1:3*(N+1)*(N+2)/2);
    bxi(mask_b) = cxi(3*(N+1)*(N+2)/2+1:3*(N+1)^2);

    f = reshape(fread(fidMemFor,nlat*nlon*3,'double'),nlat,nlon,3);

    Tubex = fread(fidSol,numNodes*numDofPerNode,'double');
    aU = fread(fidSol,3*(N+1)*(N+2)/2,'double');
    bU = fread(fidSol,3*N*(N+1)/2,'double');

    if (nstep~=1 && nstep~=NSTEPS && mod(nstep,timeStepIncrement)~=0)
        continue
    end
    nframe = nframe + 1;
    T_step(nframe) = Time(nstep);

    %% Collect velocities
    unodal(NeumannDofs) = Tubex(NeumannDofs);

    %% Collect tractions
    Telem(:,DirichletElem) = Tubex(elemDofNum(:,DirichletElem));

    %% Determine whether the target point, chi, is inside or outside
    % All the points outside for now.
    EvaluationPtsLocation = false(size(EvaluationPts,2),1);
    xi_Gauss = shsgcm(axi,bxi);
    xiGauss = reshape(xi_Gauss,nlat*nlon,3)';
    [gxi_thet,gxi_phi] = gradgcm(axi,bxi);
    Normal(:,:,1) = gxi_thet(:,:,2).*gxi_phi(:,:,3) - ...
                    gxi_thet(:,:,3).*gxi_phi(:,:,2);
    Normal(:,:,2) = gxi_thet(:,:,3).*gxi_phi(:,:,1) - ...
                    gxi_thet(:,:,1).*gxi_phi(:,:,3);
    Normal(:,:,3) = gxi_thet(:,:,1).*gxi_phi(:,:,2) - ...
                    gxi_thet(:,:,2).*gxi_phi(:,:,1);
    NormalGauss = reshape(Normal,nlat*nlon,3)';
    for n = 1 : size(EvaluationPts,2)
        chi = EvaluationPts(:,n);
        [~,IndGaussmin] = min(sqrt(sum((xiGauss-chi).^2,1)));
        XiGaussPts = xiGauss(:,IndGaussmin);

        NormalGaussPts = NormalGauss(:,IndGaussmin);
        if (chi - XiGaussPts)'*NormalGaussPts < 0
            % The target point inside
            EvaluationPtsLocation(n) = true;
        end
    end

    %% Compute influence of vessel to the velocity field inside the domain
    VelocityVessel = PostProcessing(EvaluationPts, coord, connect, ...
                                    inletelem, elemDofNum, ...
                                    DirichletElem, NeumannElem, ...
                                    FieldPts, NormalV, Weights, BasisFn, ...
                                    unodal, Telem, ...
                                    grx, grw, gtx, gtw, mu, numGaussPoints, ...
                                    numDofPerNode, numDofPerElem);

    VelocityVessel(:,EvaluationPtsLocation) = ...
                               VelocityVessel(:,EvaluationPtsLocation)/lam;

    %% Finer equally-spaced grid for the closest point computation
    axiEqSpaced = UpSampling(axi,N,N_EqSpaced);
    bxiEqSpaced = UpSampling(bxi,N,N_EqSpaced);
    xiEqSpaced_ = shsecm(axiEqSpaced,bxiEqSpaced);
    xiEqSpaced = reshape(xiEqSpaced_,nlatEqSpaced*nlonEqSpaced,3)';

    %% Compute the nearly-singular integral with Gf
    [af,bf] = shagcm(f);
    cf = [af(mask_a); bf(mask_b)];
    cxi = [axi(mask_a); bxi(mask_b)];
    Gf = NearlySingular_ComputeGf(cf, cxi, xiEqSpaced, ...
                                  EvaluationPts, EvaluationPtsLocation, ...
                                  RotationMatrix_EqSpaced, ...
                                  wg, mask_a, mask_b, mu, lam, ...
                                  N, NGSphere, NGphi);

	%% Compute the integrand for the nearly-singular integral with K kernel
    Kernel = NearlySingular_ComputeK(cxi, xiEqSpaced, ...
                                     EvaluationPts, ...
                                     RotationMatrix_EqSpaced, ...
                                     wg, mask_a, mask_b, ...
                                     N, NGSphere, NGthet, NGphi);

	%% Compute the nearly-singular integral with K(u-u_\perp)
    au = zeros(size(axi)); bu = zeros(size(bxi));
    au(mask_a) = aU; bu(mask_b) = bU;
    cu = [aU;bU];

    auEqSpaced = UpSampling(au,N,N_EqSpaced);
    buEqSpaced = UpSampling(bu,N,N_EqSpaced);
    uEqSpaced_ = shsecm(auEqSpaced,buEqSpaced);
    uEqSpaced = reshape(uEqSpaced_,nlatEqSpaced*nlonEqSpaced,3)';

    Ku = NearlySingular_ComputeKu(Kernel, cu, xiEqSpaced, uEqSpaced, ...
                                  EvaluationPts, EvaluationPtsLocation, ...
                                  RotationMatrix_EqSpaced, ...
                                  mask_a, mask_b, lam, ...
                                  N, NGSphere);

	%% Velocity at the evaluation points
    VelocityEvaluationPts = VelocityVessel(:) - Ku - Gf;

    VelEvalPts = reshape(VelocityEvaluationPts, size(EvaluationPts,1), ...
                                                size(EvaluationPts,2));

    %%
    VelEvalPtsX(nframe,:) = VelEvalPts(1,:);
    VelEvalPtsY(nframe,:) = VelEvalPts(2,:);
    VelEvalPtsZ(nframe,:) = VelEvalPts(3,:);

    %% Coordinates of vertices
    xi_equal = shsecm(axi,bxi);
    xi_Gauss = shsgcm(axi,bxi);
    xi_GaussX(:,:,nframe) = xi_Gauss(:,:,1);
    xi_GaussY(:,:,nframe) = xi_Gauss(:,:,2);
    xi_GaussZ(:,:,nframe) = xi_Gauss(:,:,3);
    xiplot = [xi_equal(1,:,:); xi_Gauss; xi_equal(end,:,:)];
    Vert = reshape(xiplot,(nlat+2)*nlon,3);
    CellX(:,:,nframe) = reshape(Vert(faces(:),1),size(faces));
    CellY(:,:,nframe) = reshape(Vert(faces(:),2),size(faces));
    CellZ(:,:,nframe) = reshape(Vert(faces(:),3),size(faces));

    %% Membrane velocity on patches vertices
    u_Gauss = shsgcm(au, bu);
    u_GaussX(:,:,nframe) = u_Gauss(:,:,1);
    u_GaussY(:,:,nframe) = u_Gauss(:,:,2);
    u_GaussZ(:,:,nframe) = u_Gauss(:,:,3);
end
save(['PostProc_',name,'.mat'], ...
      'T_step', ...
      'VelEvalPtsX', 'VelEvalPtsY', 'VelEvalPtsZ',...
      'xi_GaussX', 'xi_GaussY', 'xi_GaussZ', ...
      'u_GaussX', 'u_GaussY', 'u_GaussZ', ...
      'CellX', 'CellY', 'CellZ')