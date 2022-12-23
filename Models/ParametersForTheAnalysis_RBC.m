%% RBC set-up
%% Create grid for representation: [N+1 x 2*N+1]
[nlat, nlon, thet, phi, weight] = GridOnSphere(N);

%% Create grid for integration points on unit sphere:
[NGthet, NGphi, eta, wg] = setup_integration_grid(N, NGSphere);

%% Create geometry of RBC: undeformed geometry of an RBC
Xi = getRBCInitialGeometry(thet,phi,InitXi,InitOrient);
% Xi = getSphereGeometry(thet, phi, Radius, InitXi, InitOrient);

%% Initial conditions:
%% spherical harmonics coefficients of the undeformed RBC coordinates
[aXi,bXi] = shagcm(Xi);

%% Compute the first and second fundamental form coefficients
[~, ~, EE, FF, GG, WW, JXibrev, ~, ~, ~, LL, MM, NN] = ...
                   coefficientsOfFundamentalForm(aXi, bXi, UpSampleFactor);

%% Create rotation matrices
RotationMatrix = RotationMatrixForSHCoefficients(N, nlat, nlon, thet, phi);

%% Create fine equally-spaced grid for nearly-singular integrals
[nlatEqSpaced, nlonEqSpaced, phiEqSpaced, thetEqSpaced] = ...
                                     EquallySpacedGridOnSphere(N_EqSpaced);

%% Create rotation matrices for equally-spaced grid
RotationMatrix_EqSpaced = ...
         RotationMatrixForSHCoefficients(N, nlatEqSpaced, nlonEqSpaced, ...
                                         phiEqSpaced, thetEqSpaced);