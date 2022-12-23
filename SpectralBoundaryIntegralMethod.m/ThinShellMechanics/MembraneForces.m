function [p, viscousStress_prev, epsilbrev_prev] = ...
                           MembraneForces(EE, FF, GG, WW, JXibrev, ...
                                          LL, MM, NN, ...
                                          axi, bxi, ...
                                          UpSampleFactor, ...
                                          ES, ED, EB, ...
                                          StVenantKirchhoff, ...
                                          NeoHookean, Skalak, ...
                                          MembraneViscoelasticity, ...
                                          mu_Mem, Tau, DT, ...
                                          viscousStress_prev, ...
                                          epsilbrev_prev)
%% Input
% EE, FF, GG = the first fundamental form of coefficients
% WW, JXibrev = surface area and Jacobian determinant
% LL, MM, NN = the second fundamental form of coefficients
% axi,bxi = spherical harmonic transform of deformed geometry coordinates
% upSampleFact = upsampling factor for dealising
% ES, ED = inplane-shear and inplane-dilatational moduli of membrane
% EB = bending modulus of membrane
% mu_Mem = membrane viscoelasticity
% Tau = the relaxation time of the Maxwell element
% DT = time increment
% viscoStress_prev = viscoelastic membrane stresses from previous step
% epsilbrev_prev = strain from previous step
% p = viscoelastic and elastic forces per unit sphere area
%%
%% Compute the first and second fundamental form coefficients
[gxi_thet, gxi_phi, ee, ff, gg, ~, Jxibrev, ...
 unitNormalxi, gUnitNormalxi_thet, gUnitNormalxi_phi, ll, mm, nn] = ...
                   coefficientsOfFundamentalForm(axi, bxi, UpSampleFactor);

%% Strain metric
epsilbrev = strainMetric(EE, GG, FF, WW, ee, gg, ff);

%% Curvature metric
kappabrev = curvatureMetric(EE, FF, GG, LL, MM, NN, WW, ll, mm, nn);

%% Bending stresses
% mbrev = bendingStress(EB, kappabrev);

%% St. Venant Kirchhoff model
if StVenantKirchhoff
    nbrev = StVenantKirchoffConstitutiveEquation(ED, ES, epsilbrev);
end

%% Neo-Hookean model
if NeoHookean
    nbrev = neoHookeanConstitutiveEquation(ED, ES, epsilbrev);
end

%% Skalak constitutive equation
if Skalak
    nbrev = SkalakConstitutiveEquation(ED, ES, epsilbrev);
end

%% Membrane viscoelasticity using standard linear solid model
if MembraneViscoelasticity
    [nbrev, viscousStress_prev, epsilbrev_prev] = ...
                    StandardLinearSolidModel(nbrev, viscousStress_prev, ...
                                             epsilbrev, epsilbrev_prev, ...
                                             mu_Mem, Tau, DT);
end

%%
mbrev_Contravariant(:,:,1) = GG.*kappabrev(:,:,1) - FF.*kappabrev(:,:,3);   % mbrev_thet_thet_Contravariant
mbrev_Contravariant(:,:,2) = GG.*kappabrev(:,:,2) - FF.*kappabrev(:,:,4);   % mbrev_phi_thet_Contravariant
% mbrev_Contravariant = EE.*kappabrev(:,:,3) - FF.*kappabrev(:,:,1); % mbrev_thet_phi_Contravariant
% mbrev_phi_thet_Contravariant == mbrev_thet_phi_Contravariant
mbrev_Contravariant(:,:,3) = EE.*kappabrev(:,:,4) - FF.*kappabrev(:,:,2);   % mbrev_phi_phi_Contravariant
q_thet = EB*(gxi_thet.*mbrev_Contravariant(:,:,1) + ...
             gxi_phi.*mbrev_Contravariant(:,:,2))./JXibrev;
q_phi = EB*(gxi_thet.*mbrev_Contravariant(:,:,2) + ...
            gxi_phi.*mbrev_Contravariant(:,:,3))./JXibrev;

[br,bi,~,~] = vhagcm(q_thet,q_phi);
q = divgcm(br,bi);

Pq = q - sum(unitNormalxi.*q, 3).*unitNormalxi;

%%
gamma_theta(:,:,1) = (gxi_phi(:,:,2).*Pq(:,:,3) - ...
                      gxi_phi(:,:,3).*Pq(:,:,2))./Jxibrev;
gamma_theta(:,:,2) = (gxi_phi(:,:,3).*Pq(:,:,1) - ...
                      gxi_phi(:,:,1).*Pq(:,:,3))./Jxibrev;
gamma_theta(:,:,3) = (gxi_phi(:,:,1).*Pq(:,:,2) - ...
                      gxi_phi(:,:,2).*Pq(:,:,1))./Jxibrev;
gamma_phi(:,:,1) = (Pq(:,:,2).*gxi_thet(:,:,3) - ...
                    Pq(:,:,3).*gxi_thet(:,:,2))./Jxibrev;
gamma_phi(:,:,2) = (Pq(:,:,3).*gxi_thet(:,:,1) - ...
                    Pq(:,:,1).*gxi_thet(:,:,3))./Jxibrev;
gamma_phi(:,:,3) = (Pq(:,:,1).*gxi_thet(:,:,2) - ...
                    Pq(:,:,2).*gxi_thet(:,:,1))./Jxibrev;

%%
mu_theta = EB*(gUnitNormalxi_thet.*mbrev_Contravariant(:,:,1) + ...
               gUnitNormalxi_phi.*mbrev_Contravariant(:,:,2))./JXibrev;
mu_phi = EB*(gUnitNormalxi_thet.*mbrev_Contravariant(:,:,2) + ...
             gUnitNormalxi_phi.*mbrev_Contravariant(:,:,3))./JXibrev;

%%
nbrev_Contravariant(:,:,1) = GG.*nbrev(:,:,1) - FF.*nbrev(:,:,2);   % nbrev_thet_thet_Contravariant
nbrev_Contravariant(:,:,2) = EE.*nbrev(:,:,2) - FF.*nbrev(:,:,1);   % nbrev_phi_thet_Contravariant
% nbrev_Contravariant = GG.*nbrev(:,:,3) - FF.*nbrev(:,:,4); % nbrev_thet_phi_Contravariant
% nbrev_phi_thet_Contravariant == nbrev_thet_phi_Contravariant
nbrev_Contravariant(:,:,3) = EE.*nbrev(:,:,4) - FF.*nbrev(:,:,3);   % nbrev_phi_phi_Contravariant

eta_thet = (gxi_thet.*nbrev_Contravariant(:,:,1) + ...
            gxi_phi.*nbrev_Contravariant(:,:,2))./JXibrev;
eta_phi = (gxi_thet.*nbrev_Contravariant(:,:,2) + ...
           gxi_phi.*nbrev_Contravariant(:,:,3))./JXibrev;
       
%%
f_thet = eta_thet + mu_theta - gamma_theta;
f_phi = eta_phi + mu_phi - gamma_phi;

[br,bi,~,~] = vhagcm(f_thet,f_phi);
p = -divgcm(br,bi);

%% Downsample
[ap, bp] = shagcm(p);
apdown = DownSampling(ap,size(axi,1)-1);
bpdown = DownSampling(bp,size(axi,1)-1);
p = shsgcm(apdown, bpdown);