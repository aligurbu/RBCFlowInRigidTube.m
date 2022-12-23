function [nbrev, viscousStress_prev, epsilbrev_prev] = ...
                           MembraneTension(EE, FF, GG, WW, ...
                                           axi, bxi, ...
                                           UpSampleFactor, ...
                                           ES, ED, ...
                                           StVenantKirchhoff, ...
                                           NeoHookean, Skalak, ...
                                           MembraneViscoelasticity, ...
                                           mu_Mem, Tau, DT, ...
                                           viscousStress_prev, ...
                                           epsilbrev_prev)
%% Input
% EE, FF, GG = the first fundamental form of coefficients
% WW = surface area
% axi,bxi = spherical harmonic transform of deformed geometry coordinates
% upSampleFact = upsampling factor for dealising
% ES, ED = inplane-shear and inplane-dilatational moduli of membrane
% mu_Mem = membrane viscoelasticity
% Tau = the relaxation time of the Maxwell element
% DT = time increment
% viscoStress_prev = viscoelastic membrane stresses from previous step
% epsilbrev_prev = strain from previous step
% nbrev = in-plane membrane tensions consists of 
%         viscoelastic and elastic forces 
%%
%% Compute the first and second fundamental form coefficients
[~, ~, ee, ff, gg, ~, ~, ~, ~, ~, ~, ~, ~] = ...
                   coefficientsOfFundamentalForm(axi, bxi, UpSampleFactor);

%% Strain metric
epsilbrev = strainMetric(EE, GG, FF, WW, ee, gg, ff);

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

% %% Downsample
% [anbrev, bnbrev] = shagcm(nbrev);
% anbrev = DownSampling(anbrev,size(axi,1)-1);
% bnbrev = DownSampling(bnbrev,size(axi,1)-1);
% nbrev = shsgcm(anbrev, bnbrev);