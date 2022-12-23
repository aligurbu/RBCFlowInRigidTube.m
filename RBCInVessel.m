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