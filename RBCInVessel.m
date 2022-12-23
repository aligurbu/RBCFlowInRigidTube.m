%% Red blood cell flow inside a vessel
% Spectral Galerkin BIE for an RBC motion and deformation and
% direct BEM for the rigid vessel implementations are put together
% to analyze the stresses develop on the cell membrane while flowing
% in a vessel with a size comparable to red cell or while squeezing
% through a constriction.
%%
clear all; close all; clc;
addpath(genpath('../RBCFlowInRigidTube.m'))

%% Input the model and parameters for the analysis from Models folder
LoadElasRBC_Short_Pr4_2_Time0_75s

%% Vessel set-up
ParametersForTheAnalysis_Vessel

%% RBC set-up
ParametersForTheAnalysis_RBC