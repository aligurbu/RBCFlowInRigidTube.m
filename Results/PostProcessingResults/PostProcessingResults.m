%% Post-processing the flow around the RBC motion while flowing in the vessel
clear all; close all; clc;
addpath(genpath('../../../RBCFlowInRigidTube.m'))
verbose_Plot = false;

%% Input the model and parameters for the analysis from Models folder
LoadElasRBC_Short_Pr4_2_Time0_75s
% LoadElasRBC_RefCons_6mic_Pr8
% LoadElasRBC_LongConVes_Pr8

% LoadMemVisRBC_Short_muMem10_Pr4_2_Time0_75
% LoadMemVisRBC_RefCons_6mic_muMem_3_18_Pr8
% LoadMemVisRBC_LongConVes_muMem_3_18_P40

%% Set up output file
fidCoord = fopen(['Coord_',name,'.dat'],'r');
fidMemFor = fopen(['MemFor_',name,'.dat'],'r');
fidSol = fopen(['Sol_',name,'.dat'],'r');
fidTime = fopen(['Time_',name,'.dat'],'r');
Time = fread(fidTime,'double');
NSTEPS = length(Time);

%%
timeStepIncrement = 1;

%% Choose a view angle
viewInd = [0 90];

%% Choose a background color for visualization
blackBackground = true; % if false then white background

%% Visualization settings
VisualizeSettings
TransparencyInd = 0;

%% Set-up the evaluation points
SetUp_EvaluationPoints

%% Post-processing
PostProcessing_RBCInVessel

%%
fclose('all');