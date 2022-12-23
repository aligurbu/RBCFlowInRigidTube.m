%% Analyze the results
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

%% Analyze the results of RBCInVessel
AnalyzeResults_RBCInVessel

%% Plots
%% Visualization settings
VisualizeSettings

%% Area Error and Volume Error
Area_Volume_Change

%% Forces and moments balance errors
Force_Moment_BalanceError

%% Three axes of ellipsoid with mass and volume of RBC
AxesOfEllipsoid


%%
fclose('all');