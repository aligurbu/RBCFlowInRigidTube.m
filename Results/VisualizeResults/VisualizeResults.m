%% Visualize the results
clear all; close all; clc;
addpath(genpath('../../../RBCFlowInRigidTube.m'))
verbose_Plot = false;
WritetoGIF = true; % if false then save as MPEG-4

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
timeStepIncrement = 4;
FRAMES = floor(linspace(1,NSTEPS,5));
TimesOfFrames = Time(FRAMES)/RefShearRate % in seconds

Rendering = 2; % This upsampling for the rendering purposes

%% Choose a view angle
% viewInd = [-45 30];
% name = [name,'_3D'];
viewInd = [0 90];
name = [name,'_xy'];

%% Choose a background color for visualization
blackBackground = false; % if false then white background

%% Run Visualize m-files one at a time.

%% Visualization of RBC while flowing in a vessel

VisualizeMembraneShapeInVessel

% CaptureFramesMembraneShapeInVessel

% VisualizeIsotropicMembraneTensionInVessel

% CaptureFramesIsotropicMembraneTensionInVessel

fclose('all');