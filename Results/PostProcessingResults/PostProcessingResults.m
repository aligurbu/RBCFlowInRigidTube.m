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
% LoadMemVisRBC_LongConVes_muMem_3_18_Pr40

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
% viewInd = [0 90];
viewInd = [-45 30];

%% Choose a background color for visualization
blackBackground = true; % if false then white background

%% Visualization settings
VisualizeSettings
TransparencyInd = 0;

%% Set-up the evaluation points
SetUp_EvaluationPoints

%% Post-processing
%% Run once and load the results of post-processing.
%% Computing the velocity field at the evaluation points takes times
PostProcessing_RBCInVessel
return
if ~MembraneViscoelasticity
    if strcmp(nameVessel,'ShortMicrocapillary_16El')
        load PostProc_ElasRBC_Short_Pr4_2_Time0_75s.mat
    elseif strcmp(nameVessel,'RefCons_6mic_Ves_16El')
        load PostProc_ElasRBC_RefCons_6mic_Pr8.mat
    elseif strcmp(nameVessel,'LongConstrictedVessel_16El')
        load PostProc_ElasRBC_LongConVes_Pr8.mat
    end
else
    if strcmp(nameVessel,'ShortMicrocapillary_16El')
        load PostProc_MemVisRBC_Short_muMem10_Pr4_2_Time0_75.mat
    elseif strcmp(nameVessel,'RefCons_6mic_Ves_16El')
        load PostProc_MemVisRBC_RefCons_6mic_muMem_3_18_Pr8.mat
    elseif strcmp(nameVessel,'LongConstrictedVessel_16El')
        load PostProc_MemVisRBC_LongConVes_muMem_3_18_P40.mat
    end
end
Numframe = length(T_step);
%%
LineSpeed = zeros(size(LinePts,1), size(LinePts,2), Numframe);
MinLineSpeed = zeros(Numframe,1);
MaxLineSpeed = zeros(Numframe,1);
%%
CrossSectionSpeedmid = zeros(size(CirclePtsmid(:,:,1),1), ...
                             size(CirclePtsmid(:,:,1),2), Numframe);
MinCrossSectionSpeedmid = zeros(Numframe,1);
MaxCrossSectionSpeedmid = zeros(Numframe,1);
%%
CrossSectionSpeedinlet = zeros(size(CirclePtsinlet(:,:,1),1), ...
                               size(CirclePtsinlet(:,:,1),2), Numframe);
MinCrossSectionSpeedinlet = zeros(Numframe,1);
MaxCrossSectionSpeedinlet = zeros(Numframe,1);
%%
CrossSectionSpeedoutlet = zeros(size(CirclePtsoutlet(:,:,1),1), ...
                                size(CirclePtsoutlet(:,:,1),2), Numframe);
MinCrossSectionSpeedoutlet = zeros(Numframe,1);
MaxCrossSectionSpeedoutlet = zeros(Numframe,1);

%%
numelLinePts = numel(LinePts(:,:,1));

numelCirclePtsmid = numel(CirclePtsmid(:,:,1));
numelLinePtsCirclePtsmid = numelLinePts+numelCirclePtsmid;
rangeIndMiddle = numelLinePts+1 : numelLinePtsCirclePtsmid;

numelCirclePtsinlet = numel(CirclePtsinlet(:,:,1));
numelLinePtsCirclePtsmid_inlet = numelLinePtsCirclePtsmid + ...
                                 numelCirclePtsinlet;
rangeIndInlet = numelLinePtsCirclePtsmid+1 : numelLinePtsCirclePtsmid_inlet;

numelCirclePtsoutlet = numel(CirclePtsoutlet(:,:,1));
numelLinePtsCirclePtsmid_inlet_outlet = numelLinePtsCirclePtsmid_inlet + ...
                                        numelCirclePtsoutlet;
rangeIndOutlet = numelLinePtsCirclePtsmid_inlet+1 : ...
                                     numelLinePtsCirclePtsmid_inlet_outlet;

for nframe = 1:Numframe
    SpeedEvalPts = sqrt(VelEvalPtsX(nframe,:).^2 + ...
                        VelEvalPtsY(nframe,:).^2 + ...
                        VelEvalPtsZ(nframe,:).^2);

    %% Speed computation along the tube length
    LineSpeed_ = reshape(SpeedEvalPts(1:numelLinePts),size(LinePts,1), ...
                                                      size(LinePts,2));
    LineSpeed(:,:,nframe) = LineSpeed_;
    MinLineSpeed(nframe) = min(min(LineSpeed_));
    MaxLineSpeed(nframe) = max(max(LineSpeed_));

    %% Speed computation of middle cross-section of vessel
    CrossSectionSpeedmid_ = reshape(SpeedEvalPts(rangeIndMiddle), ...
                                    size(CirclePtsmid(:,:,1),1), ...
                                    size(CirclePtsmid(:,:,1),2));
    CrossSectionSpeedmid(:,:,nframe) = CrossSectionSpeedmid_;
    MinCrossSectionSpeedmid(nframe) = min(min(CrossSectionSpeedmid_));
    MaxCrossSectionSpeedmid(nframe) = max(max(CrossSectionSpeedmid_));

    %% Speed computation of inlet cross-section of vessel
    CrossSectionSpeedinlet_ = reshape(SpeedEvalPts(rangeIndInlet), ...
                                      size(CirclePtsinlet(:,:,1),1), ...
                                      size(CirclePtsinlet(:,:,1),2));
    CrossSectionSpeedinlet(:,:,nframe) = CrossSectionSpeedinlet_;
    MinCrossSectionSpeedinlet(nframe) = min(min(CrossSectionSpeedinlet_));
    MaxCrossSectionSpeedinlet(nframe) = max(max(CrossSectionSpeedinlet_));

    %% Speed computation of outlet cross-section of vessel
    CrossSectionSpeedoutlet_ = reshape(SpeedEvalPts(rangeIndOutlet), ...
                                       size(CirclePtsoutlet(:,:,1),1), ...
                                       size(CirclePtsoutlet(:,:,1),2));
    CrossSectionSpeedoutlet(:,:,nframe) = CrossSectionSpeedoutlet_;
    MinCrossSectionSpeedoutlet(nframe) = min(min(CrossSectionSpeedoutlet_));
    MaxCrossSectionSpeedoutlet(nframe) = max(max(CrossSectionSpeedoutlet_));

end

%% Dimensionalization
T_step = T_step/RefShearRate; % in seconds
CellX = CellX*(RefLength*10^(6)); % \mum
CellY = CellY*(RefLength*10^(6)); % \mum
CellZ = CellZ*(RefLength*10^(6)); % \mum
xi_GaussX = xi_GaussX*(RefLength*10^(6)); % \mum
xi_GaussY = xi_GaussY*(RefLength*10^(6)); % \mum
xi_GaussZ = xi_GaussZ*(RefLength*10^(6)); % \mum

u_GaussX = u_GaussX*RefVelocity*10^(3); %mm/sec
u_GaussY = u_GaussY*RefVelocity*10^(3); %mm/sec
u_GaussZ = u_GaussZ*RefVelocity*10^(3); %mm/sec

coord = coord*(RefLength*10^(6)); % \mum
EvaluationPts = EvaluationPts*(RefLength*10^(6)); % \mum
LinePts = LinePts*(RefLength*10^(6)); % \mum
CirclePtsmid = CirclePtsmid*(RefLength*10^(6)); % \mum
CirclePtsinlet = CirclePtsinlet*(RefLength*10^(6)); % \mum
CirclePtsoutlet = CirclePtsoutlet*(RefLength*10^(6)); % \mum

VelEvalPtsX = VelEvalPtsX*RefVelocity*10^(3); %mm/sec
VelEvalPtsY = VelEvalPtsY*RefVelocity*10^(3); %mm/sec
VelEvalPtsZ = VelEvalPtsZ*RefVelocity*10^(3); %mm/sec
LineSpeed = LineSpeed*RefVelocity*10^(3); %mm/sec
CrossSectionSpeedmid = CrossSectionSpeedmid*RefVelocity*10^(3); %mm/sec
CrossSectionSpeedinlet = CrossSectionSpeedinlet*RefVelocity*10^(3); %mm/sec
CrossSectionSpeedoutlet = CrossSectionSpeedoutlet*RefVelocity*10^(3); %mm/sec

%% Contour graph of speed of the fluid inside the vessel
MinCaxis = min([MinLineSpeed; MinCrossSectionSpeedmid; ...
                              MinCrossSectionSpeedinlet; ...
                              MinCrossSectionSpeedoutlet]);
MaxCaxis = max([MaxLineSpeed; MaxCrossSectionSpeedmid; ...
                              MaxCrossSectionSpeedinlet; ...
                              MaxCrossSectionSpeedoutlet]);
MinCaxis = MinCaxis*RefVelocity*10^(3); %mm/sec
MaxCaxis = MaxCaxis*RefVelocity*10^(3); %mm/sec

%%
if blackBackground
    figure('Color','black');
    ColorInd = 'white';
else
    figure('Color','white');
    ColorInd = 'black';
end
hold on

h = patch(CellX(:,:,1),CellY(:,:,1),CellZ(:,:,1),'r');
alpha(h, TransparencyInd) % to set transparency
material Dull

set(gca,'DataAspectRatio',[1 1 1])
ht = title({sprintf('Time = %4.4f sec',T_step(1))},'Color',ColorInd);
set(ht,'FontName','cambria math','FontSize',12)

Patch_Mesh(coord, connect, 0.1)

% LineSpeed
[M,cLineSpeed] = contourf(LinePts(:,:,1), LinePts(:,:,2), ...
                          LineSpeed(:,:,1), 100, ...
                          'LineColor', 'none');

% Middle section speed
ax = gca;
HGCirclePtsmid = hgtransform(ax);
[M,cCirclePtsmid] = contourf(CirclePtsmid(:,:,2), CirclePtsmid(:,:,3), ...
                             CrossSectionSpeedmid(:,:,1), 100, ...
                             'LineColor', 'none', 'Parent', HGCirclePtsmid);
HGCirclePtsmid.Matrix = makehgtform('translate',CirclePtsmid(1,1,1),0,0)* ...
                        makehgtform('yrotate', pi/2)* ...
                        makehgtform('zrotate', pi/2);

% inlet section speed
ax = gca;
HGCirclePtsinlet = hgtransform(ax);
[M,cCirclePtsinlet] = contourf(CirclePtsinlet(:,:,2), CirclePtsinlet(:,:,3), ...
                               CrossSectionSpeedinlet(:,:,1), 100, ...
                               'LineColor', 'none', 'Parent', HGCirclePtsinlet);
HGCirclePtsinlet.Matrix = makehgtform('translate',CirclePtsinlet(1,1,1),0,0)* ...
                          makehgtform('yrotate', pi/2)* ...
                          makehgtform('zrotate', pi/2);

% outlet section speed
ax = gca;
HGCirclePtsoutlet = hgtransform(ax);
[M,cCirclePtsoutlet] = contourf(CirclePtsoutlet(:,:,2), CirclePtsoutlet(:,:,3), ...
                                CrossSectionSpeedoutlet(:,:,1), 100, ...
                                'LineColor', 'none', 'Parent', HGCirclePtsoutlet);
HGCirclePtsoutlet.Matrix = makehgtform('translate',CirclePtsoutlet(1,1,1),0,0)* ...
                           makehgtform('yrotate', pi/2)* ...
                           makehgtform('zrotate', pi/2);

caxis([MinCaxis MaxCaxis])
colormap(jet)
cbar = colorbar('southoutside');
set(cbar,'FontSize',12,'Color',ColorInd)
set(get(cbar,'title'),'string','Speed (mm/s)','Color',ColorInd);
set(gca,'FontName','cambria math','FontSize',12)
view(viewInd)
axis equal
axis off
box off
rotate3d

%
v = VideoWriter(['PostProcessing',name,'.mp4'],'MPEG-4');
v.FrameRate = 10;
open(v);

% Loop through by changing XData and YData
for id = 1:length(T_step)
    %% Update graphics data. This is more efficient than recreating plots.
    set(h, 'XData', CellX(:,:,id), ...
           'YData', CellY(:,:,id), ...
           'ZData', CellZ(:,:,id), ...
           'FaceAlpha', TransparencyInd);
    set(ht, 'String', {sprintf('Time = %4.4f sec',T_step(id))})
    set(cLineSpeed, 'ZData', LineSpeed(:,:,id));
    set(cCirclePtsmid, 'ZData', CrossSectionSpeedmid(:,:,id));
    set(cCirclePtsinlet, 'ZData', CrossSectionSpeedinlet(:,:,id));
    set(cCirclePtsoutlet, 'ZData', CrossSectionSpeedoutlet(:,:,id));
    %% Get frame as an image
    frame = getframe(gcf);
    writeVideo(v,frame);
end
close(v);

fclose('all');