%% Set-up the evaluation points
%% A line along the tube length
numNodesOnCircumference = 2*16; % for 16 elements around the circimference
if strcmp(nameVessel,'ShortMicrocapillary_16El')
    numNodesAlongTubeLength = 2*16+1; % for 16 element along the length
elseif strcmp(nameVessel,'RefCons_6mic_Ves_16El')
    numNodesAlongTubeLength = 2*48+1; % for 48 element along the length
elseif strcmp(nameVessel,'LongConstrictedVessel_16El')
    numNodesAlongTubeLength = 2*48+1; % for 48 element along the length
end
wallnodeArray = reshape(wallnode, numNodesOnCircumference, ...
                                  numNodesAlongTubeLength);
if strcmp(nameVessel,'ShortMicrocapillary_16El')
    wallnodeOnlyWallLineUp = wallnodeArray(1,[3:numNodesAlongTubeLength 2]);
    wallnodeOnlyWallLineDown = wallnodeArray(17,[3:numNodesAlongTubeLength 2]);
    OffsetLinePtsX = TubeLength/(2*(numNodesAlongTubeLength-1));
elseif strcmp(nameVessel,'RefCons_6mic_Ves_16El')
    wallnodeOnlyWallLineUp = wallnodeArray(1,3:numNodesAlongTubeLength);
    wallnodeOnlyWallLineDown = wallnodeArray(17,3:numNodesAlongTubeLength);
    OffsetLinePtsX = 0;
elseif strcmp(nameVessel,'LongConstrictedVessel_16El')
    wallnodeOnlyWallLineUp = wallnodeArray(1,3:numNodesAlongTubeLength);
    wallnodeOnlyWallLineDown = wallnodeArray(17,3:numNodesAlongTubeLength);
    OffsetLinePtsX = 0;
end

%% Point grids along the vessel length; LinePts
NumPointsInsideDomainCrossSection = 21;
%% Points along the vessel length
LinePtsX = repmat(coord(1,wallnodeOnlyWallLineUp), ...
                       NumPointsInsideDomainCrossSection,1)-OffsetLinePtsX;
LinePtsY_ = repmat(coord(2,wallnodeOnlyWallLineUp), ...
                        NumPointsInsideDomainCrossSection,1);
LinePtsZ = repmat(coord(3,wallnodeOnlyWallLineUp), ...
                         NumPointsInsideDomainCrossSection,1);
DiffLinePtsY_ = (coord(2,wallnodeOnlyWallLineUp) - ...
                    coord(2,wallnodeOnlyWallLineDown));
LinePtsY = LinePtsY_ + ...
          linspace(-0.01,-0.99,NumPointsInsideDomainCrossSection)' * ...
                                                             DiffLinePtsY_;
% Note that LinePtsY points start 1% from the vessel wall
LinePts = zeros(size(LinePtsX,1), size(LinePtsX,2), 3);
LinePts(:,:,1) = LinePtsX; 
LinePts(:,:,2) = LinePtsY;
LinePts(:,:,3) = LinePtsZ;

%% Circular point grid in middle of the vessel
NumPointsInsideDomainOnRadius = 10;
if strcmp(nameVessel,'ShortMicrocapillary_16El')
    wallnodeWallCircle = wallnodeArray(:,18);
elseif strcmp(nameVessel,'RefCons_6mic_Ves_16El')
    wallnodeWallCircle = wallnodeArray(:,50);
elseif strcmp(nameVessel,'LongConstrictedVessel_16El')
    wallnodeWallCircle = wallnodeArray(:,50);
end
[CirclePtsmid, CirclePtsXmid, CirclePtsYmid, CirclePtsZmid] = ...
                EvaluationPtsAtCrossSection(coord, wallnodeWallCircle, ...
                                            NumPointsInsideDomainOnRadius);

%% Circular point grid in vicinity of inlet surface
wallnodeWallCircle = wallnodeArray(:,3);
[CirclePtsinlet, CirclePtsXinlet, CirclePtsYinlet, CirclePtsZinlet] = ...
                EvaluationPtsAtCrossSection(coord, wallnodeWallCircle, ...
                                            NumPointsInsideDomainOnRadius, ...
                                            OffsetLinePtsX);

%% Circular point grid in vicinity of outlet surface
wallnodeWallCircle = wallnodeArray(:,end);
[CirclePtsoutlet, CirclePtsXoutlet, CirclePtsYoutlet, CirclePtsZoutlet] = ...
                EvaluationPtsAtCrossSection(coord, wallnodeWallCircle, ...
                                            NumPointsInsideDomainOnRadius, ...
                                            -OffsetLinePtsX);
%% Put the points together
EvaluationPts = zeros(3, numel(LinePts(:,:,1)) + ...
                         numel(CirclePtsXmid) + ...
                         numel(CirclePtsXinlet) + ...
                         numel(CirclePtsXoutlet));
EvaluationPts(1,:) = [LinePtsX(:); CirclePtsXmid(:); ...
                                   CirclePtsXinlet(:); ...
                                   CirclePtsXoutlet(:)];
EvaluationPts(2,:) = [LinePtsY(:); CirclePtsYmid(:); ...
                                   CirclePtsYinlet(:); ...
                                   CirclePtsYoutlet(:)];
EvaluationPts(3,:) = [LinePtsZ(:); CirclePtsZmid(:); ...
                                   CirclePtsZinlet(:); ...
                                   CirclePtsZoutlet(:)];

numPts = size(EvaluationPts,2); % Total number of points

%% Depict EvaluationPts on vessel
if verbose_Plot
    figure('Color','white')
    hold on
    Patch_Mesh(coord, connect, TransparencyInd)
    plot3(EvaluationPts(1,:), EvaluationPts(2,:), EvaluationPts(3,:), ...
          'k.','MarkerSize',MarkerSizeind)
end