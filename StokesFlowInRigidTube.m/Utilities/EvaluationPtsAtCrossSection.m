function [CrossSectionPts, ...
          CrossSectionPtsX, CrossSectionPtsY, CrossSectionPtsZ] = ...
          EvaluationPtsAtCrossSection(coord, wallnode, ...
                                             NumPointsInsideDomain, Offset)
if nargin > 3 
    OffsetPts = Offset;
else
    OffsetPts = 0;
end
%% Points on the cross section of vessel
CrossSectionPtsX_ = repmat(coord(1,wallnode), NumPointsInsideDomain, 1) - ...
                                                                 OffsetPts;
CrossSectionPtsY_ = coord(2,wallnode);
meanCrossSectionPtsY_ = mean(CrossSectionPtsY_);
CrossSectionPtsY_ = linspace(0.99,0,NumPointsInsideDomain)' * ...
                          (CrossSectionPtsY_ - meanCrossSectionPtsY_) + ...
                                                     meanCrossSectionPtsY_;
CrossSectionPtsZ_ = coord(3,wallnode);
meanCrossSectionPtsZ_ = mean(CrossSectionPtsZ_);
CrossSectionPtsZ_ = linspace(0.99,0,NumPointsInsideDomain)' * ...
                          (CrossSectionPtsZ_ - meanCrossSectionPtsZ_) + ...
                                                     meanCrossSectionPtsZ_;
CrossSectionPtsX = [CrossSectionPtsX_ CrossSectionPtsX_(:,1)];
CrossSectionPtsY = [CrossSectionPtsY_ CrossSectionPtsY_(:,1)];
CrossSectionPtsZ = [CrossSectionPtsZ_ CrossSectionPtsZ_(:,1)];
CrossSectionPts = zeros(size(CrossSectionPtsX,1), size(CrossSectionPtsX,2), 3);
CrossSectionPts(:,:,1) = CrossSectionPtsX;
CrossSectionPts(:,:,2) = CrossSectionPtsY;
CrossSectionPts(:,:,3) = CrossSectionPtsZ;