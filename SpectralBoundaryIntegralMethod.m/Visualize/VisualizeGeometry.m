function VisualizeGeometry(nlat, nlon, axi, bxi, ColorInd, PointInd, ...
                           Transparency)
%%
% Visualize the shape of an object for given spherical harmonics 
% coefficients axi, bxi
% Specify ColorInd for the patch color
% Set PointInd = true to plot the nodes of patches

%%
if nargin > 4
    PatchColor = ColorInd;
else
    PatchColor = 'r';
end

%% Visualization settings 
VisualizeSettings

%% Patch faces (connectivity matrix)
faces = PatchFaces(nlat, nlon);

%% Coordinates of vertices
% xiplot = shsecm(axi,bxi);
% Vert = reshape(xiplot,nlat*nlon,3)';
xi_equal = shsecm(axi,bxi);
xi_Gauss = shsgcm(axi,bxi);
xiplot = [xi_equal(1,:,:); xi_Gauss; xi_equal(end,:,:)];
Vert = reshape(xiplot,(nlat+2)*nlon,3)';

x = reshape(Vert(1,faces(:)),size(faces));
y = reshape(Vert(2,faces(:)),size(faces));
z = reshape(Vert(3,faces(:)),size(faces));

if nargin > 5 && PointInd
    plot3(x,y,z,'k.')%,'MarkerSize',MarkerSizeind)
end
if nargin > 6
    TransparencyInd = Transparency;
end
h = patch(x,y,z,PatchColor);
alpha(h, TransparencyInd) % to set transparency
% set(h, 'linestyle', 'none') % no lines showning element edges
axis off
material Dull
view(3)
rotate3d
axis equal