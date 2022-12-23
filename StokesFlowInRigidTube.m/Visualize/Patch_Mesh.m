function Patch_Mesh(coord, connect, Transparency)
%% 
% Patch plot the mesh from given coordinates of nodes and connectivity 
% matrix of elements. 

%%
if nargin > 2
    TransparencyInd = Transparency;
else
    TransparencyInd = 0.8;
end
%%
connect_rearrange = connect([1 5 2 6 3 7 4 8 1],:);

x = reshape(coord(1,connect_rearrange),size(connect_rearrange));
y = reshape(coord(2,connect_rearrange),size(connect_rearrange));
z = reshape(coord(3,connect_rearrange),size(connect_rearrange));

h = patch(x,y,z,'w');
alpha(h, TransparencyInd) % to set transparency
% set(h, 'linestyle', 'none') % no lines showning element edges
material dull
axis off
% xlabel('X'), ylabel('Y'), zlabel('Z')
view(3)
rotate3d
axis equal