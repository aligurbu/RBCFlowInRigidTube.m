%% Capture frames at visualization of
%% membrane shape while flowing in a vessel

%% Create grid for representation: [N+1 x 2*N+1]
nlat = N+1;
nlon = 2*N+1;

%% Visualization settings
VisualizeSettings

%%
connect_rearrange = connect([1 5 2 6 3 7 4 8],:);

%% Patch faces (connectivity matrix)
faces = PatchFaces(nlat, nlon);

%% Masks to go between Spherepack and vector representations of SH coeff
mask_a = repmat(triu(true(N+1),0),1,1,3);
mask_b = mask_a;
mask_b(1,:,:) = false;

%%
Numframe = 0;
for nstep = 1:NSTEPS
    if ~ismembertol(nstep,FRAMES)
        continue
    end
    Numframe = Numframe + 1;
end

%%
Vertices = zeros((nlat+2)*nlon,3,Numframe);

nframe = 0;
T_step = zeros(Numframe,1);
for nstep = 1:NSTEPS
    %% Read from file
    cxi = fread(fidCoord,3*(N+1)^2,'double');
    axi = zeros(size(mask_a));  bxi = zeros(size(mask_b));
    axi(mask_a) = cxi(1:3*(N+1)*(N+2)/2);
    bxi(mask_b) = cxi(3*(N+1)*(N+2)/2+1:3*(N+1)^2);

    if ~ismembertol(nstep,FRAMES)
        continue
    end
    nframe = nframe + 1;
    T_step(nframe) = Time(nstep);

    %% Coordinates of vertices
    xi_equal = shsecm(axi,bxi);
    xi_Gauss = shsgcm(axi,bxi);
    xiplot = [xi_equal(1,:,:); xi_Gauss; xi_equal(end,:,:)];
    Vertices(:,:,nframe) = reshape(xiplot,(nlat+2)*nlon,3);
end
%% Dimensionalization
T_step = T_step/RefShearRate; % in seconds
coord = coord*RefLength*10^(6); % \mum
Vertices = Vertices*RefLength*10^(6); % \mum
%%
for kk = 1:nframe
    fig = figure(kk);
    set(fig,'Color','white');
    ColorInd = 'black';
    hold on

    htube = patch('Faces',connect_rearrange','Vertices',coord');
    set(htube,'FaceColor','white','FaceAlpha',0.1,...
              'EdgeColor','none','EdgeAlpha','0.1','LineWidth',1.1,...
              'FaceLighting','gouraud',...
              'SpecularColorReflectance',0.1)
%     htube.BackFaceLighting = 'unlit';

    hRBC = patch('Faces',faces','Vertices',Vertices(:,:,kk),...
                 'FaceColor','r','FaceAlpha',0.6,...
                 'EdgeColor','none','EdgeAlpha','0.5', 'FaceLighting','gouraud',...
                 'SpecularColorReflectance',0.1);
    material Dull
    if ~verbose_Plot
        ht = title({sprintf('Time = %4.4f sec',T_step(kk))},'Color',ColorInd);
        set(ht,'FontName','cambria math','FontSize',12)
    end
    set(gca,'DataAspectRatio',[1 1 1])
    material dull
    view(viewInd)
    axis equal
    axis off
    box off
    camlight('infinite')

    if verbose_Plot
        axis([-5 5 -5 5 -5 5])
        axis square
        set(gcf, 'InvertHardCopy', 'off');
        set(gcf,'Color','white');
        set(gcf,'PaperPositionMode','auto')
        exportgraphics(gca,['Frame_',sprintf('%d',kk), ...
                            '_MembraneProfile',name,'.png'], ...
                            'ContentType','image');%,'Resolution',600);
    end
end