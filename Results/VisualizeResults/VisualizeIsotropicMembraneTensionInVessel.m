%% Visualize isotropic membrane tension of RBC while flowing in a vessel
%% Create grid for representation: [N+1 x 2*N+1]
Nup = Rendering*N;
nlatup = Nup + 1;
nlonup = 2*Nup + 1;

numDofPerNode = size(coord,1); % Number of velocity components

%% Visualization settings 
VisualizeSettings

%% Patch faces (connectivity matrix)
faces = PatchFaces(nlatup, nlonup);

%% Masks to go between Spherepack and vector representations of SH coeff
mask_a = repmat(triu(true(N+1),0),1,1,3);
mask_b = mask_a;
mask_b(1,:,:) = false;

%%
Numframe = 0;
for nstep = 1:NSTEPS
    if (nstep~=1 && nstep~=NSTEPS && mod(nstep,timeStepIncrement)~=0)
        continue
    end
    Numframe = Numframe + 1;
end

%%
viscousStress_prev = zeros(UpSampleFactor*N+1, 2*UpSampleFactor*N+1, 4);
epsilbrev_prev = zeros(UpSampleFactor*N+1, 2*UpSampleFactor*N+1, 4);

%%
CellX = zeros(size(faces,1),size(faces,2),Numframe);
CellY = zeros(size(faces,1),size(faces,2),Numframe);
CellZ = zeros(size(faces,1),size(faces,2),Numframe);
isotropic_nbrevVert = zeros((nlatup+2)*nlonup,Numframe);
Field = zeros(size(faces,1),size(faces,2),Numframe);

nframe = 0;
T_step = zeros(Numframe,1);
for nstep = 1:NSTEPS
    %% Read from file
    cxi = fread(fidCoord,3*(N+1)^2,'double');
    axi = zeros(size(mask_a));  bxi = zeros(size(mask_b));
    axi(mask_a) = cxi(1:3*(N+1)*(N+2)/2);
    bxi(mask_b) = cxi(3*(N+1)*(N+2)/2+1:3*(N+1)^2);
    
    if nstep == 1
        %% Compute the first and second fundamental form coefficients
        [~, ~, EE, FF, GG, WW, ~, ~, ~, ~, ~, ~, ~] = ...
                   coefficientsOfFundamentalForm(axi, bxi, Rendering);
    end
    if (nstep~=1 && nstep~=NSTEPS && mod(nstep,timeStepIncrement)~=0)
        continue
    end
    nframe = nframe + 1;
    T_step(nframe) = Time(nstep);
    
    %% Compute in-plane membrane tensions
    [nbrev, viscousStress_prev, epsilbrev_prev] = ...
                           MembraneTension(EE, FF, GG, WW, ...
                                           axi, bxi, ...
                                           UpSampleFactor, ...
                                           ES, ED, ...
                                           StVenantKirchhoff, ...
                                           NeoHookean, Skalak, ...
                                           MembraneViscoelasticity, ...
                                           mu_Mem, Tau, DT, ...
                                           viscousStress_prev, ...
                                           epsilbrev_prev);
    isotropic_nbrev = zeros(size(nbrev,1),size(nbrev,2));
    for kk = 1:size(nbrev,1)
        for mm = 1:size(nbrev,2)
            eig_nbrev = eig(reshape(nbrev(kk,mm,:),2,2));
            isotropic_nbrev(kk,mm) = 0.5*sum(eig_nbrev);
        end
    end

    axiup = UpSampling(axi,N,Nup);
    bxiup = UpSampling(bxi,N,Nup);
    %% Coordinates of vertices
    xi_equal = shsecm(axiup,bxiup);
    xi_Gauss = shsgcm(axiup,bxiup);
    xiplot = [xi_equal(1,:,:); xi_Gauss; xi_equal(end,:,:)];
    Vert = reshape(xiplot,(nlatup+2)*nlonup,3);
    CellX(:,:,nframe) = reshape(Vert(faces(:),1),size(faces));
    CellY(:,:,nframe) = reshape(Vert(faces(:),2),size(faces));
    CellZ(:,:,nframe) = reshape(Vert(faces(:),3),size(faces));
    %% Fields values on patches
    [aisotropic_nbrev,bisotropic_nbrev] = shagcm(isotropic_nbrev); 
    isotropic_nbrev_equal = shsecm(aisotropic_nbrev,bisotropic_nbrev);
    isotropic_nbrev_Gauss = shsgcm(aisotropic_nbrev,bisotropic_nbrev);
    isotropic_nbrevplot = [isotropic_nbrev_equal(1,:,:); ...
                           isotropic_nbrev_Gauss; ...
                           isotropic_nbrev_equal(end,:,:)];
    isotropic_nbrevVert(:,nframe) = ...
                          reshape(isotropic_nbrevplot,(nlatup+2)*nlonup,1);
    Field(:,:,nframe) = ...
                reshape(isotropic_nbrevVert(faces(:),nframe), size(faces));
end
%% Dimensionalization 
T_step = (T_step*1000)/RefShearRate; % in milli seconds
coord = coord*RefLength*10^(6); % \mum
CellX = CellX*RefLength*10^(6); % \mum
CellY = CellY*RefLength*10^(6); % \mum
CellZ = CellZ*RefLength*10^(6); % \mum
isotropic_nbrevVert = isotropic_nbrevVert*(RefElasticModulus*10^(6)); % \muN/m
Field = Field*(RefElasticModulus*10^(6)); % \muN/m

%%
figure('Color','white');
ColorInd = 'black';
Patch_Mesh(coord, connect, 0.1)
hold on

h = patch(CellX(:,:,1),CellY(:,:,1),CellZ(:,:,1),Field(:,:,1));
alpha(h, TransparencyInd) % to set transparency
set(h,'FaceColor','interp','EdgeColor','none')
h.CDataMapping = 'scaled';
caxis([min(min(isotropic_nbrevVert)) max(max(isotropic_nbrevVert))])
colormap(jet)
cbar = colorbar('southoutside');
set(cbar,'FontSize',12,'Color',ColorInd)
set(get(cbar,'title'),'string','\muN/m','FontSize',12,'Color',ColorInd);
set(gca,'FontName','cambria math','FontSize',12)
material Dull
axis off
view(viewInd)
set(gca,'DataAspectRatio',[1 1 1])
ht = title({sprintf('Time = %4.4f ms',T_step(1))},'Color',ColorInd);
set(ht,'FontName','cambria math','FontSize',12)
axis equal
axis off
box off

%% Get figure size
if WritetoGIF
    frame = getframe(gcf);
    height = size(frame.cdata, 1);
    width = size(frame.cdata, 2);

    %% Preallocate data (for storing frame data)
    mov = zeros(height, width, 1, length(T_step), 'uint8');

    %% Loop through by changing XData and YData
    for id = 1:length(T_step)
        %% Update graphics data. This is more efficient than recreating plots.
        set(h, 'XData', CellX(:,:,id), ...
               'YData', CellY(:,:,id), ...
               'ZData', CellZ(:,:,id), ...
               'CData', Field(:,:,id), ...
               'FaceAlpha', TransparencyInd, ...
               'FaceColor','interp', 'EdgeColor','none')
        set(ht, 'String', {sprintf('Time = %4.4f ms',T_step(id))})
        %% Get frame as an image
        frame = getframe(gcf);

        %% Create a colormap for the first frame. For the rest of the frames,
        %% use the same colormap
        if id == 1
            [mov(:,:,1,id), map] = rgb2ind(frame.cdata, 1024, 'dither');
        else
            mov(:,:,1,id) = rgb2ind(frame.cdata, map, 'dither');
        end
    end
    %% Create animated GIF
    imwrite(mov, map, ['isotropicTensionInVessel',name,'.gif'], 'DelayTime', 0, 'LoopCount', inf)
else % Write to MP4 file 
    v = VideoWriter(['isotropicTensionInVessel',name,'.mp4'],'MPEG-4');
    v.FrameRate = 10;
    open(v);

    %% Loop through by changing XData and YData
    for id = 1:length(T_step)
        %% Update graphics data. This is more efficient than recreating plots.
        set(h, 'XData', CellX(:,:,id), ...
               'YData', CellY(:,:,id), ...
               'ZData', CellZ(:,:,id), ...
               'CData', Field(:,:,id), ...
               'FaceAlpha', TransparencyInd, ...
               'FaceColor','interp', 'EdgeColor','none')
        set(ht, 'String', {sprintf('Time = %4.4f ms',T_step(id))})
        
        %% Get frame as an image
        frame = getframe(gcf);
        writeVideo(v,frame);
    end
    close(v);
end