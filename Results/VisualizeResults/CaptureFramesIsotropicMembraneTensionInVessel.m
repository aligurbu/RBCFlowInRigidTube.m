%% Capture frames at visualize isotropic membrane tension of RBC
%% while flowing in a vessel

%% Create grid for representation: [N+1 x 2*N+1]
Nup = Rendering*N;
nlatup = Nup + 1;
nlonup = 2*Nup + 1;

%% Visualization settings
VisualizeSettings

%%
connect_rearrange = connect([1 5 2 6 3 7 4 8],:);

%% Patch faces (connectivity matrix)
faces = PatchFaces(nlatup, nlonup);

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
    if ~ismembertol(nstep,FRAMES)
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
T_step = T_step/RefShearRate; % in seconds
coord = coord*RefLength*10^(6); % \mum
CellX = CellX*RefLength*10^(6); % \mum
CellY = CellY*RefLength*10^(6); % \mum
CellZ = CellZ*RefLength*10^(6); % \mum
isotropic_nbrevVert = isotropic_nbrevVert*(RefElasticModulus*10^(6)); % \muN/m
Field = Field*(RefElasticModulus*10^(6)); % \muN/m

%%
for kk = 1:nframe
    fig = figure(kk);
    set(fig,'Color','white');
    hold on

    h = patch(CellX(:,:,kk),CellY(:,:,kk),CellZ(:,:,kk),Field(:,:,kk));
    set(h,'FaceColor','interp','EdgeColor','none', 'FaceAlpha', 0.6,...
          'EdgeAlpha','0.5', 'FaceLighting','gouraud',...
          'SpecularColorReflectance',0.1);
    h.CDataMapping = 'scaled';
    caxis([min(min(isotropic_nbrevVert)) max(max(isotropic_nbrevVert))])
    colormap(jet)
    material dull
    set(gca,'DataAspectRatio',[1 1 1])

    htube = patch('Faces',connect_rearrange','Vertices',coord');
    set(htube,'FaceColor','white','FaceAlpha',0.1,...
              'EdgeColor','none','EdgeAlpha','0.1','LineWidth',1.1,...
              'FaceLighting','gouraud',...
              'SpecularColorReflectance',0.1)
%     htube.BackFaceLighting = 'unlit';

    view(viewInd)
    if ~verbose_Plot
        ht = title({sprintf('Time = %4.4f sec',T_step(kk))},'Color','black');
        set(ht,'FontName','cambria math','FontSize',12)
        cbar = colorbar('southoutside');
        set(cbar,'FontSize',12,'Color','black')
        set(get(cbar,'title'),'string','\muN/m','FontSize',12,'Color','black');
        set(gca,'FontName','cambria math','FontSize',12)
    end
    axis equal
    axis off
    box off
    camlight('infinite')

    if verbose_Plot
        axis([-5 5 -5 5 -5 5])
        axis square
        set(gcf, 'InvertHardCopy', 'off');
        set(gcf,'PaperPositionMode','auto')
        exportgraphics(gca,['Frame_',sprintf('%d',kk), ...
                                 '_isotropicTensionInVessel',name,'.png']);
    end
end
if verbose_Plot
    fig = figure(kk+1);
    set(fig,'Color','white');
    hold on
    caxis([min(min(isotropic_nbrevVert)) max(max(isotropic_nbrevVert))])
    colormap(jet)
    cbar = colorbar('southoutside');
    set(cbar,'FontSize',12,'Color','black')
    set(get(cbar,'title'),'string','\muN/m','FontSize',12,'Color','black');
    set(gca,'FontName','cambria math','FontSize',12)
    view(viewInd)
    set(gca,'DataAspectRatio',[1 1 1])
    axis equal
    axis off
    box off
    set(gcf, 'InvertHardCopy', 'off');
    set(gcf,'PaperPositionMode','auto')
    exportgraphics(gca,['ColorBar_isotropicTensionInVessel',name,'.png']);
end