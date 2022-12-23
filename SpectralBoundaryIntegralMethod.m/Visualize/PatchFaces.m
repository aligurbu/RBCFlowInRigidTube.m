function faces = PatchFaces(nlat, nlon)
%% Patch faces (connectivity matrix)
faces = zeros(4,(nlat+1)*nlon);
% do all but the last column of patches
for lon = 1:nlon-1
    for lat = 1:nlat+1
        %(lat,lon),(lat+1,lon),(lat+1,lon+1),(lat,lon+1)
        faces(:,(lon-1)*(nlat+1)+lat) = [(lon-1)*(nlat+2)+lat ...
                                         (lon-1)*(nlat+2)+lat+1 ...
                                          lon*(nlat+2)+lat+1 ...
                                          lon*(nlat+2)+lat]; 
    end
end
% do the last column of patches
for lat = 1:nlat+1
    %(lat,nlon),(lat+1,nlon),(lat+1,1),(lat,1)
    faces(:,(nlon-1)*(nlat+1)+lat) = [(nlon-1)*(nlat+2)+lat ...
                                      (nlon-1)*(nlat+2)+lat+1 ...
                                       lat+1 ...
                                       lat];
end