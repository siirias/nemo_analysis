% Set paths for velocity and lon and lat
pathVel = [path,scenario,'/SS-GOB_1d_'];
pathNav = [path,scenario,'/SS-GOB_1m_'];

fileName = ['all_',scenario,'_',num2str(startyear),'_',num2str(endyear),'.mat']
navlat=ncread([pathNav,num2str(startyear),'0101_',num2str(startyear),'1231_ptrc_T.nc'],'nav_lat'); navlat=double(navlat); navlat(navlat==0)=NaN;
navlon=ncread([pathNav,num2str(startyear),'0101_',num2str(startyear),'1231_ptrc_T.nc'],'nav_lon'); navlon=double(navlon); navlon(navlon==0)=NaN;
% This is to ensure that there is no nan and that nalat anc navlon is
% covered as far aout as needed
navlatAll = zeros(size(navlat));
navlonAll = zeros(size(navlon));
navlatValues = max(navlat);
navlonValues = max(navlon,[],2);
for i=1:size(navlat,1)
    for j=1:size(navlat,2)
        navlatAll(i,j) = navlatValues(j);
        navlonAll(i,j) = navlonValues(i);
    end
end

% Read the mask
tmask=ncread('/nobackup/fouo6/sm_aho/smartsea/drivning/mesh_mask_GoB_1nm.nc','tmask'); tmask=double(tmask);

% Plot
% Define plot area
%addpath('m_map'); % use matlab m_map
 
% samfr 190327 since the mesh does not extend to the end west and southward
[Mlat,Ilat] = min(min(navlat,[],1)); [Mlon,Ilon] = min(min(navlon,[],2));
dlat = max(max(diff(navlat,1,2))); dlon = max(max(diff(navlon,1,1))); 
lat1 =  Mlat - (Ilat - 1)*dlat; lon1 = Mlon - (Ilon -1)*dlon;
lat = (lat1:dlat:lat1+(latsize)*dlat); lon = (lon1:dlon:lon1+(lonsize)*dlon);
lat2=max(max(lat)); lon2=max(max(lon));
latarea=[lat1 lat2]; lonarea=[lon1 lon2];
%m_proj('mercator','lon',lonarea,'lat',latarea);


% Remember that some years have 366 days. Plot each day but calculate mean
% and max as with only 28 days in February
% float uos(time_counter, depthu, y, x) ; in nc file
xDir = 340;
yDir = 360;
zDir = 36;

u = zeros(xDir, yDir, zDir, 366);
v = zeros(xDir, yDir, zDir, 366);
normVel = zeros(xDir, yDir, zDir, 366);
normVelMax = zeros(xDir,yDir,zDir,4);
meanCount = 1;
for year=startyear:endyear
    year
    u=ncread([pathVel,num2str(year),'0101_',num2str(year),'1231_grid_U.nc'],'uos');
    v=ncread([pathVel,num2str(year),'0101_',num2str(year),'1231_grid_V.nc'],'vos');
%end

% This might need to be inside the time loop as well

% Use the mask to put NaN outside the active cells
    for i=1:size(u,1)
        for j=1:size(u,2)
            for k=1:size(u,3)
                if(tmask(i,j,k)==0)
                    u(i,j,k,:)=NaN;
                    v(i,j,k,:)=NaN;
                end
            end
        end
    end

    normVel = sqrt(u.^2+v.^2);
    if(year==1960 || year==1964 || year==1968 || year==1972 || year==1976 || year==1980 || year==1984 || year==1988 || year==1992 ...
                        || year==1996 || year==2000 || year==2004 || year==2008 || year==2012 || year==2016 || year==2020 || year==2024 ...
                        || year==2028 || year==2032 || year==2036 || year==2040 || year==2044 || year==2048 || year==2052 || year==2056 ...
                        || year==2060 || year==2064 || year==2068 || year==2072 || year==2076 || year==2080 || year==2084 || year==2088 ...
                        || year==2092 || year==2096 || year==2100 || year==2104 || year==2108 || year==2112 || year==2116 || year==2120)
                    for day = 1:366
                        if meanCount == 1
                            M = squeeze(normVel(:,:,:,1));
                            Mold = zeros(size(M));
                            S = zeros(size(M));
                            SPrime2 = zeros(size(S));
                            meanCount = meanCount +1
                        else
                            Mold = M;
                            M = Mold + (squeeze(normVel(:,:,:,day))-Mold)./meanCount;
                            S = S + (squeeze(normVel(:,:,:,day))-Mold).*(squeeze(normVel(:,:,:,day))-M);
                            SPrime2 = S./(k-1);
                            meanCount = meanCount +1
                        end
                    end
    else
                    for day = 1:365
                        if meanCount == 1
                            M = squeeze(normVel(:,:,:,1));
                            Mold = zeros(size(M));
                            S = zeros(size(M));
                            SPrime2 = zeros(size(S));
                            meanCount = meanCount +1
                        else
                            Mold = M;
                            M = Mold + (squeeze(normVel(:,:,:,day))-Mold)./meanCount;
                            S = S + (squeeze(normVel(:,:,:,day))-Mold).*(squeeze(normVel(:,:,:,day))-M);
                            SPrime2 = S./(k-1);
                            meanCount = meanCount +1
                        end
                    end
    end

% Here it would be possible to calulate SPrime2 as well but is more
% efficient to do it in the end.

% For leap years the 28th and 29th of february will be saved on the same day
    for k = 1:size(normVelMax,3)
        k
        for j = 1:size(normVelMax,2)
            for i = 1:size(normVelMax,1)
                if ~isnan(normVel(i,j,k))% Detect and dont analyse nan cells
                    if(year==1960 || year==1964 || year==1968 || year==1972 || year==1976 || year==1980 || year==1984 || year==1988 || year==1992 ...
                        || year==1996 || year==2000 || year==2004 || year==2008 || year==2012 || year==2016 || year==2020 || year==2024 ...
                        || year==2028 || year==2032 || year==2036 || year==2040 || year==2044 || year==2048 || year==2052 || year==2056 ...
                        || year==2060 || year==2064 || year==2068 || year==2072 || year==2076 || year==2080 || year==2084 || year==2088 ...
                        || year==2092 || year==2096 || year==2100 || year==2104 || year==2108 || year==2112 || year==2116 || year==2120)
                        for day = 1:366
                        % day
                            if day < 60
                                if normVelMax(i,j,k,1) < normVel(i,j,k,day)
                                    normVelMax(i,j,k,1) = normVel(i,j,k,day);
                                    normVelMax(i,j,k,2) = atan2(v(i,j,k,day),u(i,j,k,day));
                                    normVelMax(i,j,k,3) = year;
                                    normVelMax(i,j,k,4) = day;
                                end
                            else
                                if normVelMax(i,j,k,1) < normVel(i,j,k,day)
                                    normVelMax(i,j,k,1) = normVel(i,j,k,day);
                                    normVelMax(i,j,k,2) = atan2(v(i,j,k,day),u(i,j,k,day));
                                    normVelMax(i,j,k,3) = year;
                                    normVelMax(i,j,k,4) = day-1;
                                end
                            end
                        end
                    else
                        for day = 1:365
                            %day
                            if normVelMax(i,j,k,1) < normVel(i,j,k,day)
                                normVelMax(i,j,k,1) = normVel(i,j,k,day);
                                normVelMax(i,j,k,2) = atan2(v(i,j,k,day),u(i,j,k,day));
                                normVelMax(i,j,k,3) = year;
                                normVelMax(i,j,k,4) = day;
                            end
                        end
                    end
                end
            end
        end  
    end


end

save(fileName,'-v7.3')
