close all; clear all;
load all_A001_1975_2005.mat
normVelA = zeros(340,360,36,6,5);
normVelA(:,:,:,1:4,1) = normVelMax;
normVelA(:,:,:,5,1) = M;
normVelA(:,:,:,6,1) = SPrime2;

% Use the mask
for i=1:size(normVelA,1)
    for j=1:size(normVelA,2)
        for k=1:size(normVelA,3)
            if(tmask(i,j,k)==0)
                normVelA(i,j,k,:,:)=NaN;
            end
        end
    end
end

save normVelA001_1975_2005.mat normVelA -v7.3

% Plot
% Define plot area
addpath('m_map'); % use matlab m_map

figure(1);
clf
set(gcf,'Color','w')
% This sets the units of the current figure (gcf = get current figure) on paper to centimeters.
set(gcf,'PaperUnits','centimeters')
% These are my size variables, width of 20 and a height of 20 for film in
% ppt
xSize = 20; ySize = 20;
% Additional coordinates to center the figure on A4-paper lying
xLeft = (30-xSize)/2; yTop = (21-ySize)/2;
% This command sets the position and size of the figure on the paper to the desired values
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
% This additional command resizes and repositions the figure window on the computer screen, so that it has
% the same aspect ratio as the final figure.
X = 100; Y = 50;
set(gcf,'Position',[X Y xSize*35 ySize*35])

width = 15/xSize;
height = 15/ySize;
px = 3.0/xSize;
py = 3.0/ySize;   % Distance to the first row
axes('position',[px py width height])

m_proj('azimuthal Equidistant','lon',20,'lat',63.3,'radius',5.','rectbox','on');
%m_pcolor(lon,lat,layers_allDepthMean_Phy3_4year'./Scaling(1,1)); shading flat;
m_pcolor(lon,lat,squeeze(normVelA(:,:,1,1,1))'); shading flat;
m_gshhs('ic','patch','w');
m_gshhs('ir','line','color',[0.7 0.7 0.7]);
axis xy;
m_grid('xtick',[15 18 21 24],'ytick',[60 62 64 66],'FontSize',20);
title(['A001 - Max surface speed [m s^-^1]'],'FontSize',16);
str = ['Year 19750101 - 20051231'];
m_text(16.,58.7,str,'FontSize',20)
colorbar('FontSize',12);
caxis([0 0.75]);
print('-dpng',['surfNormVelMax_A001_1975_2005.png'])
close all

figure(1);
clf
set(gcf,'Color','w')
% This sets the units of the current figure (gcf = get current figure) on paper to centimeters.
set(gcf,'PaperUnits','centimeters')
% These are my size variables, width of 20 and a height of 20 for film in
% ppt
xSize = 20; ySize = 20;
% Additional coordinates to center the figure on A4-paper lying
xLeft = (30-xSize)/2; yTop = (21-ySize)/2;
% This command sets the position and size of the figure on the paper to the desired values
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
% This additional command resizes and repositions the figure window on the computer screen, so that it has
% the same aspect ratio as the final figure.
X = 100; Y = 50;
set(gcf,'Position',[X Y xSize*35 ySize*35])

width = 15/xSize;
height = 15/ySize;
px = 3.0/xSize;
py = 3.0/ySize;   % Distance to the first row
axes('position',[px py width height])

m_proj('azimuthal Equidistant','lon',20,'lat',63.3,'radius',5.','rectbox','on');
%m_pcolor(lon,lat,layers_allDepthMean_Phy3_4year'./Scaling(1,1)); shading flat;
m_pcolor(lon,lat,squeeze(normVelA(:,:,1,5,1))'); shading flat;
m_gshhs('ic','patch','w');
m_gshhs('ir','line','color',[0.7 0.7 0.7]);
axis xy;
m_grid('xtick',[15 18 21 24],'ytick',[60 62 64 66],'FontSize',20);
title(['A001 - Mean surface speed [m s^-^1]'],'FontSize',16);
str = ['Year 19750101 - 20051231'];
m_text(16.,58.7,str,'FontSize',20)
colorbar('FontSize',12);
caxis([0 0.2]);
print('-dpng',['surfNormVelMean_A001_1975_2005.png'])
close all

figure(1);
clf
set(gcf,'Color','w')
% This sets the units of the current figure (gcf = get current figure) on paper to centimeters.
set(gcf,'PaperUnits','centimeters')
% These are my size variables, width of 20 and a height of 20 for film in
% ppt
xSize = 20; ySize = 20;
% Additional coordinates to center the figure on A4-paper lying
xLeft = (30-xSize)/2; yTop = (21-ySize)/2;
% This command sets the position and size of the figure on the paper to the desired values
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
% This additional command resizes and repositions the figure window on the computer screen, so that it has
% the same aspect ratio as the final figure.
X = 100; Y = 50;
set(gcf,'Position',[X Y xSize*35 ySize*35])

width = 15/xSize;
height = 15/ySize;
px = 3.0/xSize;
py = 3.0/ySize;   % Distance to the first row
axes('position',[px py width height])

m_proj('azimuthal Equidistant','lon',20,'lat',63.3,'radius',5.','rectbox','on');
%m_pcolor(lon,lat,layers_allDepthMean_Phy3_4year'./Scaling(1,1)); shading flat;
m_pcolor(lon,lat,squeeze(normVelA(:,:,1,6,1))'.^0.5); shading flat;
m_gshhs('ic','patch','w');
m_gshhs('ir','line','color',[0.7 0.7 0.7]);
axis xy;
m_grid('xtick',[15 18 21 24],'ytick',[60 62 64 66],'FontSize',20);
title(['A001 - Var^{0.5} surface speed [m s^-^1]'],'FontSize',16);
str = ['Year 19750101 - 20051231'];
m_text(16.,58.7,str,'FontSize',20)
colorbar('FontSize',12);
caxis([0 0.75]);
print('-dpng',['surfNormVelVarSquareRoot_A001_1975_2005.png'])
close all

figure(1);
clf
set(gcf,'Color','w')
% This sets the units of the current figure (gcf = get current figure) on paper to centimeters.
set(gcf,'PaperUnits','centimeters')
% These are my size variables, width of 20 and a height of 20 for film in
% ppt
xSize = 20; ySize = 20;
% Additional coordinates to center the figure on A4-paper lying
xLeft = (30-xSize)/2; yTop = (21-ySize)/2;
% This command sets the position and size of the figure on the paper to the desired values
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
% This additional command resizes and repositions the figure window on the computer screen, so that it has
% the same aspect ratio as the final figure.
X = 100; Y = 50;
set(gcf,'Position',[X Y xSize*35 ySize*35])

width = 15/xSize;
height = 15/ySize;
px = 3.0/xSize;
py = 3.0/ySize;   % Distance to the first row
axes('position',[px py width height])

m_proj('azimuthal Equidistant','lon',20,'lat',63.3,'radius',5.','rectbox','on');
%m_pcolor(lon,lat,layers_allDepthMean_Phy3_4year'./Scaling(1,1)); shading flat;
m_pcolor(lon,lat,squeeze(normVelA(:,:,1,3,1))'); shading flat;
m_gshhs('ic','patch','w');
m_gshhs('ir','line','color',[0.7 0.7 0.7]);
axis xy;
m_grid('xtick',[15 18 21 24],'ytick',[60 62 64 66],'FontSize',20);
title(['A001 - Year of max surface speed'],'FontSize',16);
str = ['Year 19750101 - 20051231'];
m_text(16.,58.7,str,'FontSize',20)
colorbar('Ticks',[1975,1980,1985,1990,1995,2000,2005],'FontSize',12);
caxis([1975 2005]);
print('-dpng',['surfNormVelMaxYear_A001_1975_2005.png'])
close all

figure(1);
clf
set(gcf,'Color','w')
% This sets the units of the current figure (gcf = get current figure) on paper to centimeters.
set(gcf,'PaperUnits','centimeters')
% These are my size variables, width of 20 and a height of 20 for film in
% ppt
xSize = 20; ySize = 20;
% Additional coordinates to center the figure on A4-paper lying
xLeft = (30-xSize)/2; yTop = (21-ySize)/2;
% This command sets the position and size of the figure on the paper to the desired values
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
% This additional command resizes and repositions the figure window on the computer screen, so that it has
% the same aspect ratio as the final figure.
X = 100; Y = 50;
set(gcf,'Position',[X Y xSize*35 ySize*35])

width = 15/xSize;
height = 15/ySize;
px = 3.0/xSize;
py = 3.0/ySize;   % Distance to the first row
axes('position',[px py width height])

m_proj('azimuthal Equidistant','lon',20,'lat',63.3,'radius',5.','rectbox','on');
%m_pcolor(lon,lat,layers_allDepthMean_Phy3_4year'./Scaling(1,1)); shading flat;
m_pcolor(lon,lat,squeeze(normVelA(:,:,1,4,1))'); shading flat;
m_gshhs('ic','patch','w');
m_gshhs('ir','line','color',[0.7 0.7 0.7]);
axis xy;
m_grid('xtick',[15 18 21 24],'ytick',[60 62 64 66],'FontSize',20);
title(['A001 - Day of max surface speed'],'FontSize',16);
str = ['Year 19750101 - 20051231'];
m_text(16.,58.7,str,'FontSize',20)
colorbar('FontSize',12);
caxis([1 365]);
print('-dpng',['surfNormVelMaxDay_A001_1975_2005.png'])
close all

