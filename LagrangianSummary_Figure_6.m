%%%Plot Particle Distribution After 12 or 24 Hours
%Used for Plotting Figure 6
%------------------------------------------------------------------------
%Author: R Alan Mason
%Created: 10/25/25 
%Edited: 10/25/25
%% ----------------------------------------------------------------------

% - Import Data

%Length of Hours to Plot
timelength = 12;

%Use Spring Tide or Neap Tide Data?
springtide = 0;
neaptide = 1;

if springtide == 1      %Spring Tide - 1 km
    %Full Flow
    load('\SpringTidePartTrack_FullFlow\PartTrackData.mat')
elseif neaptide == 1    %Neap Tide - 1 km
    %Full Flow
    load('\NeapTidePartTrack_FullFlow\PartTrackData.mat');
end

%Set-Up Needed Parameters
tims = 1:6:858;         %Time Start    
parts = 1:100:2400;     %First index point of a particle release
partf = 100:100:2400;   %Last index point of a particle release 
if timelength == 12
    timf = 145:6:858;
elseif timelength == 24
    timf = 289:6:858;
end
%Set-Up Parameters for ending points of particles
finalposone = zeros(24,100);
for ii = 1:24
    finalposone(ii,:) = squeeze(ys(timf(ii),parts(ii):partf(ii)));
end
finalpos = cat(1,finalposone,finalposone(1,:));

%% - Plot Particle Locations at Time Denoted by "timelength"

pxs = ones(25,100).*[0:24]';

figure
plot(pxs,finalpos,'k.')
xticks(0:3:24)
xticklabels({'0','1.5','3','4.5','6','7.5','9','10.5','12'})
xlabel('Release Time (hour)')
ylim([0 1])
ylabel('End Cross-Channel Location (km)')
if timelength == 12
    title('12 hours')
elseif timelength == 24
    title('24 hours')
end

%% - Plot Particle End Point Concentrations

%Run Histogram to bin data
figure
cmap = cmocean('dense',20);
%histvec = 0.0025:0.005:1.005;
histvec = 0.005:0.01:1.005;
%histvec = 0.025:0.05:1;
parthist = zeros(size(finalpos,1),(size(histvec,2)-1));

for iii = 1:size(finalpos,1)
    datasplit = squeeze(finalpos(iii,:));
    h = histogram(datasplit,histvec);
    parthist(iii,:) = h.Values;
end

%Set up needed arrays to plot histogram data
pxpc = ones(size(parthist,1),size(parthist,2)).*[0:24]';
ydif = (diff(histvec)/2)+histvec(1:(length(histvec)-1));
pypc = (ones(size(parthist,1),size(parthist,2)).*ydif)-(ydif(1,1)/2);
pypcnew = (ones(size(parthist,1),size(parthist,2)).*ydif);
partpercentorg = parthist/100;
partpercent = partpercentorg;
partpercent(partpercent == 0) = nan;

%Plot histogram data with ending points for each release time
figure
pcolor(pxpc,pypc,partpercent)
shading('flat')
xticks(0:3:24)
xticklabels({'0','1.5','3','4.5','6','7.5','9','10.5','12'})
xlabel('Release Time (hour)')
ylim([0 1])
ylabel('End Cross-Channel Location (km)')
if timelength == 12
    title('12 hours')
elseif timelength == 24
    title('24 hours')
end
c1 = colorbar;
colormap(cmap)
caxis([0 1])
c1.Label.String = 'Particle Fraction';         %Colorbar label
grid on

%% - Plot Average Position of Particles

mpp = mean(squeeze(partpercentorg(1:24,:)),1);

figure
plot(pypc(1,:),mpp,'k-')
xlim([0 1])
xlabel('Cross-Channel End Location (km)')
ylim([0 1])
ylabel('Particle Percentage')
if timelength == 12
    title('12 hours')
elseif timelength == 24
    title('24 hours')
end
