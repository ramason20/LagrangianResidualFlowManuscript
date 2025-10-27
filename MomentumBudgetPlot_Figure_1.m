%%%Idealized Model Momentum Budget Plot
%Plot Code for Figure 1
%------------------------------------------------------------------------
%Author: R Alan Mason
%Created: 10/24/25 
%Edited: 10/24/25
%% ------------------------------------------------------------------------
tic
% - Flags and Parameters

dt = 0.1;               %Model Time Step
min15res = 0;           %Is the data saved in 15minute resolution
min10res = 1;           %Is the data saved in 15minute resolution
starthr = 96;          %Starting hour for the time series
endhr = 120;            %Ending hour for the time series
isv = 0;                %Save plots (1=Y 0=N)
cori = 0;               %Coriolis Data? (1=Y 0=N)
datacheck = 0;          %Extra plot code to ensure the data looks right
cnt = 1000;             %Save numbering variable
savdir = '\FolderLocation\';
fpath = [savdir,'FileNameMinusTimeStep'];    %Base to files for the for loop below
savdir = [savdir,'MomentumBudget\'];
timeplot = 0;

%Global Variables
depthlev = 70;  %Depth Level to Depth --> 70=~1m 65=~2m

%Set-Up File Name Time Array
if min15res == 1
    hrts = 3600/dt;
    starttime = starthr*hrts;   %25 hrs = 120000;  37 hrs = 177600;  48 hrs = 230400;   61 hrs = 292800;  108 hrs = 518400; 
    endtime = endhr*hrts;     %36 hrs = 172800;  48 hrs = 230400;  60 hours = 288000; 72 hrs = 345600;  120 hrs = 576000;
    timesteps = starttime:(hrts/4):endtime;
elseif min10res == 1
    hrts = 3600/dt;
    starttime = starthr*hrts;   %25 hrs = 120000;  37 hrs = 177600;  48 hrs = 230400;   61 hrs = 292800;  108 hrs = 518400; 
    endtime = endhr*hrts;     %36 hrs = 172800;  48 hrs = 230400;  60 hours = 288000; 72 hrs = 345600;  120 hrs = 576000;
    timesteps = starttime:(hrts/6):endtime;
else
    hrts = 3600/dt;
    starttime = starthr*hrts;   %25 hrs = 120000;  37 hrs = 177600;  48 hrs = 230400;   61 hrs = 292800;  108 hrs = 518400; 
    endtime = endhr*hrts;     %36 hrs = 172800;  48 hrs = 230400;  60 hours = 288000; 72 hrs = 345600;  120 hrs = 576000;
    timesteps = starttime:hrts:endtime;
end

%% - Plot Code Lerczak & Geyer 2004 Plots

%Load first time step to set-up preallocated arrays
ffp = [fpath,num2str(timesteps(1)),'.mat'];
load(ffp);

%Set-Up the Arrays
ufull = zeros(length(timesteps),size(u,1),size(u,2));
uaccelfull = zeros(length(timesteps),size(u,1),size(u,2));
uptropfull = zeros(length(timesteps),size(u,1),size(u,2));
upclinfull = zeros(length(timesteps),size(u,1),size(u,2));
uavecyfull = zeros(length(timesteps),size(u,1),size(u,2));
uaveczfull = zeros(length(timesteps),size(u,1),size(u,2));
uturbyfull = zeros(length(timesteps),size(u,1),size(u,2));
uturbzfull = zeros(length(timesteps),size(u,1),size(u,2));

for itloc = 1:length(timesteps)

    %Load in the data
    ts = timesteps(itloc);
    ffp = [fpath,num2str(ts),'.mat'];
    load(ffp);

    %Concatenate to 1 array
    ufull(itloc,:,:) = u;
    uaccelfull(itloc,:,:) = uaccel;
    uptropfull(itloc,:,:) = uptrop+utide;
    upclinfull(itloc,:,:) = upclin;
    uavecyfull(itloc,:,:) = uavecy;
    uaveczfull(itloc,:,:) = uavecz;
    uturbyfull(itloc,:,:) = uturby;
    uturbzfull(itloc,:,:) = uturbz;
end

%%
% Vertical Budget - Center Channel
crossloc = 54;      %One side = 54; Middle = 76
depthar = -h0:dz:0;     depthar = depthar(1,2:end);

uamean = squeeze(uaccelfull(:,crossloc,:));
    uamean = mean(uamean,1,'omitnan');
upgmean = squeeze(uptropfull(:,crossloc,:))+squeeze(upclinfull(:,crossloc,:));
    upgmean = mean(upgmean,1,'omitnan');
uavmean = squeeze(uavecyfull(:,crossloc,:))+squeeze(uaveczfull(:,crossloc,:));
    uavmean = mean(uavmean,1,'omitnan');
utbmean = squeeze(uturbyfull(:,crossloc,:))+squeeze(uturbzfull(:,crossloc,:));
    utbmean = mean(utbmean,1,'omitnan');

figure
plot(-uamean,depthar,'k-')
hold on
plot(upgmean,depthar,'r-')
plot(utbmean,depthar,'g-')
plot(uavmean,depthar,'b-')
legend('Acceleration','Pressure Gradient','Turbulence','Advection','Location','northwest')
title('Time Averaged Momentum Budget Profile')
ylabel('Depth (m)')
xlabel('Acceleration (m s^{-2})')
xlim([-0.00001 0.00001])
grid on
if isv == 1
    print([savdir,'PaperPlot_OffCenter'],'-dpng');
end

toc