%%%Idealized Model Plot u, v, and s Snapshots
%Used for Figures 2, 3, 10, and 11
%------------------------------------------------------------------------
%Author: R Alan Mason
%Created: 10/25/25 
%Edited: 10/25/25
%% ----------------------------------------------------------------------

% - Set Flags and Parameters
timesteps = 108:0.5:120;
isv = 0;
filepath = '\FilePath\SpringTide\';
filename = 'ModelNumerics\it';

%% - Load Data

%Set-Up Blank Arrays
maxu = zeros(length(timesteps),1);
maxv = zeros(length(timesteps),1);
maxs = zeros(length(timesteps),1);
surfu = zeros(length(timesteps),151);
depthavgu = zeros(length(timesteps),151);
ct = 1000;

%Plot u, v, and s data for each time step
for i = 1:length(timesteps)

    tval = timesteps(i).*36000;
    load([filepath,filename,num2str(tval),'.mat']);

    %Set-Up Needed Arrays for Plotting
    halfwidth = B/2;
    depthplotr = ones(size(u,1),size(u,2)).*((-(nz-1):0).*dz);
    widthplotr = ones(size(u,1),size(u,2)).*(((0:(ny-1)).*dy)');
    depthplotv = ones(size(v,1),size(v,2)).*((-(nz-1):0).*dz);
    widthplotv = ones(size(v,1),size(v,2)).*(((0:(ny)).*dy)');
    depthplotw = ones(size(w,1),size(w,2)).*((-(nz):0).*dz);
    widthplotw = ones(size(w,1),size(w,2)).*(((0:(ny-1)).*dy)');
    widthplotrn = widthplotr-halfwidth;
    widthplotvn = widthplotv-halfwidth;
    widthplotwn = widthplotw-halfwidth;

    figure
    set(gcf,'Position',[150   50   750   900])
    cmap = cmocean('balance');
    %Plot u data
    subplot(3,1,1);
    pcolor(widthplotrn',depthplotr',-u'); shading('flat')
    maxu(i,1) = max(max(abs(u)));
    xlim([-halfwidth halfwidth])
    %xticks([-500 -400 -300 -200 -100 0 100 200 300 400 500])
    %xticks([-2000 -1500 -1000 -500 0 500 1000 1500 2000])
    xticks([-5000 -4000 -3000 -2000 -1000 0 1000 2000 3000 4000 5000])
    xlabel('Distance (m)')
    ylabel('Depth (m)')
    c1 = colorbar;
    colormap(cmap)
    caxis([-1 1]);
    c1.Label.String = 'Along-Channel Velocity (m s^{-1})';
    hold on
    conval = -1:0.1:1;
    %conval = -1:0.01:1;
    contour(widthplotrn',depthplotr',u',conval,'k-','ShowText','on')
    if length(timesteps)==4
        if i == 1
            title('High Slack')
        elseif i == 2
            title('Max Ebb')
        elseif i == 3
            title('Low Slack')
        elseif i == 4
            title('Max Flood')
        end
    else
        title(['Hour ',num2str(timesteps(i))])
    end
    ax = gca; 
    ax.FontSize = 14;

    surfu(i,:) = u(:,76);
    bottomu(i,:) = u(76,1);
    depthavgu(i,:) = mean(u,2,'omitnan');
    salttime(i,:,:) = s;
    
    %Plot v data
    subplot(3,1,2);
    pcolor(widthplotvn',depthplotv',v'); shading('flat')
    maxv(i,1) = max(max(abs(v)));
    xlim([-halfwidth halfwidth])
    %xticks([-500 -400 -300 -200 -100 0 100 200 300 400 500])
    %xticks([-2000 -1500 -1000 -500 0 500 1000 1500 2000])
    xticks([-5000 -4000 -3000 -2000 -1000 0 1000 2000 3000 4000 5000])
    xlabel('Distance (m)')
    ylabel('Depth (m)')
    c1 = colorbar;
    %caxis([-0.0045 0.0045]);
    %caxis([-0.018 0.018]);
    caxis([-0.1 0.1]);
    c1.Label.String = 'Across-Channel Velocity (m s^{-1})';
    %conval = -1:0.005:1;
    conval = -1:0.01:1;
    hold on
    contour(widthplotvn',depthplotv',v',conval,'k-','ShowText','on')
    ax = gca; 
    ax.FontSize = 14;

    %Plot s data
    subplot(3,1,3);
    %subplot(4,1,4);
    pcolor(widthplotrn',depthplotr',s'); shading('flat')
    maxs(i,1) = max(max(abs(s)));
    xlim([-halfwidth halfwidth])
    %xticks([-500 -400 -300 -200 -100 0 100 200 300 400 500])
    %xticks([-2000 -1500 -1000 -500 0 500 1000 1500 2000])
    xticks([-5000 -4000 -3000 -2000 -1000 0 1000 2000 3000 4000 5000])
    xlabel('Distance (m)')
    ylabel('Depth (m)')
    c1 = colorbar;
    caxis([-2 2]);
    c1.Label.String = 'Salinity (PSU)';
    conval = -2:0.1:2;
    hold on
    contour(widthplotrn',depthplotr',s',conval,'k-','ShowText','on')
    ax = gca; 
    ax.FontSize = 14;

    if isv == 1
        print([filepath,'TidalPhaseFigs\PhaseFig',num2str(ct+i)],'-dpng');
        close all
    end

    %Plot Center u-profile
    figure
    plot(-u(76,:)',depthplotr')
    xlim([-1 1])
    xlabel('Along-Channel Velocity (m s^{-1}')
    ylabel('Depth (m)')
    title(['Hour ',num2str(timesteps(i))])
    grid on
    if isv == 1
        print([filepath,'TidalPhaseFigs\CenterProfile',num2str(ct+i)],'-dpng');
        close all
    end

    close all

end


%% - Kukulka and Chant 2025 Parameters

kpgamma = ((h0^2)/(B^2))*(Kmy/Km);
advect = ((g^2)*(bet^2)*(sx^2)*(h0^10))/((Km^4)*(B^2));
