%Code that Plots Particle Track Data over the Surface Velocity Field
%Used for Figures 5, 7, 8, and 9 
%------------------------------------------------------------------------
%Author: R Alan Mason
%Created: 10/25/25 
%Edited: 10/25/25
%% ----------------------------------------------------------------------

tic

%Set Flags
IncremMajor = 1;    %Plot High Slack, Max Ebb, Low Slack, Max Flood
IncremPartial = 0;  %Plot every 1.5 hour release
showall = 1;        %Show all releases against the Eulerian Mean Flow
minuslang = 1;      %Plot Lagrangian Residual Flow (Full Tidal minus Eulerian Mean)
isv = 0;            %Save Plots?

%% - Load Data

%Load Model Data for needed plotting variables
folderpath = '\ParticleTrack\1000mWidth\';
load([folderpath,'ModelNumerics\it6000.mat']);

%Set-Up Needed Values for Particle Subset Plotting
if IncremMajor == 1
    pcut = [1:100;601:700;1201:1300;1801:1900];
    tst = [1,181,361,541,721];
    starttime = 1:36:145;
    endtime = starttime+576;
elseif IncremPartial == 1
    pcut = [1:100;301:400;601:700;901:1000;1201:1300;1501:1600;1801:1900;2101:2200];
    starttime = 1:18:127;
    endtime = starttime+144;
end

%Download Particles Transported by the Eulerian Mean Current
load([folderpath,'PartTrack_EulerianAvg\PartTrackData.mat'])
xsea = xs;              ysea = ys;
xsea(xsea==0)=nan;      ysea(ysea==0)=nan;

%Download Particles transported by either the Full Flow field or Lagrangian
%Residual Flow depending on flags set
if minuslang == 1
    load([folderpath,'PartTrack_MinusEulAvg\PartTrackData.mat']);
else
    load([folderpath,'PartTrack_FullFlow\PartTrackData.mat']);
end
xs(xs==0)=nan;          ys(ys==0)=nan;

%% - Plot Full Path

%Set-Up needed Plotting Arrays
distarrho = ones(size(u,1),1).*(((0:(ny-1)).*dy)');
plotwidth = ones(size(vsurf,1),151).*(1:size(vsurf,1))';
plottime = ones(size(vsurf,1),151).*(distarrho'./1000);
plottimep = ones(size(ys,1),size(ys,2)).*(1:(1/6):(143+(5/6)))';
plottimef = ones(size(ys,1),size(ys,2)).*(1:(1/6):(143+(5/6)))';

%Find particle position after one tidal cycle
lowval = (distarrho(60)/1000);   %74
highval = (distarrho(85)/1000);  %78
for iphase = 1:5
    if iphase == 5
        phaseend = endtime(1);
        ycut(iphase,:) = ysea(phaseend,pcut(1,:));
    else
        phaseend = endtime(iphase);
        ycut(iphase,:) = ys(phaseend,pcut(iphase,:));
    end
    for np = 1:length(ycut)
        unoparto = ycut(iphase,np);
        if unoparto >= lowval && unoparto <= highval
            conpart(iphase,np) = 1;
        else
            conpart(iphase,np) = 0;
        end
    end
end

%Find percent of particles that converge to the channel center
for ip = 1:5
    percentcon(ip,1) = sum(conpart(ip,:))/100;
end

%Create velocity field for plotting for just Lagrangian Residual Flow
if minuslang == 1
    vsurf = vsurf-vsurfavg;
end

%Plot
figure
set(gcf,'Position',[150   150   1400   800])
%Particles with Eulerian Mean at the End
subplot(2,1,1)
pcolor(plotwidth',plottime',vsurf');    shading('flat');
hold on
if IncremMajor == 1
    plot(plottimep(tst(1):6:end,pcut(1,:)),ys(starttime(1):6:end,pcut(1,:)),'k-');
    plot(plottimep(tst(2):6:end,pcut(2,:)),ys(starttime(2):6:709,pcut(2,:)),'r-');
    plot(plottimep(tst(3):6:end,pcut(3,:)),ys(starttime(3):6:565,pcut(3,:)),'b-');
    plot(plottimep(tst(4):6:end,pcut(4,:)),ys(starttime(4):6:421,pcut(4,:)),'g-');
    plot(plottimep(tst(5):6:end,pcut(1,:)),ysea(starttime(1):6:133,pcut(1,:),:),'Color',[0.9290 0.6940 0.1250]);
elseif IncremPartial == 1
    plot(plottimep(starttime(1):6:end,pcut(1,:)),ys(starttime(1):6:end,pcut(1,:)),'k-');
    plot(plottimep(starttime(2):6:end,pcut(2,:)),ys(starttime(2):6:end,pcut(2,:)),'k-.');
    plot(plottimep(starttime(3):6:end,pcut(3,:)),ys(starttime(3):6:end,pcut(3,:)),'r-');
    plot(plottimep(starttime(4):6:end,pcut(4,:)),ys(starttime(4):6:end,pcut(4,:)),'r-.');
    plot(plottimep(starttime(5):6:end,pcut(5,:)),ys(starttime(5):6:end,pcut(5,:)),'b-');
    plot(plottimep(starttime(6):6:end,pcut(6,:)),ys(starttime(6):6:end,pcut(6,:)),'b-.');
    plot(plottimep(starttime(7):6:end,pcut(7,:)),ys(starttime(7):6:end,pcut(7,:)),'g-');
    plot(plottimep(starttime(8):6:end,pcut(8,:)),ys(starttime(8):6:end,pcut(8,:)),'g-.');
end
for ip2 = 1:5
    for numpart = 1:length(ycut)
        conyn = conpart(ip2,numpart);
        if ip2 == 5
            if conyn == 1
                plot(plottimep(tst(ip2),1),ysea(starttime(1),pcut(1,numpart)),'rs','MarkerFaceColor','k');
            end
        else
            if conyn == 1
                plot(plottimep(tst(ip2),1),ys(starttime(ip2),pcut(ip2,numpart)),'rs','MarkerFaceColor','k');
            end
        end
    end
end
xlabel('Time (hours)')
xticks(1:12:145)
xticklabels({'0','6','12','18','24','30','36','42','48','54','60','66','72'})
ylabel('Width (km)')
title('Release Time Comparison - Across-Channel Velocity')
c1 = colorbar;
cmap = cmocean('balance',16);
colormap(cmap)
caxis([-0.075 0.075])
c1.Label.String = 'Across-Channel Velocity (m s^{-1})';

%Particles with Eulerian Mean overlaid
subplot(2,1,2)
%set(gcf,'Position',[150   150   1400   400])
pcolor(plotwidth',plottime',vsurf');    shading('flat');
hold on
if showall == 1
    plot(plottimep(tst(1):6:end,pcut(1,:)),ys(starttime(1):6:end,pcut(1,:)),'k-');
    plot(plottimep(tst(1):6:end,pcut(1,:)),ysea(starttime(1):6:end,pcut(1,:),:),'Color',[0.9290 0.6940 0.1250]);
    plot(plottimep(tst(2):6:end,pcut(2,:)),ys(starttime(2):6:709,pcut(2,:)),'r-');
    plot(plottimep(tst(2):6:end,pcut(1,:)),ysea(starttime(2):6:709,pcut(2,:),:),'Color',[0.9290 0.6940 0.1250]);
    plot(plottimep(tst(3):6:end,pcut(3,:)),ys(starttime(3):6:565,pcut(3,:)),'b-');
    plot(plottimep(tst(3):6:end,pcut(3,:)),ysea(starttime(3):6:565,pcut(3,:),:),'Color',[0.9290 0.6940 0.1250]);
    plot(plottimep(tst(4):6:end,pcut(4,:)),ys(starttime(4):6:421,pcut(4,:)),'g-');
    plot(plottimep(tst(4):6:end,pcut(4,:)),ysea(starttime(4):6:421,pcut(4,:),:),'Color',[0.9290 0.6940 0.1250]);
else
    plot(plottimep(starttime(1):6:end,pcut(1,:)),ys(starttime(1):6:end,pcut(1,:),:),'k-');
    plot(plottimep(starttime(1):6:end,pcut(1,:)),ysea(starttime(1):6:end,pcut(1,:),:),'r-');
end
xlabel('Time (hours)')
xticks(1:12:145)
xticklabels({'0','6','12','18','24','30','36','42','48','54','60','66','72'})
ylabel('Width (km)')
title('Average vs Full Flow Release - Across-Channel Velocity')
c1 = colorbar;
cmap = cmocean('balance',16);
colormap(cmap)
caxis([-0.075 0.075])
c1.Label.String = 'Across-Channel Velocity (m s^{-1})';

if isv == 1
    if minuslang == 1
        print([folderpath,'PartTrackComparisons_LangMinusEul'],'-dpng');
    else
        print([folderpath,'PartTrackComparisons'],'-dpng');
    end
end

%% - Paper Plots
if minuslang == 1 %Code for Panels B and D for Figures 5, 7, 8 and 9
    plottimemid = (plottime-0.5)*1000;
    ysmid = (ys-0.5)*1000;
    yseamid = (ysea-0.5)*1000;
    %plottimemid = (plottime-2.5)*1000;
    %ysmid = (ys-2.5)*1000;
    %yseamid = (ysea-2.5)*1000;
    vsurfnew = vsurf;
    vsurfnew(121:144,:)= vsurfavg.*ones(24,1);
    figure
    set(gcf,'Position',[150   150   1400   400])
    pcolor(plotwidth',plottimemid',vsurfnew');    shading('flat');
    hold on
    plot(plottimep(tst(1):6:end,pcut(1,:)),ysmid(starttime(1):6:end,pcut(1,:)),'k-');
    plot(plottimep(tst(2):6:end,pcut(2,:)),ysmid(starttime(2):6:709,pcut(2,:)),'r-');
    plot(plottimep(tst(3):6:end,pcut(3,:)),ysmid(starttime(3):6:565,pcut(3,:)),'b-');
    plot(plottimep(tst(4):6:end,pcut(4,:)),ysmid(starttime(4):6:421,pcut(4,:)),'g-');
    plot(plottimep(tst(5):6:end,pcut(1,:)),yseamid(starttime(1):6:133,pcut(1,:),:),'Color',[0.9290 0.6940 0.1250]);
    for ip2 = 1:4
        for numpart = 1:length(ycut)
            conyn = conpart(ip2,numpart);
            if ip2 == 5
                if conyn == 1
                    plot(plottimep(tst(ip2),1),yseamid(starttime(1),pcut(1,numpart)),'rs','MarkerFaceColor','k');
                end
            else
                if conyn == 1
                    plot(plottimep(tst(ip2),1),ysmid(starttime(ip2),pcut(ip2,numpart)),'rs','MarkerFaceColor','k');
                end
            end
        end
    end
    xlabel('Time (hours)','FontSize',16)
    xticks(1:12:145)
    xticklabels({'0','6','12','18','24','30','36','42','48','54','60','66','72'})
    xlim([1 145])
    ylabel('Width (m)','FontSize',16)
    ylim([-500 500])
    %ylim([-2500 2500])
    title('Release Time Comparison - Across-Channel Velocity','FontSize',16)
    c1 = colorbar;
    %cmap = cmocean('balance',16);
    cmap = cmocean('balance');
    colormap(cmap)
    caxis([-0.05 0.05])
    %caxis([-0.065 0.065])
    c1.Label.String = 'Across-Channel Velocity (m s^{-1})';
    ax = gca; 
    ax.FontSize = 16;
else %Code for Panels A and C for Figures 5, 7, 8 and 9
    plottimemid = (plottime-0.5)*1000;
    ysmid = (ys-0.5)*1000;
    yseamid = (ysea-0.5)*1000;
    %plottimemid = (plottime-2.5)*1000;
    %ysmid = (ys-2.5)*1000;
    %yseamid = (ysea-2.5)*1000;
    figure
    set(gcf,'Position',[150   150   1400   400])
    pcolor(plotwidth',plottimemid',vsurf');    shading('flat');
    hold on
    if IncremMajor == 1
        plot(plottimep(tst(1):6:end,pcut(1,:)),ysmid(starttime(1):6:end,pcut(1,:)),'k-');
        plot(plottimep(tst(2):6:end,pcut(2,:)),ysmid(starttime(2):6:709,pcut(2,:)),'r-');
        plot(plottimep(tst(3):6:end,pcut(3,:)),ysmid(starttime(3):6:565,pcut(3,:)),'b-');
        plot(plottimep(tst(4):6:end,pcut(4,:)),ysmid(starttime(4):6:421,pcut(4,:)),'g-');
        %plot(plottimep(tst(5):6:end,pcut(1,:)),yseamid(starttime(1):6:133,pcut(1,:),:),'Color',[0.9290 0.6940 0.1250]);
    elseif IncremPartial == 1
        plot(plottimep(starttime(1):6:end,pcut(1,:)),ys(starttime(1):6:end,pcut(1,:)),'k-');
        plot(plottimep(starttime(2):6:end,pcut(2,:)),ys(starttime(2):6:end,pcut(2,:)),'k-.');
        plot(plottimep(starttime(3):6:end,pcut(3,:)),ys(starttime(3):6:end,pcut(3,:)),'r-');
        plot(plottimep(starttime(4):6:end,pcut(4,:)),ys(starttime(4):6:end,pcut(4,:)),'r-.');
        plot(plottimep(starttime(5):6:end,pcut(5,:)),ys(starttime(5):6:end,pcut(5,:)),'b-');
        plot(plottimep(starttime(6):6:end,pcut(6,:)),ys(starttime(6):6:end,pcut(6,:)),'b-.');
        plot(plottimep(starttime(7):6:end,pcut(7,:)),ys(starttime(7):6:end,pcut(7,:)),'g-');
        plot(plottimep(starttime(8):6:end,pcut(8,:)),ys(starttime(8):6:end,pcut(8,:)),'g-.');
    end
    for ip2 = 1:4
        for numpart = 1:length(ycut)
            conyn = conpart(ip2,numpart);
            if ip2 == 5
                if conyn == 1
                    plot(plottimep(tst(ip2),1),yseamid(starttime(1),pcut(1,numpart)),'rs','MarkerFaceColor','k');
                end
            else
                if conyn == 1
                    plot(plottimep(tst(ip2),1),ysmid(starttime(ip2),pcut(ip2,numpart)),'rs','MarkerFaceColor','k');
                end
            end
        end
    end
    xlabel('Time (hours)','FontSize',16)
    xticks(1:12:145)
    xticklabels({'0','6','12','18','24','30','36','42','48','54','60','66','72'})
    xlim([1 121])
    %xlim([1 145])
    ylabel('Width (m)','FontSize',16)
    ylim([-500 500])
    %ylim([-2500 2500])
    title('Release Time Comparison - Across-Channel Velocity','FontSize',16)
    c1 = colorbar;
    cmap = cmocean('balance');
    %cmap = cmocean('balance',16);
    colormap(cmap)
    caxis([-0.05 0.05])
    %caxis([-0.065 0.065])
    c1.Label.String = 'Across-Channel Velocity (m s^{-1})';
    ax = gca; 
    ax.FontSize = 16;
end
if isv == 1
    if minuslang == 1
        print([folderpath,'PartTrack_PaperFig_LangMinusEul'],'-dpng');
    else
        print([folderpath,'PartTrack_PaperFig'],'-dpng');
    end
end

toc

