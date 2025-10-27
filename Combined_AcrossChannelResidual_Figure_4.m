%%%Code to Plot the Surface Velocity Field of the Across-Channel Velocities
%%Can Plot both Full Flow Field and Combined Tidaly Varying/Mean Field
%%Used for Figure 4
%------------------------------------------------------------------------
%Author: R Alan Mason
%Created: 10/24/25 
%Edited: 10/24/25
%------------------------------------------------------------------------

%Initial Flags
isv=0;              %Save Figures?
spring = 1;         %Plot Spring Tide (1)? Or Neap Tide (0)?
resid = 1;          %Plot Tidal Varying Flow (1)? Or Full Flow Field (0)?

%Load In Data
if spring == 1
    load('.\BaseCases\HighVerticalEddyViscosity\FullFlowField_PartTrackData.mat');
else
    load('.\BaseCases\LowVerticalEddyViscosity\FullFlowField_PartTrackData.mat');
end
%%
%Create Variables Needed for Plotting
plotwidth = ((ones(size(vsurf,1),151).*(1:size(vsurf,1))'));
distarrho = ones(151,1).*(((0:(151-1)).*6.622516556291391)');
plottime = ones(size(vsurf,1),151).*(distarrho')-500;   %Correct width so that center of channel = 0
vaary = ones(108,151).*vsurfavg;    
vsurfnew = vsurf-vsurfavg;          %Solve for Tidally Varying Portion
vsurfnew(37:144,:) = vaary;         %Stitch Tidally Varing and Tidal Mean

%Plot Data
figure
set(gcf,'Position',[150   150   1000   400])
if resid == 1
    pcolor(plotwidth',plottime',vsurfnew');    shading('flat');
else
    pcolor(plotwidth',plottime',vsurf');    shading('flat');
end
hold on
%plot(plottimep(1:6:end,:),ys(1:6:end,:),'k-');
xlim([1 43])
xticks([1:6:145])
xticklabels({'0','3','6','9','12','15','18','21','24','27','30','33','36','39','42','45','48','51','54','57','60','63','66','69','72'})
xlabel('Time (hours)')
ylim([-500 500])
ylabel('Across-Channel Distance (km)')
%title('Residual Velocity')
c1 = colorbar;
cmap = cmocean('balance');
colormap(cmap)
%caxis([-0.125 0.125])
caxis([-0.06 0.06])
c1.Label.String = 'Across-Channel Velocity (m s^{-1})';
if isv == 1
    if resid == 1
        if spring == 1
            print('.\SpringTide\AcrossChannelResidualVel','-dpng');
        else
            print('.\NeapTide\\AcrossChannelResidualVel','-dpng');
        end
    else
        if spring == 1
            print('.\SpringTide\AcrossChannelFullVel','-dpng');
        else
            print('.\NeapTide\AcrossChannelFullVel','-dpng');
        end
    end
end
