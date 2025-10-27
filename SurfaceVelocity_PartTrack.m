%%%Idealized Model Particle Track
%------------------------------------------------------------------------
%Author: R Alan Mason
%Created: 10/25/25 
%Edited: 10/25/25
%% ----------------------------------------------------------------------

% - Set Parameters
numpart = 100;          %Number of particles to be released
timestep = 300;         %Particle time step in seconds
isv = 0;                %Save Plot? 1=Yes; 0=No;
isd = 0;                %Save Data? 1=Yes; 0=No;
euleravg = 0;           %Use an Eulerian Average Flow Field
lmeavg = 0;             %Use full flow field minus the Eulerian Average
numtidcyc = 6;          %Number of Tidal cycles to use (1 is just 12 hours)
halfhour = 1;           %Use 30 minute model output (halfhour=1) or full hour (halfhour=0)
min15out = 0;           %Use 15 minute model output
addpart = 1;

%File Path Set-Up
savloc = '\FilePath\SpringTide\';
fpath = [savloc,'ModelNumerics\it'];    %Base to files for the for loop below
casename = 'Spring Tide Base Case';

%Set which Hours of the Model Output to Track Particles Over
if halfhour == 1
    if addpart == 1
        %hourrange = (108:0.5:119.5).*36000;
        hourrange = (96:0.5:119.5).*36000;
        %hourrange = (72:0.5:119.5).*36000;
        %hourrange = (228:0.5:239.5).*36000;
    else
        hourrange = (108:0.5:119.5).*36000;  %For 10km Estuary - High Slack
        %hourrange = (72:0.5:119.5).*36000;
        %hourrange = (105:0.5:116.5).*36000;  %For 10km Estuary - Max Flood
        %hourrange = (102:0.5:113.5).*36000;  %For 10km Estuary - Low Slack
        %hourrange = (99:0.5:110.5).*36000;  %For 10km Estuary - Max Ebb
    end
elseif min15out == 1
    %hourrange = (108:0.5:119.5).*36000;  %For 10km Estuary
    hourrange = (48:0.25:59.5).*36000;  %For 1km Estuary
else
    hourrange = (48:59).*36000;  %For 1km Estuary
    %hourrange = (108:119).*36000;  %For 10km Estuary
end

%Plotting Flags
di = 3;
cnt=1000;
bounce=0;
buffer=0;
plotu = 1;
plotv = 0;
plotsalt = 0;

%% - Load Data and Format

tic

%Load first time step to set-up preallocated arrays
ffp = [fpath,num2str(hourrange(1)),'.mat'];
load(ffp);

%Set-Up the Arrays
usurf = zeros(length(hourrange),size(u,1));
vsurf = zeros(length(hourrange),size(u,1));
ssurf = zeros(length(hourrange),size(s,1));
distar = ones(size(v,1),1).*(((0:(ny)).*dy)');
distar = distar./1000;
distarrho = ones(size(u,1),1).*(((0:(ny-1)).*dy)');

for itloc = 1:length(hourrange)

    %Load in the data
    ts = hourrange(itloc);
    ffp = [fpath,num2str(ts),'.mat'];
    load(ffp);

    %Put v on the rho grid
    vugrid = nan(size(u,1),size(u,2));
    for i1 = 1:(size(u,1))
        for i2 = 1:size(u,2)
            vugrid(i1,i2) = (v(i1,i2)+v(i1+1,i2))/2;
        end
    end

    %Concatenate to 1 array
    usurf(itloc,:) = squeeze(u(:,76));
    vsurf(itloc,:) = squeeze(vugrid(:,76));
    ssurf(itloc,:) = squeeze(s(:,76));
    
end

%Concantenate the tidal cycle together to lengthen time particles can be
%tracked over
if numtidcyc == 1
    %Nothing to do
elseif numtidcyc == 2
    usb = cat(1,usurf,usurf);
    vsb = cat(1,vsurf,vsurf);
    ssb = cat(1,ssurf,ssurf);
    usurf = usb;  vsurf = vsb;  ssurf = ssb;
elseif numtidcyc == 3
    usb = cat(1,usurf,usurf,usurf);
    vsb = cat(1,vsurf,vsurf,vsurf);
    ssb = cat(1,ssurf,ssurf,ssurf);
    usurf = usb;  vsurf = vsb;  ssurf = ssb;
elseif numtidcyc == 4
    usb = cat(1,usurf,usurf,usurf,usurf);
    vsb = cat(1,vsurf,vsurf,vsurf,vsurf);
    ssb = cat(1,ssurf,ssurf,ssurf,ssurf);
    usurf = usb;  vsurf = vsb;  ssurf = ssb;
elseif numtidcyc == 5
    usb = cat(1,usurf,usurf,usurf,usurf,usurf);
    vsb = cat(1,vsurf,vsurf,vsurf,vsurf,vsurf);
    ssb = cat(1,ssurf,ssurf,ssurf,ssurf,ssurf);
    usurf = usb;  vsurf = vsb;  ssurf = ssb;
elseif numtidcyc == 6
    usb = cat(1,usurf,usurf,usurf,usurf,usurf,usurf);
    vsb = cat(1,vsurf,vsurf,vsurf,vsurf,vsurf,vsurf);
    ssb = cat(1,ssurf,ssurf,ssurf,ssurf,ssurf,ssurf);
    usurf = usb;  vsurf = vsb;  ssurf = ssb;
end

%Set-Up Needed Arrays for Particle Tracking
xr = 1/1000;   %Cut down distance array into 1D
yr = ones(size(v,1),1).*(((0:(ny)).*dy)')./1000;   %Cut down distance array into 1D

numhours = (size(vsurf,1));               %Number of hours of run
crosslength = max(max(distar));           %Length of cross-section in km
crosssize = size(vsurf,2);                %Number of indices in cross-section

%Set up velocity fields for tidal mean or tidal variance
if euleravg == 1
    usurfavg = mean(usurf,1);
    vsurfavg = mean(vsurf,1);
    ssurfavg = mean(ssurf,1);
end

if lmeavg == 1
    usurfavg = mean(usurf,1);
    vsurfavg = mean(vsurf,1);
    ssurfavg = mean(ssurf,1);
end

%% - Set Up Needed Values/Arrays

if halfhour == 1
    timlen = ((numhours-1).*1800)/timestep;
elseif min15out == 1
    timlen = ((numhours-1).*1800)/timestep;
else
    timlen = ((numhours-1).*3600)/timestep;
end

%Create Particle Arrays and Initial Locations
xs = zeros(timlen,numpart);
xsadd = xs;
ys = zeros(timlen,numpart);
ysadd = ys;
initialpos = linspace(0.2,0.8,100);
ys(1,:)= initialpos;
ysaddin = initialpos;
yit = numpart;

%Set Up Needed Arrays
upre = ones(1,size(xs,2));
vpre = ones(1,size(ys,2));
yupbrd = max(distar);
ydownbrd = min(distar);
disval = 0.75;

fig1 = figure;
set(gcf,'Position',[150   150   1400   400])

%% - Solve Particle Track
for it=1:(timlen-1)
    %---Time Interpolation for Velocity Field------------------------------
    if halfhour == 1
        itcur = (((it-1)*timestep)/1800)+1;
    elseif min15out == 1
        itcur = (((it-1)*timestep)/9000)+1;
    else
        itcur = (((it-1)*timestep)/3600)+1;
    end
    if itcur == round(itcur)
        if euleravg == 1
            ubar = usurfavg;
            vbar = vsurfavg;
            sbar = ssurfavg;
        else
            ubar = squeeze(usurf(itcur,:,:));
            vbar = squeeze(vsurf(itcur,:,:));
            sbar = squeeze(ssurf(itcur,:,:));
        end
        plotpartflag = 1;
    else
        if euleravg == 1
            ubar = usurfavg;
            vbar = vsurfavg;
            sbar = ssurfavg;
        else
            t1 = floor(itcur);  t2 = ceil(itcur);
    
            ubar1=squeeze(usurf(t1,:,:));   ubar2=squeeze(usurf(t2,:,:));
            vbar1=squeeze(vsurf(t1,:,:));   vbar2=squeeze(vsurf(t2,:,:));
            sbar1=squeeze(ssurf(t1,:,:));   sbar2=squeeze(ssurf(t2,:,:));
    
            wt = itcur-t1;
            ubar = (1-wt)*ubar1 + wt*ubar2;
            vbar = (1-wt)*vbar1 + wt*vbar2;
            sbar = (1-wt)*sbar1 + wt*sbar2;
            
            if lmeavg == 1
                ubar=ubar-usurfavg;
                vbar=vbar-vsurfavg;
                sbar=sbar-ssurfavg;
            end
        end
    end

    %Adds Particles at times steps after t=1
    if addpart == 1
        if it > 1 && it < 144
            if itcur == round(itcur)
                xs = cat(2,xs,xsadd);
                ys = cat(2,ys,ysadd);
                yit = (100.*(itcur-1))+1;
                ys(it,(yit:end)) = ysaddin;
            end
        end
        numpart = size(xs,2);
    end
    
    for ip=1:numpart
 
        %-- compute x-> x+dx ----------------------------------------------
        x=xs(it,ip);    y=ys(it,ip);   %current position of particle ip   
        %iy=iys(ip);    %approximate closest grid point indices
        
        %-- get u, v at (x,y) for each particle ---------------------------
        % start with u
        %Find closest index value
        [d,ix] = min(abs(distar-y));
        resval = distar(ix)-y;
        %Interpolate vbar value
        if resval == 0
            %u = ubar(y);
            u = ubar(ix);
        elseif resval > 0
            uu = ix;    ubu = ubar(uu);
            ud = ix-1;  ubd = ubar(ud);
            tdis = distar(uu)-distar(ud);
            wtv = resval/tdis;
            u = (1-wtv)*ubu + wtv*ubd;
        else
            resval = resval.*(-1);
            uu = ix+1;  ubu = ubar(uu);
            ud = ix;    ubd = ubar(ud);
            tdis = distar(uu)-distar(ud);
            wtv = resval/tdis;
            u = (1-wtv)*ubd + wtv*ubu;
        end

        % do same for v--------------------------------------------------
        %Find closest index value
        [d,ix] = min(abs(distar-y));
        resval = distar(ix)-y;
        %Interpolate vbar value
        if resval == 0
            %v = vbar(y);
            v = vbar(ix);
        elseif resval > 0
            vu = ix;    vbu = vbar(vu);
            vd = ix-1;  vbd = vbar(vd);
            tdis = distar(vu)-distar(vd);
            wtv = resval/tdis;
            v = (1-wtv)*vbu + wtv*vbd;
        else
            resval = resval.*(-1);
            vu = ix+1;  vbu = vbar(vu);
            vd = ix;    vbd = vbar(vd);
            tdis = distar(vu)-distar(vd);
            wtv = resval/tdis;
            v = (1-wtv)*vbd + wtv*vbu;
        end

        %---Integrate with Adams-Bashforth step-----------------------
            if it == 1          %Forward step for it=1
                x = x + u*timestep/1000;
                y = y + v*timestep/1000;
            else                %Adams-Bashforth step
                if ip < yit
                    x = x + (timestep/(2*1000))*((3*u) - upre(1,ip));
                    y = y + (timestep/(2*1000))*((3*v) - vpre(1,ip));
                else
                    x = x + u*timestep/1000;
                    y = y + v*timestep/1000;
                end
            end
            %---Remove Particles Pushed Beyond the Cross-Section----------
            if y > yupbrd
                x = 0;
                y = yupbrd+disval;
                u = 0; v = 0;
            elseif y < ydownbrd
                x = 0;
                y = ydownbrd-disval;
                u = 0; v = 0;
            end
            %---Save Values for Later Time Steps for Particle ip----------
            xs(it+1,ip) = x;    ys(it+1,ip) = y;
            upre(1,ip) = u;     vpre(1,ip) = v;
    end
    
    if plotpartflag == 1
        figure(fig1)
        %lenar = (-2:.5:9);
        lenar = (-25:.5:25);
        %lenar = (-10:.5:50);
        palong = ones(length(lenar),151).*lenar';
        pwidth = ones(length(lenar),151).*(distarrho'./1000);
        plotu = ones(length(lenar),151).*ubar;
        pcolor(palong,pwidth,plotu)
        shading('flat')
        hold on
        plot(xs(it,:),ys(it,:),'r.')
        %xlim([0 400])
        ylim([0 max(max(yr))])
        xlabel('Distance (km)')
        ylabel('Distance (km)')
        cmap = cmocean('balance',36);
        c1 = colorbar;
        colormap(cmap)
        caxis([-1.0 1.0])
        c1.Label.String = 'Along-Channel Velocity (m s^{-1})';
        if halfhour == 1
            title(['Particle Track - ',casename,' - Hour: ',num2str(itcur/2)]);
        elseif min15out == 1
            title(['Particle Track - ',casename,' - Hour: ',num2str(itcur/4)]);
        else
            title(['Particle Track - ',casename,' - Hour: ',num2str(itcur)]);
        end
        if isv == 1
            cnt = cnt+1;
            print([savloc,'\PartTrack\',num2str(cnt)],'-dpng');
        end
        plotpartflag = 0;
    end
    
end %End time loop---------------------------------------------------------

%% - Plot Center of Mass and 1 Standard Deviation

%Calculate Mean Position
xsquee = squeeze(xs(:,12:91));    centeravgx = mean(xsquee,2);
ysquee = squeeze(ys(:,12:91));    centeravgy = mean(ysquee,2);

%Calculate Standard Deviation
npcut = size(xsquee,2);
summationx = zeros(1,npcut);             summationy = zeros(1,npcut);
stadevalong = zeros(size(xsquee,1),1);   stadevacross = zeros(size(xsquee,1),1);
for istdvt = 1:size(xsquee,1)
    for istdvp = 1:npcut
        summationx(1,istdvp) = (xsquee(istdvt,istdvp)-centeravgx(istdvt,1))^2;
        summationy(1,istdvp) = (ysquee(istdvt,istdvp)-centeravgy(istdvt,1))^2;
    end
    sumsumx = sum(summationx);  sumsumy = sum(summationy);
    stadevalong(istdvt,1) = sqrt((sumsumx/npcut));
    stadevacross(istdvt,1) = sqrt((sumsumy/npcut));
end
onestdxdw = centeravgx(end,1)-stadevalong(end,1);   onestdxup = centeravgx(end,1)+stadevalong(end,1);
onestdx = cat(1,onestdxdw,onestdxup);
onestdydw = centeravgy(end,1)-stadevacross(end,1);  onestdyup = centeravgy(end,1)+stadevacross(end,1);
onestdy = cat(1,onestdydw,onestdyup);

saltmean = mean(ssurf,1);
pasalt = ones(length(lenar),151).*saltmean;

%% - Saves PartTrack Data
if isd == 1
    savname = [savloc,'\PartTrack\PartTrackData'];  % Name for save file
    usurfavg = mean(usurf,1);
    vsurfavg = mean(vsurf,1);
    % Saves essential variables to be saved off to avoid having to rerun
    % this multiple times
    save(savname,'numpart','xs','ys','centeravgx','centeravgy','stadevalong','stadevacross','onestdx','onestdy',...
        'usurf','vsurf','usurfavg','vsurfavg');
end

toc