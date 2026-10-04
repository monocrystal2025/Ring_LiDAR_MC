function result = tmp_2(mode,options)
%TMP_2 Direction sensitivity at local-crossing and complete-ROE levels.
% The local angle is the RELATIVE angular track vs line long axis.
% The ROE angle fixes transverse heading at entry LOS; it is not the
% instantaneous relative target-beam angle later in the scan.
arguments
    mode (1,1) string = "preview"
    options.N (1,1) double {mustBeInteger,mustBeNonnegative} = 0
end
cfg = discussion.config(mode);
if options.N>0
    cfg.N = options.N;
end
cfg.outputRoot = fullfile(cfg.outputRoot,"tmp_2");
width = 120;
angles = 0:15:180;
localN = max(1024,cfg.N);
if mode=="smoke"
    localN = cfg.N;
    angles = 0:45:180;
end
localBank = discussion.bank(localN,cfg.seed);
localRaw = cell(numel(angles),2);
localHit = false(localN,numel(angles),2);
localCore = localHit;
for j = 1:numel(angles)
    for k = 1:2
        names = ["line","ring"];
        localRaw{j,k} = discussion.localCrossing(cfg,discussion.beam(names(k),width), ...
            localBank,angles(j));
        localHit(:,j,k) = isfinite(localRaw{j,k}.metrics.firstTime);
        localCore(:,j,k) = isfinite(localRaw{j,k}.metrics.coreEncounterTime);
    end
end
headings = 0:45:315;
speeds = [30,90];
bank = discussion.bank(cfg.N,cfg.seed+1);
headingRaw = cell(numel(headings),numel(speeds),2);
delta = zeros(numel(headings),numel(speeds));
low = delta; high = delta; exposure = delta;
rows = cell(numel(headings)*numel(speeds),1);
for v = 1:numel(speeds)
    for j = 1:numel(headings)
        caseCfg = cfg;
        caseCfg.heading_deg = headings(j);
        caseCfg.system.targetSpeed_m_s = speeds(v);
        for k = 1:2
            names = ["line","ring"];
            headingRaw{j,v,k} = discussion.runCase(caseCfg, ...
                discussion.beam(names(k),width),bank,"heading_"+j+"_speed_"+v);
        end
        line = headingRaw{j,v,1}.metrics;
        ring = headingRaw{j,v,2}.metrics;
        s = discussion.paired(double(isfinite(line.firstTime)), ...
            double(isfinite(ring.firstTime)),cfg);
        delta(j,v) = s.delta; low(j,v) = s.low; high(j,v) = s.high;
        exposure(j,v) = mean(line.exposure_s);
        rows{j+(v-1)*numel(headings)} = struct('EntryHeading_deg',headings(j), ...
            'Speed_m_s',speeds(v),'LineP',s.lineMean,'RingP',s.ringMean, ...
            'DeltaP',s.delta,'Low95',s.low,'High95',s.high, ...
            'MeanAvailableTime_s',exposure(j,v),'N',cfg.N);
    end
end
summary = struct2table(vertcat(rows{:}));
% Rotate the line field relative to the SAME scan and target trajectories.
% This isolates field orientation; the nominal line orientation is included.
lineRotations = 0:30:150;
rotationRaw = cell(numel(lineRotations),2);
rotationHit = false(cfg.N,numel(lineRotations),2);
for j = 1:numel(lineRotations)
    caseCfg = cfg;
    caseCfg.lineRotation_deg = lineRotations(j);
    rotationRaw{j,1} = discussion.runCase(caseCfg,discussion.beam("line",width), ...
        bank,"field_rotation_"+j);
    if j==1
        rotationRaw{j,2} = discussion.runCase(cfg,discussion.beam("ring",width),bank,"field_rotation_ring");
    else
        rotationRaw{j,2} = rotationRaw{1,2};
    end
    rotationHit(:,j,1) = isfinite(rotationRaw{j,1}.metrics.firstTime);
    rotationHit(:,j,2) = isfinite(rotationRaw{j,2}.metrics.firstTime);
end
[fig,layout,colors] = discussion.figure("Direction sensitivity",2,2);
ax = nexttile(layout);
handles = discussion.probabilityPlot(ax,angles,localCore,colors);
legend(ax,handles,["Line","Annular"],'Location','south','Box','off');
% For uniform impact b in [-D/2,D/2], a zero-width line has a projected
% crossing probability |sin(theta)|; a closed circle has probability one.
hold(ax,'on');
plot(ax,angles,100*abs(sind(angles)),':','Color',colors(4,:),'HandleVisibility','off');
discussion.styleAxes(ax,"\theta_{rel} (deg)","Core encounter (%)", ...
    "(a) Local geometry (dotted: thin line)");
ax = nexttile(layout);
discussion.probabilityPlot(ax,angles,localHit,colors);
discussion.styleAxes(ax,"\theta_{rel} (deg)","P_{cross} (%)", ...
    "(b) Local SNR feasibility at 2 km");
ax = nexttile(layout);
hold(ax,'on');
for v = 1:numel(speeds)
    errorbar(ax,headings,100*delta(:,v),100*(delta(:,v)-low(:,v)), ...
        100*(high(:,v)-delta(:,v)),'-o','Color',colors(v,:), ...
        'LineWidth',1.1,'MarkerSize',3,'DisplayName',speeds(v)+" m/s");
end
yline(ax,0,':','HandleVisibility','off');
legend(ax,'Location','best','Box','off');
discussion.styleAxes(ax,"\theta_{entry} (deg)","\Delta P_T (pp)", ...
    "(c) Heading and dwell effects; w_D=120 mrad");
ax = nexttile(layout);
plot(ax,headings,exposure,'-o','LineWidth',1.1,'MarkerSize',3);
discussion.styleAxes(ax,"\theta_{entry} (deg)","Available observation time (s)", ...
    "(d) min(15 s, ROE path length / speed)");
result = struct('localAngles_deg',angles,'localN',localN, ...
    'localRaw',{localRaw},'localHit',localHit,'localCore',localCore, ...
    'headings_deg',headings,'speeds_m_s',speeds, ...
    'headingRaw',{headingRaw},'summary',summary, ...
    'targetScanSpeedRatio',speeds/cfg.system.maxRange_m/cfg.system.scanRate_rad_s);
[fig2,layout2,colors2] = discussion.figure("Field orientation under the identical ROE scan",1,1);
ax = nexttile(layout2);
handles = discussion.probabilityPlot(ax,lineRotations,rotationHit,colors2);
legend(ax,handles,["Line","Annular"],'Box','off','Location','best');
discussion.styleAxes(ax,"Line-field rotation relative to the original frame (deg)", ...
    "P_T (%)","w_D=120 mrad for both; same scan and target tracks");
result.lineRotations_deg = lineRotations;
result.rotationRaw = rotationRaw;
result.rotationHit = rotationHit;
result.angleDefinitions = struct('entry','Transverse velocity azimuth in the entry tangent plane', ...
    'rotation','Rotation of the line long axis within each scan-point local frame', ...
    'local','Instantaneous relative angular track angle to the long axis');
fig3 = discussion2.directionDiagram();
discussion.finish(cfg,"tmp_2",result,[fig,fig2,fig3],summary);
discussion2.exportFigures(cfg,"tmp_2",[fig,fig2,fig3]);
end
