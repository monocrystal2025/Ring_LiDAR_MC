function output = para_predraw(options)
%PARA_PREDRAW Convert raw sensitivity MAT files into ready-to-plot statistics.
%   para_predraw
%   para_predraw(Requirements=[15 .6;30 .7;60 .9],RangeWidthMrad=100)
%   para_predraw(RangeAtCrossing=true)
%   para_predraw(RebuildCache=true)
%   para_draw                         % draw the saved snapshot, without raw I/O
%
% Reads para_manifest.mat and raw/ in DataDirectory, default para_sensative_N800_range2.
% Parameter values come from the manifest cases, including independent energy
% and alpha grids; no shared multiplier grid is assumed during preprocessing.
% Reuses para_statistics_cache.mat after checking raw file size/timestamps.
% Saves para_plot_data.mat (self-contained plotting input), para_summary.mat
% (full processed curves and provenance), and para_summary.csv. No figures
% are opened. Re-run this function when raw data or analysis settings change.
%
% Row 1: equal_wD/minimum_thetaD cummax + first linear crossing;
% kappa=min(wD_ring)/min(wD_line). Three [time_s,probability] requirements.
% Row 2: blind_snr_EA static SNR-blind fraction, first nonzero meeting after
% strict annular advantage, mean relative blind reduction below it.
% Row 3: advantage_data successful-window counts, positive-pulse clusters,
% mean first-detection range at RangeWidthMrad (or the blind crossing), and
% maxima of ring/line mean pulse and cluster ratios across divergences.
% Unattainable requirements and absent crossings remain NaN.
% LinePulseDisplayFactor=0.6 reproduces advantage_draw's display adjustment;
% default 1 leaves the measured line pulse counts unchanged.
arguments
    options.DataDirectory (1,1) string = ""
    options.Requirements (3,2) double {mustBeFinite,mustBePositive} = [15 .6;30 .7;60 .9]
    options.RangeWidthMrad (1,1) double {mustBePositive,mustBeFinite} = 100
    options.RangeAtCrossing (1,1) logical = false
    options.LinePulseDisplayFactor (1,1) double {mustBePositive,mustBeFinite} = 1
    options.RebuildCache (1,1) logical = false
end
totalTimer = tic;
if options.DataDirectory == ""
    options.DataDirectory = fullfile(fileparts(mfilename('fullpath')),'para_sensative_N3000_range2');
end
assert(all(options.Requirements(:,2) <= 1),'para_predraw:Probability', ...
    'Required probabilities must be <=1.');
file = fullfile(options.DataDirectory,'para_manifest.mat');
assert(isfile(file),'para_predraw:MissingManifest','Run para_data first: %s',file);
saved = load(file,'manifest');
m = saved.manifest;
assert(m.complete,'para_predraw:Incomplete','Run para_data again to complete this sweep.');
assert(~m.truncated,'para_predraw:Truncated','Smoke-test data cannot be used for the scientific figure.');
assert(numel(m.widthsMrad) >= 2 && all(arrayfun(@(p) nnz([m.cases.parameterIndex] == p) >= 2,1:4)), ...
    'para_predraw:Grid','At least two parameter values and two divergences are required.');
if ~options.RangeAtCrossing
    assert(any(m.widthsMrad == options.RangeWidthMrad),'para_predraw:RangeWidth', ...
        'RangeWidthMrad must be an actually simulated divergence.');
end
fprintf('[1/3] Checking raw files and loading compact statistics...\n');
[rawStatistics,fileIndex,cacheInfo] = loadStatistics(m,options);
fprintf('[2/3] Calculating three advantages for %d cases...\n',numel(m.cases));
curves = cell(numel(m.cases),1);
rows = cell(numel(m.cases),1);
for j = 1:numel(m.cases)
    curves{j} = loadCurves(m,m.cases(j),options.Requirements,rawStatistics,fileIndex);
    rows{j} = summarizeCase(m,m.cases(j),curves{j},options);
end
summary = vertcat(rows{:});
summaryTable = makeTable(summary);
% The plotting snapshot deliberately contains no raw file lists, source code
% or raw photon arrays. It can be copied and plotted without the raw folder.
plotManifest = struct('N',m.N,'baseline',m.baseline,'parameters',m.parameters, ...
    'labels',m.labels,'displayScale',m.displayScale);
analysisOptions = rmfield(options,{'DataDirectory','RebuildCache'});
plotData = struct('schemaVersion',2,'manifest',plotManifest,'summary',summary, ...
    'analysisOptions',analysisOptions);
statistics = struct('manifest',m,'summary',summary,'curves',{curves}, ...
    'summaryTable',summaryTable,'options',options,'cacheInfo',cacheInfo);
fprintf('[3/3] Saving plot-ready MAT and analysis summaries...\n');
plotDataFile = fullfile(options.DataDirectory,'para_plot_data.mat');
temporaryFile = plotDataFile + ".partial.mat";
save(temporaryFile,'plotData','-v7');
movefile(temporaryFile,plotDataFile,'f');
save(fullfile(options.DataDirectory,'para_summary.mat'),'-struct','statistics','-v7.3');
writetable(summaryTable,fullfile(options.DataDirectory,'para_summary.csv'));
output = statistics;
output.plotData = plotData;
output.plotDataFile = plotDataFile;
fprintf('Preprocessing finished in %.1f s. Raw files read: %d; cached files reused: %d.\n', ...
    toc(totalTimer),cacheInfo.readCount,cacheInfo.reusedCount);
fprintf('Saved %s\nRun para_draw to plot this snapshot.\n',plotDataFile);
end

function [entries,fileIndex,info] = loadStatistics(m,options)
% Store only event times and scalar statistics, never the photon windows.
dynamicFiles = vertcat(m.cases.dynamicFiles);
staticFiles = vertcat(m.cases.staticFiles);
files = unique([dynamicFiles(:);staticFiles(:)],'stable');
fileIndex = containers.Map(cellstr(files),num2cell(1:numel(files)));
[bytes,modified] = rawFileStamps(files);
cacheFile = fullfile(options.DataDirectory,'para_statistics_cache.mat');
entries = cell(numel(files),1);
valid = false(numel(files),1);
% Increment when the meaning/algorithm of a cached statistic changes.
statisticsVersion = 1;
if ~options.RebuildCache && isfile(cacheFile)
    cached = load(cacheFile,'cache');
    if isfield(cached,'cache') && isfield(cached.cache,'version') ...
            && cached.cache.version == statisticsVersion ...
            && isequaln(cached.cache.manifestSignature,m.signature) ...
            && isequal(cached.cache.files,files)
        entries = cached.cache.entries;
        valid = cached.cache.bytes == bytes & cached.cache.modified == modified ...
            & ~cellfun(@isempty,entries);
    end
end
info = struct('readCount',nnz(~valid),'reusedCount',nnz(valid),'file',cacheFile);
fprintf('  %d unique files: %d cached, %d to read.\n', ...
    numel(files),info.reusedCount,info.readCount);
if all(valid)
    return
end
% Clear invalid entries before checkpointing, so an interruption cannot
% mark old statistics as valid with the new raw-file timestamps.
entries(~valid) = {[]};
cache = struct('version',statisticsVersion,'manifestSignature',m.signature, ...
    'files',files,'bytes',bytes,'modified',modified,'entries',{entries});
pending = find(~valid);
readTimer = tic;
lastReport = tic;
for k = 1:numel(pending)
    index = pending(k);
    [~,name] = fileparts(files(index));
    isStatic = startsWith(name,'PARA_STATIC_');
    if contains(name,'_LINE_')
        code = 'L';
    else
        code = 'R';
    end
    entries{index} = readStatistics(files(index),code,m.N,isStatic);
    if mod(k,100) == 0 || k == numel(pending)
        cache.entries = entries;
        writeStatisticsCache(cacheFile,cache);
    end
    if toc(lastReport) >= 2 || k == numel(pending)
        elapsed = toc(readTimer);
        fprintf('  Read %d/%d files (%.0f%%), %.1f s elapsed, about %.1f s remaining.\n', ...
            k,numel(pending),100*k/numel(pending),elapsed,elapsed*(numel(pending)-k)/k);
        lastReport = tic;
    end
end
end

function [bytes,modified] = rawFileStamps(files)
% List each directory once, instead of 6400 separate existence checks.
folders = strings(size(files));
for k = 1:numel(files)
    folders(k) = fileparts(files(k));
end
uniqueFolders = unique(folders);
bytes = nan(numel(files),1);
modified = bytes;
for k = 1:numel(uniqueFolders)
    listing = dir(fullfile(uniqueFolders(k),'*.mat'));
    listedFiles = fullfile(uniqueFolders(k),string({listing.name}));
    members = find(folders == uniqueFolders(k));
    [found,location] = ismember(files(members),listedFiles);
    missing = members(~found);
    assert(isempty(missing),'para_predraw:MissingRaw', ...
        'Missing raw results in %s. Run para_data to complete them.',uniqueFolders(k));
    bytes(members) = [listing(location).bytes];
    modified(members) = [listing(location).datenum];
end
end

function writeStatisticsCache(file,cache)
temporaryFile = file + ".partial.mat";
save(temporaryFile,'cache','-v7');
movefile(temporaryFile,file,'f');
end

function c = loadCurves(m,item,requirements,entries,fileIndex)
nD = numel(m.widthsMrad);
c = struct('pd',nan(nD,2,3),'blind',nan(nD,2),'meanRange',nan(nD,2), ...
    'meanPulses',nan(nD,2),'meanClusters',nan(nD,2),'successfulN',zeros(nD,2));
for d = 1:nD
    for b = 1:2
        moving = entries{fileIndex(char(item.dynamicFiles(d,b)))};
        stationary = entries{fileIndex(char(item.staticFiles(d,b)))};
        validateStatistics(moving,item,m.widthsMrad(d),m.N,false);
        validateStatistics(stationary,item,m.widthsMrad(d),m.N,true);
        c.blind(d,b) = stationary.blind;
        c.successfulN(d,b) = moving.successfulN;
        for r = 1:3
            c.pd(d,b,r) = nnz(moving.successfulTimes <= requirements(r,1))/m.N;
        end
        c.meanPulses(d,b) = moving.meanPulses;
        c.meanClusters(d,b) = moving.meanClusters;
        c.meanRange(d,b) = moving.meanRange;
    end
end
end

function r = readStatistics(file,code,n,isStatic)
names = {['detect_' code],'meta'};
if ~isStatic
    names = [names,{['first_time_' code],['effective_pulses_' code], ...
        ['first_detection_range_' code]}];
end
% Static blindness needs no photon structures or detection ranges.
raw = load(file,names{:});
assert(all(isfield(raw,names)),'para_predraw:Schema','Missing MC_EA variables: %s',file);
for k = [1,3:numel(names)]
    assert(isequal(size(raw.(names{k})),[n 1]),'para_predraw:Size', ...
        '%s must be %d-by-1: %s',names{k},n,file);
end
assert(islogical(raw.(names{1})),'para_predraw:Types','Invalid detect data: %s',file);
hit = raw.(names{1});
r = struct('meta',raw.meta,'successfulTimes',zeros(0,1), ...
    'successfulN',nnz(hit),'blind',1-mean(hit),'meanPulses',nan, ...
    'meanClusters',nan,'meanRange',nan);
if isStatic
    return
end
firstTime = raw.(names{3});
windows = raw.(names{4});
ranges = raw.(names{5});
assert(isstruct(windows) && isfield(windows,'signal_photons'), ...
    'para_predraw:Types','Invalid effective-pulse windows: %s',file);
assert(all(isfinite(firstTime(hit))) && all(isfinite(ranges(hit))), ...
    'para_predraw:Detection','Detected samples require finite time and range: %s',file);
r.successfulTimes = firstTime(hit);
if any(hit)
    successfulWindows = windows(hit);
    signal = {successfulWindows.signal_photons};
    r.meanPulses = mean(cellfun(@positivePulses,signal));
    r.meanClusters = mean(cellfun(@positiveClusters,signal));
    r.meanRange = mean(ranges(hit));
end
end

function validateStatistics(r,item,width,n,isStatic)
assert(isequal(r.meta.physics,item.physics) && r.meta.widthMrad == width ...
    && r.meta.isStatic == isStatic && r.meta.N == n, ...
    'para_predraw:Metadata','Raw metadata does not match case %s at %g mrad.',item.id,width);
end

function n = positivePulses(signal)
n = nnz(isfinite(signal(:)) & signal(:)>0);
end

function n = positiveClusters(signal)
positive = isfinite(signal(:)) & signal(:)>0;
if isempty(positive)
    n = 0;
else
    n = nnz(positive & [true;~positive(1:end-1)]);
end
end

function s = summarizeCase(m,item,c,options)
w = m.widthsMrad;
s = struct('id',item.id,'parameterIndex',item.parameterIndex, ...
    'factorIndex',item.factorIndex,'value',item.value, ...
    'minimumWidths',nan(3,2),'kappa',nan(1,3),'crossingMrad',nan, ...
    'crossingBlind',nan,'allCrossingsMrad',[],'meanBlindReduction',nan, ...
    'meanRange',nan(1,2),'rangeWidthMrad',nan, ...
    'maxPulseRatio',nan,'pulseRatioWidthMrad',nan, ...
    'maxClusterRatio',nan,'clusterRatioWidthMrad',nan);
for r = 1:3
    for b = 1:2
        s.minimumWidths(r,b) = firstMonotoneCrossing(w,c.pd(:,b,r),options.Requirements(r,2));
    end
    s.kappa(r) = s.minimumWidths(r,2)/s.minimumWidths(r,1);
end
[s.crossingMrad,s.crossingBlind,s.allCrossingsMrad] = blindCrossing(w,c.blind);
if isfinite(s.crossingMrad)
    below = w < s.crossingMrad & c.blind(:,1)>0;
    relative = nan(size(w));
    relative(below) = (c.blind(below,1)-c.blind(below,2))./c.blind(below,1);
    % Equal-weight arithmetic mean over simulated widths strictly below the
    % crossing. Zero line blind fractions are undefined and excluded; keep
    % any negative reductions rather than selecting only annular advantages.
    validRelative = relative(isfinite(relative));
    if ~isempty(validRelative)
        s.meanBlindReduction = mean(validRelative);
    end
end
if options.RangeAtCrossing
    s.rangeWidthMrad = s.crossingMrad;
    if isfinite(s.crossingMrad)
        s.meanRange = interp1(w,c.meanRange,s.crossingMrad,'linear',nan);
    end
else
    s.rangeWidthMrad = options.RangeWidthMrad;
    s.meanRange = c.meanRange(w == options.RangeWidthMrad,:);
end
pulseRatio = safeRatio(c.meanPulses(:,2),c.meanPulses(:,1)*options.LinePulseDisplayFactor);
clusterRatio = safeRatio(c.meanClusters(:,2),c.meanClusters(:,1));
[s.maxPulseRatio,s.pulseRatioWidthMrad] = finiteMaximum(pulseRatio,w);
[s.maxClusterRatio,s.clusterRatioWidthMrad] = finiteMaximum(clusterRatio,w);
end

function width = firstMonotoneCrossing(w,p,target)
% Identical inverse-envelope algorithm to minimum_thetaD.m/equal_wD.m.
p = cummax(p(:));
index = find(p >= target,1,'first');
if isempty(index)
    width = nan;
elseif index == 1
    width = w(1);
elseif abs(diff(p(index-1:index))) <= eps(max(abs(p(index-1:index))))
    width = w(index);
else
    width = interp1(p(index-1:index),w(index-1:index),target,'linear');
end
end

function [selected,blind,allCrossings] = blindCrossing(w,B)
% No smoothing/cummax is applied to the measured static blind fractions.
selected = nan;
blind = nan;
allCrossings = zeros(0,1);
delta = B(:,2)-B(:,1);
for k = 2:numel(w)
    if delta(k-1)*delta(k)<0
        cross = interp1(delta(k-1:k),w(k-1:k),0,'linear');
        allCrossings(end+1,1) = cross; %#ok<AGROW>
        if delta(k-1)<0 && ~isfinite(selected)
            selected = cross;
        end
    elseif delta(k)==0 && delta(k-1)~=0 && B(k,1)>0
        allCrossings(end+1,1) = w(k); %#ok<AGROW>
        if delta(k-1)<0 && ~isfinite(selected)
            selected = w(k);
        end
    end
end
if isfinite(selected)
    blind = mean(interp1(w,B,selected,'linear'));
end
end

function ratio = safeRatio(numerator,denominator)
ratio = nan(size(numerator));
valid = isfinite(numerator) & isfinite(denominator) & denominator>0;
ratio(valid) = numerator(valid)./denominator(valid);
end

function [value,width] = finiteMaximum(values,w)
valid = find(isfinite(values));
value = nan;
width = nan;
if ~isempty(valid)
    [value,k] = max(values(valid));
    width = w(valid(k));
end
end

function t = makeTable(s)
t = table(string({s.id}).',[s.value].','VariableNames',{'Case','ParameterValueSI'});
kappa = vertcat(s.kappa);
range = vertcat(s.meanRange);
for k = 1:3
    t.(sprintf('Kappa%d',k)) = kappa(:,k);
end
t.CrossingWidth_mrad = [s.crossingMrad].';
t.CommonBlindFraction = [s.crossingBlind].';
t.MeanRelativeBlindReduction = [s.meanBlindReduction].';
t.RangeWidth_mrad = [s.rangeWidthMrad].';
t.LineMeanRange_m = range(:,1);
t.RingMeanRange_m = range(:,2);
t.MaxPulseRatio = [s.maxPulseRatio].';
t.PulseRatioWidth_mrad = [s.pulseRatioWidthMrad].';
t.MaxClusterRatio = [s.maxClusterRatio].';
t.ClusterRatioWidth_mrad = [s.clusterRatioWidthMrad].';
end
