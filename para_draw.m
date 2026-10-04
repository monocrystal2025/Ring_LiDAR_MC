function output = para_draw(options)
%PARA_DRAW Draw the preprocessed sensitivity snapshot without raw-data access.
%   para_predraw                      % once after simulation/analysis changes
%   para_draw                         % plot and export PNG/FIG
%   para_draw(Export=false)            % display only
%   para_draw(DataDirectory="para_sensative_N1000")
%
% Reads ONLY para_plot_data.mat in DataDirectory. It does not load or scan
% raw/, para_manifest.mat, or para_statistics_cache.mat, and never invokes
% preprocessing automatically. To change Requirements, RangeWidthMrad,
% RangeAtCrossing or LinePulseDisplayFactor, run para_predraw with those
% options first. Figure labels use the settings saved in the snapshot.
% Axis labels are shared by row/column: bottom X, first-column left Y,
% last-column right Y. Each Y-axis type has common limits/ticks within a row.
% Parameter samples use equally spaced slots labeled with physical values,
% keeping the middle baseline sample centered even for nonuniform grids.
arguments
    options.DataDirectory (1,1) string = ""
    options.Visible (1,1) string {mustBeMember(options.Visible,["on","off"])} = "on"
    options.Export (1,1) logical = true
end
totalTimer = tic;
if options.DataDirectory == ""
    options.DataDirectory = fullfile(fileparts(mfilename('fullpath')),'para_sensative_N800_range2');
end
file = fullfile(options.DataDirectory,'para_plot_data.mat');
assert(isfile(file),'para_draw:MissingPlotData', ...
    'Plot-ready data not found: %s. Run para_predraw for this DataDirectory first.',file);
fprintf('[1/3] Loading plot-ready data: %s\n',file);
saved = load(file,'plotData');
assert(isfield(saved,'plotData'),'para_draw:InvalidPlotData', ...
    'Missing plotData in %s. Run para_predraw again.',file);
plotData = saved.plotData;
validatePlotData(plotData);
drawingOptions = plotData.analysisOptions;
drawingOptions.Visible = options.Visible;
fprintf('[2/3] Drawing the 3-by-4 figure...\n');
figureHandle = drawFigure(plotData.manifest,plotData.summary,drawingOptions);
output = plotData;
output.figure = figureHandle;
output.dataFile = file;
if options.Export
    fprintf('[3/3] Exporting PNG (220 dpi)...\n');
    set(figureHandle,'PaperPositionMode','auto');
    print(figureHandle,char(fullfile(options.DataDirectory,'para_sensitivity.png')),'-dpng','-r220');
    fprintf('[3/3] Saving editable FIG...\n');
    savefig(figureHandle,fullfile(options.DataDirectory,'para_sensitivity.fig'));
else
    fprintf('[3/3] Export skipped (Export=false).\n');
end
fprintf('Plotting finished in %.1f s. No raw files were read.\n',toc(totalTimer));
end

function validatePlotData(data)
assert(isstruct(data) && isscalar(data) && all(isfield(data, ...
    {'schemaVersion','manifest','summary','analysisOptions'})), ...
    'para_draw:InvalidPlotData','Invalid plotting snapshot. Run para_predraw again.');
assert(isequal(data.schemaVersion,2),'para_draw:SchemaVersion', ...
    'Plotting snapshot must contain mean blind reduction. Run para_predraw again.');
assert(all(isfield(data.manifest,{'N','baseline','parameters','labels','displayScale'})) ...
    && numel(data.manifest.baseline)==4 && numel(data.manifest.parameters)==4 ...
    && numel(data.manifest.labels)==4 && numel(data.manifest.displayScale)==4, ...
    'para_draw:Manifest','The plotting snapshot must describe four parameters.');
assert(all(isfield(data.analysisOptions,{'Requirements','RangeWidthMrad', ...
    'RangeAtCrossing','LinePulseDisplayFactor'})) ...
    && isequal(size(data.analysisOptions.Requirements),[3 2]), ...
    'para_draw:AnalysisSettings','Missing analysis settings. Run para_predraw again.');
assert(~isempty(data.summary) && all(isfield(data.summary, ...
    {'parameterIndex','value','kappa','crossingMrad','crossingBlind', ...
    'meanBlindReduction','meanRange','maxPulseRatio','maxClusterRatio'})), ...
    'para_draw:Summary','Incomplete plotting statistics. Run para_predraw again.');
for p = 1:4
    cases = data.summary([data.summary.parameterIndex] == p);
    assert(numel(cases)>=2 && isequal(size(vertcat(cases.kappa)),[numel(cases) 3]) ...
        && isequal(size(vertcat(cases.meanRange)),[numel(cases) 2]), ...
        'para_draw:SummaryShape','Invalid statistics for parameter %d. Run para_predraw again.',p);
end
end

function fig = drawFigure(m,summary,options)
fig = figure('Color','w','Position',[30 50 1060 700], ...
    'Name','Parameter sensitivity','Visible',options.Visible);
layout = tiledlayout(fig,3,4,'TileSpacing','compact','Padding','compact');
% Reserve space only at the outer edge for the last column's third ruler.
layout.Units = 'normalized';
layout.OuterPosition = [0 0 .92 1];
lineColor = [0 .447 .741];
ringColor = [.85 .325 .098];
green = [.15 .55 .25];
purple = [.55 .25 .65];
colors = lines(3);
mainAxes = gobjects(3,4);
outerRightAxes = gobjects(3,4);
for p = 1:4
    s = summary([summary.parameterIndex] == p);
    [physicalValues,order] = sort([s.value]);
    s = s(order);
    tickLabels = compose('%g',physicalValues*m.displayScale(p));
    x = 1:numel(s);
    baseline = interp1(physicalValues,x,m.baseline(p),'linear');
    % Slots are uniformly spaced; labels retain the actual parameter values.
    % Symmetric limits place the manuscript baseline at the axis center.
    radius = max(abs(x-baseline))*1.22;
    limits = baseline+[-radius radius];
    for row = 1:3
        ax = nexttile(layout,(row-1)*4+p);
        set(ax,'FontName','Helvetica','Units','normalized', ...
            'FontSize',13,'LabelFontSizeMultiplier',1,'TitleFontSizeMultiplier',1, ...
            'LineWidth',1.5,'Box','on','XColor','k','YColor','k', ...
            'Tag',sprintf('para_main_r%d_c%d',row,p));
        mainAxes(row,p) = ax;
        hold(ax,'on');
        if row == 1
            kappa = vertcat(s.kappa);
            handles = gobjects(1,3);
            labels = strings(1,3);
            for r = 1:3
                handles(r) = plot(ax,x,kappa(:,r),'-o','Color',colors(r,:), ...
                    'LineWidth',2.5,'MarkerSize',4);
                labels(r) = sprintf('t<=%g s, P_D>=%g%%',options.Requirements(r,1),100*options.Requirements(r,2));
            end
            yline(ax,1,'k--','LineWidth',2.5,'HandleVisibility','off');
            ylabel(ax,'\kappa = w_{D,R}^{min}/w_{D,L}^{min}');
            legend(ax,handles,labels,'Location','best','FontSize',13,'Box','off');
        else
            if row == 2
                % Both beams share the same divergence at the crossing.
                values = [s.crossingMrad].';
                barLabels = {'Crossing w_D'};
                right1 = 100*[s.crossingBlind];
                right2 = 100*[s.meanBlindReduction];
                leftLabel = 'Equal-blind w_D (mrad)';
                rightLabel1 = 'Common blind (%)';
                rightLabel2 = 'Mean blind reduction (%)';
            else
                values = vertcat(s.meanRange);
                barLabels = {'Line','Annular'};
                right1 = [s.maxPulseRatio];
                right2 = [s.maxClusterRatio];
                leftLabel = 'Mean first-detection range (m)';
                rightLabel1 = 'Max pulse ratio (R/L)';
                rightLabel2 = 'Max cluster ratio (R/L)';
            end
            bars = bar(ax,x,values,.8,'grouped','EdgeColor','none');
            bars(1).FaceColor = lineColor;
            if row == 3
                bars(2).FaceColor = ringColor;
            end
            ylabel(ax,leftLabel);
            finiteBars = values(isfinite(values));
            if ~isempty(finiteBars) && max(finiteBars)>0
                ylim(ax,[0,1.15*max(finiteBars)]);
            end
            yyaxis(ax,'right');
            curve1 = plot(ax,x,right1,'-s','Color',green,'LineWidth',2.5,'MarkerSize',4);
            ylabel(ax,rightLabel1,'Color',green);
            ax.YColor = green;
            yyaxis(ax,'left');
            ax.YColor = 'k';
            % A wider transparent axes places the THIRD ruler outside the
            % first right ruler. Its expanded XLim preserves data alignment.
            extra = axes(fig,'Units','normalized','Position',ax.Position, ...
                'Color','none','YAxisLocation','right','XColor','none', ...
                'YColor',purple,'FontSize',13,'FontName','Helvetica', ...
                'LabelFontSizeMultiplier',1,'TitleFontSizeMultiplier',1, ...
                'LineWidth',1.5,'Box','off','HitTest','off', ...
                'Tag',sprintf('para_outer_right_r%d_c%d',row,p));
            outerRightAxes(row,p) = extra;
            hold(extra,'on');
            curve2 = plot(extra,x,right2,'-^','Color',purple,'LineWidth',2.5,'MarkerSize',4);
            extra.XLim = limits;
            extra.XTick = [];
            ylabel(extra,rightLabel2,'Color',purple);
            % Legend proxies keep one legend even though the third axis is separate.
            proxy = plot(ax,nan,nan,'-^','Color',purple,'LineWidth',2.5);
            legend(ax,[bars(:);curve1;proxy],[barLabels,{rightLabel1,rightLabel2}], ...
                'Location','northoutside','NumColumns',2,'FontSize',13,'Box','off');
            curve2.HandleVisibility = 'off';
            if row == 3
                if options.RangeAtCrossing
                    rangeText = 'Range evaluated at blind crossing';
                else
                    rangeText = sprintf('Range at w_D=%g mrad',options.RangeWidthMrad);
                end
                text(ax,.03,.96,rangeText,'Units','normalized', ...
                    'VerticalAlignment','top','FontSize',13,'Color',[.3 .3 .3]);
            end
        end
        xlim(ax,limits);
        xticks(ax,x);
        if row == 3
            xticklabels(ax,tickLabels);
            xlabel(ax,m.labels(p),'Interpreter','tex');
        else
            xticklabels(ax,[]);
        end
        grid(ax,'on');
        xline(ax,baseline,':','Color',[.4 .4 .4],'LineWidth',2.5,'HandleVisibility','off');
    end
end
% Include labels, legends, annotations and both right rulers explicitly.
set(findall(fig,'-property','FontName'),'FontName','Helvetica');
set(findall(fig,'-property','FontSize'),'FontSize',13);
formatSharedYAxes(mainAxes,outerRightAxes);
positionOverlayAxes(mainAxes,outerRightAxes);
fig.SizeChangedFcn = @(~,~) positionOverlayAxes(mainAxes,outerRightAxes);
drawnow;
end

function formatSharedYAxes(mainAxes,outerRightAxes)
% Synchronize like quantities only: left, inner right, and outer right each
% keep their own units and scale, shared across all four columns of a row.
drawnow;
for row = 1:3
    synchronizeRulers(mainAxes(row,:),1);
    if row > 1
        synchronizeRulers(mainAxes(row,:),2);
        synchronizeRulers(outerRightAxes(row,:),1);
    end
    for column = 1:4
        if column ~= 1
            mainAxes(row,column).YAxis(1).Label.String = '';
            mainAxes(row,column).YAxis(1).TickLabels = {};
        end
        if row > 1 && column ~= 4
            % Keep only the main axes' black right edge. The third-axis
            % curve stays visible while its ruler and ticks are hidden.
            mainAxes(row,column).YAxis(2).Color = 'k';
            mainAxes(row,column).YAxis(2).Label.String = '';
            mainAxes(row,column).YAxis(2).TickLabels = {};
            outerRightAxes(row,column).YAxis(1).Label.String = '';
            outerRightAxes(row,column).YAxis(1).TickLabels = {};
            outerRightAxes(row,column).YAxis(1).Visible = 'off';
        end
    end
end
end

function synchronizeRulers(axesHandles,rulerIndex)
limits = zeros(numel(axesHandles),2);
for k = 1:numel(axesHandles)
    limits(k,:) = axesHandles(k).YAxis(rulerIndex).Limits;
end
sharedLimits = [min(limits(:,1)),max(limits(:,2))];
for k = 1:numel(axesHandles)
    axesHandles(k).YAxis(rulerIndex).Limits = sharedLimits;
    axesHandles(k).YAxis(rulerIndex).LimitsMode = 'manual';
end
% Let MATLAB choose readable ticks on the shared range, then freeze exactly
% the same tick vector on every corresponding ruler, including hidden labels.
reference = axesHandles(1).YAxis(rulerIndex);
reference.TickValuesMode = 'auto';
drawnow;
sharedTicks = reference.TickValues;
for k = 1:numel(axesHandles)
    axesHandles(k).YAxis(rulerIndex).TickValues = sharedTicks;
    axesHandles(k).YAxis(rulerIndex).TickValuesMode = 'manual';
end
end

function positionOverlayAxes(mainAxes,outerRightAxes)
% Follow the compact tiles, including when the figure is resized.
if ~all(isgraphics(mainAxes),'all')
    return
end
drawnow;
for row = 2:3
    for column = 1:4
        ax = mainAxes(row,column);
        extra = outerRightAxes(row,column);
        % Axes positions are normalized to the tiled layout. Apply its
        % reserved-width transform explicitly for figure-level overlays;
        % getpixelposition does not apply this transform consistently at print.
        tilePosition = ax.Parent.OuterPosition;
        position = [tilePosition(1:2)+tilePosition(3:4).*ax.Position(1:2), ...
            tilePosition(3:4).*ax.Position(3:4)];
        extra.Units = 'normalized';
        limits = ax.XLim;
        if column == 4
            % Expand width and X range equally to preserve data alignment.
            offset = .043;
            extra.Position = position+[0 0 offset 0];
            extra.XLim = [limits(1),limits(1)+diff(limits)*(position(3)+offset)/position(3)];
        else
            extra.Position = position;
            extra.XLim = limits;
        end
    end
end
end
