function output = lubang_data(varargin)
%LUBANG_DATA Run all four independent robustness-data generators in order.
%
%   LUBANG_DATA updates all four MAT files in the LUBANG_data folder using
%   the default preview grid and N=200 for every normalized parameter case.
%
%   Use one common Monte Carlo sample count:
%       lubang_data('N', 3000, 'Mode', 'paper')
%
%   Or set the four per-case sample counts independently:
%       lubang_data('Mode', 'paper', ...
%           'N_PRF', 3000, ...
%           'N_Energy', 2000, ...
%           'N_Atmosphere', 2500, ...
%           'N_Receiver', 1500)
%
%   Parameter order and meaning:
%       N_PRF         scan-sampling-density cases
%       N_Energy      laser-return-energy-budget cases
%       N_Atmosphere  atmospheric-quality cases
%       N_Receiver    receiver-noise-quality cases
%
%   Optional settings:
%       Mode             'smoke', 'preview' or 'paper' (default 'preview')
%       N                common fallback sample count (default 200)
%       UseParallel      enable simulator parallelism (default true)
%       Seed             common random seed (default 20260817)
%       OutputDirectory  MAT output folder (default ./LUBANG_data)
%       MakeFigure       call lubang after all data finish (default false)
%
%   The four generators remain independent and do not call one another;
%   this file is only a sequential convenience entry point.

scriptDir = fileparts(mfilename('fullpath'));
if isempty(scriptDir)
    scriptDir = pwd;
end

parser = inputParser;
parser.FunctionName = mfilename;
addParameter(parser, 'Mode', 'preview', @isTextScalar);
addParameter(parser, 'N', 200, @positiveInteger);
addParameter(parser, 'N_PRF', [], @optionalPositiveInteger);
addParameter(parser, 'N_Energy', [], @optionalPositiveInteger);
addParameter(parser, 'N_Atmosphere', [], @optionalPositiveInteger);
addParameter(parser, 'N_Receiver', [], @optionalPositiveInteger);
addParameter(parser, 'UseParallel', true, ...
    @(x) islogical(x) && isscalar(x));
addParameter(parser, 'Seed', 20260817, @nonnegativeInteger);
addParameter(parser, 'OutputDirectory', ...
    fullfile(scriptDir, 'LUBANG_data'), @isTextScalar);
addParameter(parser, 'MakeFigure', false, ...
    @(x) islogical(x) && isscalar(x));
parse(parser, varargin{:});
options = parser.Results;

mode = lower(string(options.Mode));
assert(any(mode == ["smoke", "preview", "paper"]), ...
    'lubang_data:UnknownMode', ...
    'Mode must be smoke, preview or paper.');

sampleCounts = struct( ...
    'prf', chooseCount(options.N_PRF, options.N), ...
    'energy', chooseCount(options.N_Energy, options.N), ...
    'atmosphere', chooseCount(options.N_Atmosphere, options.N), ...
    'receiver', chooseCount(options.N_Receiver, options.N));

outputDirectory = char(options.OutputDirectory);
if ~isfolder(outputDirectory)
    mkdir(outputDirectory);
end
dataFiles = struct( ...
    'prf', fullfile(outputDirectory, 'prf_omiga_data.mat'), ...
    'energy', fullfile(outputDirectory, 'energy_budget_data.mat'), ...
    'atmosphere', fullfile(outputDirectory, 'atmosphere_data.mat'), ...
    'receiver', fullfile(outputDirectory, 'receiver_noise_data.mat'));

fprintf('\n============================================================\n');
fprintf('LUBANG four-study data update\n');
fprintf('Mode: %s\n', mode);
fprintf(['N per normalized case: PRF=%d, Energy=%d, ', ...
    'Atmosphere=%d, Receiver=%d\n'], ...
    sampleCounts.prf, sampleCounts.energy, ...
    sampleCounts.atmosphere, sampleCounts.receiver);
fprintf('Output folder: %s\n', outputDirectory);
fprintf('============================================================\n\n');

commonArguments = {'Mode', char(mode), ...
    'UseParallel', options.UseParallel, 'Seed', options.Seed};
totalTimer = tic;

fprintf('Study 1/4: scan-sampling density\n');
studyTimer = tic;
prfResult = prf_omiga_data(commonArguments{:}, ...
    'N', sampleCounts.prf, 'OutputFile', dataFiles.prf);
elapsedByStudy.prf_s = toc(studyTimer);

fprintf('\nStudy 2/4: laser-return energy budget\n');
studyTimer = tic;
energyResult = energy_budget_data(commonArguments{:}, ...
    'N', sampleCounts.energy, 'OutputFile', dataFiles.energy);
elapsedByStudy.energy_s = toc(studyTimer);

fprintf('\nStudy 3/4: atmospheric quality\n');
studyTimer = tic;
atmosphereResult = atmosphere_data(commonArguments{:}, ...
    'N', sampleCounts.atmosphere, 'OutputFile', dataFiles.atmosphere);
elapsedByStudy.atmosphere_s = toc(studyTimer);

fprintf('\nStudy 4/4: receiver-noise quality\n');
studyTimer = tic;
receiverResult = receiver_noise_data(commonArguments{:}, ...
    'N', sampleCounts.receiver, 'OutputFile', dataFiles.receiver);
elapsedByStudy.receiver_s = toc(studyTimer);

output = struct();
output.mode = mode;
output.sampleCounts = sampleCounts;
output.dataFiles = dataFiles;
output.results = struct('prf', prfResult, 'energy', energyResult, ...
    'atmosphere', atmosphereResult, 'receiver', receiverResult);
output.elapsedByStudy = elapsedByStudy;
output.totalElapsed_s = toc(totalTimer);
output.figureOutput = [];

if options.MakeFigure
    orderedFiles = [string(dataFiles.prf), string(dataFiles.energy), ...
        string(dataFiles.atmosphere), string(dataFiles.receiver)];
    output.figureOutput = lubang('DataFiles', orderedFiles);
end

fprintf('\n============================================================\n');
fprintf('All four data studies finished in %.1f s.\n', ...
    output.totalElapsed_s);
fprintf('Run lubang to draw the updated 2-by-2 figure.\n');
fprintf('============================================================\n');
end

function count = chooseCount(individualCount, commonCount)
if isempty(individualCount)
    count = double(commonCount);
else
    count = double(individualCount);
end
end

function tf = isTextScalar(value)
tf = ischar(value) || (isstring(value) && isscalar(value));
end

function tf = positiveInteger(value)
tf = isnumeric(value) && isscalar(value) && isfinite(value) && ...
    value >= 1 && fix(value) == value;
end

function tf = optionalPositiveInteger(value)
tf = isempty(value) || positiveInteger(value);
end

function tf = nonnegativeInteger(value)
tf = isnumeric(value) && isscalar(value) && isfinite(value) && ...
    value >= 0 && fix(value) == value;
end
