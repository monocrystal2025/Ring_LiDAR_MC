%% Generate and save the LiDAR beam-vector path
f = 5e3;
R = 1000;
omega = 2 * pi;
fasan_D_list = 5e-3:5e-3:200e-3;

load_dir = "G:\BeamVEC_NEW";

if omega > 1.5 * pi
    omi_type = "fast";
else
    omi_type = "slow";
end

if ~isfolder(load_dir)
    mkdir(load_dir);
end

% Suppress and close the visualization created by lidarpath during batch
% generation so that figures do not accumulate in the loop.
originalFigureVisibility = get(groot, "DefaultFigureVisible");
figureCleanup = onCleanup(@() set(groot, ...
    "DefaultFigureVisible", originalFigureVisibility));
set(groot, "DefaultFigureVisible", "off");

numberOfPaths = numel(fasan_D_list);
for pathIndex = 1:numberOfPaths
    fasan_D = fasan_D_list(pathIndex);
    D = 2 * R * tan(fasan_D / 2);

    figuresBefore = findall(groot, "Type", "figure");
    beam_vectors = lidarpath(f, D, omega, R);
    % Store paths in the legacy top-to-bottom order used by the MC runs.
    beam_vectors = flipud(beam_vectors);
    figuresAfter = findall(groot, "Type", "figure");
    close(setdiff(figuresAfter, figuresBefore));

    file_name = "beam_vector_" + omi_type + "_" + num2str(f) + ...
        "Hz_h3000m_D" + num2str(fasan_D * 1e6) + "u.npy";
    file_path = fullfile(load_dir, file_name);

    writeNpyDouble(file_path, beam_vectors);
    fprintf("[%d/%d] Saved %d-by-3 array: %s\n", ...
        pathIndex, numberOfPaths, size(beam_vectors, 1), file_path);
end

clear figureCleanup
fprintf("Finished generating %d NPY files in %s\n", ...
    numberOfPaths, load_dir);

function writeNpyDouble(filePath, data)
%WRITENPYDOUBLE Write a real double matrix as a NumPy .npy v1.0 file.
arguments
    filePath (1, 1) string
    data (:, :) double {mustBeReal, mustBeFinite}
end

shapeText = sprintf("%d, %d", size(data, 1), size(data, 2));
header = sprintf( ...
    '{''descr'': ''<f8'', ''fortran_order'': False, ''shape'': (%s), }', ...
    shapeText);

% NPY v1.0 requires the magic bytes, version, two-byte header length, and
% header (including its newline) to end on a 16-byte boundary.
preambleLength = 10;
paddingLength = mod(16 - mod(preambleLength + numel(header) + 1, 16), 16);
header = [header, repmat(' ', 1, paddingLength), newline];
headerBytes = uint8(header);

if numel(headerBytes) > intmax("uint16")
    error("path_init:NpyHeaderTooLong", ...
        "The generated NPY header exceeds the version 1.0 size limit.");
end

fileId = fopen(filePath, "w", "ieee-le");
if fileId < 0
    error("path_init:FileOpenFailed", ...
        "Unable to open the output file: %s", filePath);
end
cleanupObject = onCleanup(@() fclose(fileId));

fwrite(fileId, [uint8(147), uint8('NUMPY')], "uint8");
fwrite(fileId, [1, 0], "uint8");
fwrite(fileId, numel(headerBytes), "uint16");
fwrite(fileId, headerBytes, "uint8");

% Transpose before fwrite so MATLAB column-major storage is emitted in
% NumPy C-order while preserving the matrix's N-by-3 shape.
fwrite(fileId, data.', "double");
end
