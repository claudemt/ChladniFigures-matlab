function result = run_chladni_generation(project_root, params)
%RUN_CHLADNI_GENERATION Unified backend for GUI-triggered figure generation.

if nargin < 1 || isempty(project_root)
    project_root = pwd;
end

output_root = fullfile(project_root, 'output');
if ~exist(output_root, 'dir')
    mkdir(output_root);
end

folder_name = sprintf('%s_%s_nu_%s', ...
    lower(string(params.type)), lower(string(params.boundary)), local_num_tag(params.nu));
output_folder = fullfile(output_root, char(folder_name));
prepare_output_folder(output_folder);

switch char(lower(string(params.type)))
    case {'rect', 'square'}
        chladni_formula_rect(params.nu, params.k, params.n, params.normalize, output_folder, params.boundary);
    case {'circ', 'circle'}
        chladni_formula_circ(params.nu, params.k, params.n, params.normalize, output_folder, params.boundary);
    otherwise
        error('Unknown domain type: %s', params.type);
end

result = struct();
result.output_folder = output_folder;
result.files = find_output_images(output_folder);
end

function tag = local_num_tag(x)
tag = sprintf('%.6g', x);
end

function prepare_output_folder(output_folder)
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
    return;
end

png_info = dir(fullfile(output_folder, '*.png'));
for i = 1:numel(png_info)
    delete(fullfile(output_folder, png_info(i).name));
end
end
