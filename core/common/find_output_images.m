function files = find_output_images(folder_path)
%FIND_OUTPUT_IMAGES Return all generated PNG files in lexical order.

info = dir(fullfile(folder_path, '*.png'));
if isempty(info)
    files = {};
    return;
end

names = sort({info.name});
files = cell(size(names));
for i = 1:numel(names)
    files{i} = fullfile(folder_path, names{i});
end
end
