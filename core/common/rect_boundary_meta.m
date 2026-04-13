function meta = rect_boundary_meta(boundary)
%RECT_BOUNDARY_META Normalize rectangular boundary aliases.

bc = char(lower(string(boundary)));
switch bc
    case {'free', 'ffff'}
        meta = struct('solver_key', 'ffff', 'title_tag', 'F', 'file_tag', 'F', 'is_levy', false);
    case {'simply', 'ssss'}
        meta = struct('solver_key', 'ssss', 'title_tag', 'S', 'file_tag', 'S', 'is_levy', true);
    case {'clamped', 'cccc'}
        meta = struct('solver_key', 'cccc', 'title_tag', 'C', 'file_tag', 'C', 'is_levy', false);
    case 'sscc'
        meta = struct('solver_key', 'sscc', 'title_tag', 'SSCC', 'file_tag', 'SSCC', 'is_levy', true);
    case 'ssff'
        meta = struct('solver_key', 'ssff', 'title_tag', 'SSFF', 'file_tag', 'SSFF', 'is_levy', true);
    case 'sssc'
        meta = struct('solver_key', 'sssc', 'title_tag', 'SSSC', 'file_tag', 'SSSC', 'is_levy', true);
    case 'sssf'
        meta = struct('solver_key', 'sssf', 'title_tag', 'SSSF', 'file_tag', 'SSSF', 'is_levy', true);
    case 'sscf'
        meta = struct('solver_key', 'sscf', 'title_tag', 'SSCF', 'file_tag', 'SSCF', 'is_levy', true);
    otherwise
        error('Unknown rectangular boundary condition: %s', boundary);
end
end
