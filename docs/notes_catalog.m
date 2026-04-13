function notes = notes_catalog(domain_type, boundary)
%NOTES_CATALOG Provide concise GUI notes for the current setup.

domain_type = char(lower(string(domain_type)));
boundary = char(lower(string(boundary)));

switch domain_type
    case 'rect'
        domain_line = 'rect: FFFF and CCCC stay on the legacy solver path; SSSS uses Navier; SSCC / SSFF / SSSC / SSSF / SSCF use analytic Levy separation.';
    case 'circ'
        domain_line = 'circ: circular-plate backend uses Bessel and modified-Bessel characteristic equations.';
    otherwise
        domain_line = 'unknown domain.';
end

switch boundary
    case {'free', 'ffff'}
        bc_line = 'free / FFFF: keeps the existing sparse free-edge branch.';
    case {'clamped', 'cccc'}
        bc_line = 'clamped / CCCC: keeps the existing clamped finite-difference branch.';
    case {'simply', 'ssss'}
        bc_line = 'SSSS: keeps the explicit Navier mode formula.';
    case 'sscc'
        bc_line = 'SSCC: Levy analytic branch with x-pair simply supported and y-pair clamped.';
    case 'ssff'
        bc_line = 'SSFF: Levy analytic branch with x-pair simply supported and y-pair free.';
    case 'sssc'
        bc_line = 'SSSC: Levy analytic branch with y=0 simply supported and y=b clamped.';
    case 'sssf'
        bc_line = 'SSSF: Levy analytic branch with y=0 simply supported and y=b free.';
    case 'sscf'
        bc_line = 'SSCF: Levy analytic branch with y=0 clamped and y=b free.';
    otherwise
        bc_line = 'boundary preset not recognized.';
end

notes = { ...
    'This GUI preserves the original image composition by previewing generated PNG files directly.', ...
    domain_line, ...
    bc_line, ...
    'The GUI layer only manages parameters, run flow, preview selection, and PNG export.', ...
    'Visualization helpers are separated from the solver folders, so color bar and display formatting are no longer stored in core.'};
end
