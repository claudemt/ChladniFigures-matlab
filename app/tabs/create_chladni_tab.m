function tab = create_chladni_tab(tab_group, project_root)
%CREATE_CHLADNI_TAB Build the Chladni GUI tab.

app_figure = ancestor(tab_group, 'figure');

tab = uitab(tab_group, 'Title', 'chladni figures');
root = uigridlayout(tab, [1 2]);
root.ColumnWidth = {292, '1x'};
root.RowHeight = {'1x'};
root.Padding = [8 8 8 8];
root.ColumnSpacing = 10;

left_panel = uipanel(root, 'Title', 'controls');
left_panel.Layout.Row = 1;
left_panel.Layout.Column = 1;
left_grid = uigridlayout(left_panel, [5 1]);
left_grid.RowHeight = {'fit', 'fit', 'fit', 140, '1x'};
left_grid.ColumnWidth = {'1x'};
left_grid.Padding = [8 8 8 8];
left_grid.RowSpacing = 8;

physical_panel = uipanel(left_grid, 'Title', 'physical parameters');
physical_panel.Layout.Row = 1;
physical_panel.Layout.Column = 1;
physical_grid = uigridlayout(physical_panel, [4 1]);
physical_grid.RowHeight = {'fit', 'fit', 'fit', 'fit'};
physical_grid.Padding = [8 8 8 8];
physical_grid.RowSpacing = 5;

type_dd = create_dropdown_control(physical_grid, 'domain', {'rect', 'circ'}, 'rect', 'Square plate or circular plate.');
boundary_dd = create_dropdown_control(physical_grid, 'boundary', rect_boundary_options(), 'ffff', 'Boundary condition preset for the selected plate.');
nu_field = create_numeric_control(physical_grid, 'nu', 0.225, 'Poisson ratio in (0, 0.5).');
mode_count = create_numeric_control(physical_grid, 'number of modes', 10, 'How many output mode figures to generate.');

numerical_panel = uipanel(left_grid, 'Title', 'numerical / display parameters');
numerical_panel.Layout.Row = 2;
numerical_panel.Layout.Column = 1;
numerical_grid = uigridlayout(numerical_panel, [3 1]);
numerical_grid.RowHeight = {'fit','fit','fit'};
numerical_grid.Padding = [8 8 8 8];
numerical_grid.RowSpacing = 5;

grid_n = create_numeric_control(numerical_grid, 'grid size', 200, 'Grid resolution passed to the plotting backend.');
normalize_dd = create_dropdown_control(numerical_grid, 'normalize display', {'on', 'off'}, 'on', 'Use the backend display normalization switch.');
auto_preview_dd = create_dropdown_control(numerical_grid, 'preview on run', {'first mode', 'last mode'}, 'first mode', 'Choose which generated image to preview immediately after a run.');

action_panel = uipanel(left_grid, 'Title', 'actions');
action_panel.Layout.Row = 3;
action_panel.Layout.Column = 1;
action_grid = uigridlayout(action_panel, [1 1]);
action_grid.RowHeight = {28};
action_grid.Padding = [8 8 8 8];
button_block = create_button_row(action_grid, @run_simulation, @reset_defaults, @export_selected);
button_block.run.Tooltip = 'Generate Chladni figures with the current settings.';
button_block.export.Tooltip = 'Export the currently previewed PNG to a location you choose.';

status_box = uitextarea(left_grid, 'Editable', 'off', 'Value', {'status: ready'}, 'FontName', 'Courier New');
status_box.Layout.Row = 4;
status_box.Layout.Column = 1;

right_grid = uigridlayout(root, [2 1]);
right_grid.Layout.Row = 1;
right_grid.Layout.Column = 2;
right_grid.RowHeight = {'1x', 122};
right_grid.ColumnWidth = {'1x'};
right_grid.Padding = [0 0 0 0];
right_grid.RowSpacing = 8;

preview_panel = uipanel(right_grid, 'Title', 'preview');
preview_panel.Layout.Row = 1;
preview_panel.Layout.Column = 1;
preview_grid = uigridlayout(preview_panel, [1 2]);
preview_grid.RowHeight = {'1x'};
preview_grid.ColumnWidth = {186, '1x'};
preview_grid.Padding = [6 6 6 6];
preview_grid.RowSpacing = 4;
preview_grid.ColumnSpacing = 6;

mode_list = uilistbox(preview_grid, 'Items', {'(no results yet)'}, 'ValueChangedFcn', @on_select_mode);
mode_list.Layout.Row = 1;
mode_list.Layout.Column = 1;

preview_axes = uiaxes(preview_grid);
preview_axes.Layout.Row = 1;
preview_axes.Layout.Column = 2;
apply_axes_style(preview_axes);
preview_axes.XTick = [];
preview_axes.YTick = [];
title(preview_axes, '$\mathrm{preview}$', 'Interpreter', 'latex');
text(preview_axes, 0.5, 0.5, '$\mathrm{run\ the\ model\ to\ generate\ figures}$', ...
    'Interpreter', 'latex', 'HorizontalAlignment', 'center');
preview_axes.XLim = [0 1];
preview_axes.YLim = [0 1];

notes_box = uitextarea(right_grid, 'Editable', 'off');
notes_box.Layout.Row = 2;
notes_box.Layout.Column = 1;
notes_box.Value = notes_catalog('rect', 'ffff');

state = struct();
state.generated_files = {};
state.current_file = '';
state.current_folder = '';

install_boundary_items();
type_dd.ValueChangedFcn = @(~,~) on_domain_changed();
boundary_dd.ValueChangedFcn = @(~,~) on_boundary_changed();

    function run_simulation(~, ~)
        dlg = create_progress_dialog(app_figure, 'running chladni generator');
        try
            params = read_params();
            update_progress_dialog(dlg, 0.15, 'validating parameters');
            notes_box.Value = notes_catalog(params.type, params.boundary);

            update_progress_dialog(dlg, 0.35, 'generating figures');
            result = run_chladni_generation(project_root, params);
            state.generated_files = result.files;
            state.current_folder = result.output_folder;

            update_progress_dialog(dlg, 0.78, 'updating preview');
            mode_items = make_mode_items(result.files);
            if isempty(mode_items)
                error('No PNG output files were generated.');
            end
            mode_list.Items = mode_items;
            if strcmp(auto_preview_dd.Value, 'last mode')
                mode_list.Value = mode_items{end};
            else
                mode_list.Value = mode_items{1};
            end
            refresh_preview_from_list();

            status_box.Value = { ...
                sprintf('status            : complete'), ...
                sprintf('domain            : %s', params.type), ...
                sprintf('boundary          : %s', params.boundary), ...
                sprintf('nu                : %.6f', params.nu), ...
                sprintf('number of modes   : %d', params.k), ...
                sprintf('grid size         : %d', params.n), ...
                sprintf('normalize display : %s', normalize_dd.Value), ...
                sprintf('generated files   : %d', numel(result.files)), ...
                sprintf('output folder     : %s', result.output_folder)};

            update_progress_dialog(dlg, 1.0, 'run complete');
            pause(0.05);
            close_progress_dialog(dlg);
        catch ME
            close_progress_dialog(dlg);
            status_box.Value = {'status: failed', ['message: ' ME.message]};
            uialert(app_figure, sprintf('Chladni run failed:\n%s', ME.message), 'Run failed', 'Icon', 'error');
        end
    end

    function params = read_params()
        params = struct();
        params.type = char(lower(string(type_dd.Value)));
        params.boundary = char(lower(string(boundary_dd.Value)));
        params.nu = nu_field.Value;
        params.k = round(mode_count.Value);
        params.n = round(grid_n.Value);
        params.normalize = strcmp(normalize_dd.Value, 'on');

        if ~isfinite(params.nu) || params.nu <= 0 || params.nu >= 0.5
            error('nu must be in (0, 0.5).');
        end
        if ~isfinite(params.k) || params.k < 1
            error('number of modes must be a positive integer.');
        end
        if ~isfinite(params.n) || params.n < 32
            error('grid size must be at least 32.');
        end
    end

    function reset_defaults(~, ~)
        type_dd.Value = 'rect';
        install_boundary_items();
        boundary_dd.Value = 'ffff';
        nu_field.Value = 0.225;
        mode_count.Value = 10;
        grid_n.Value = 200;
        normalize_dd.Value = 'on';
        auto_preview_dd.Value = 'first mode';
        notes_box.Value = notes_catalog('rect', 'ffff');
        status_box.Value = {'status: ready'};
        mode_list.Items = {'(no results yet)'};
        cla(preview_axes);
        apply_axes_style(preview_axes);
        preview_axes.Visible = 'on';
        preview_axes.XTick = [];
        preview_axes.YTick = [];
        title(preview_axes, '$\mathrm{preview}$', 'Interpreter', 'latex');
        text(preview_axes, 0.5, 0.5, '$\mathrm{run\ the\ model\ to\ generate\ figures}$', ...
            'Interpreter', 'latex', 'HorizontalAlignment', 'center');
        preview_axes.XLim = [0 1];
        preview_axes.YLim = [0 1];
        state.generated_files = {};
        state.current_file = '';
        state.current_folder = '';
    end

    function export_selected(~, ~)
        if isempty(state.current_file) || ~isfile(state.current_file)
            uialert(app_figure, 'No preview image is selected yet.', 'Nothing to export', 'Icon', 'warning');
            return;
        end
        [~, default_name, default_ext] = fileparts(state.current_file);
        [file_name, file_path] = uiputfile({'*.png', 'PNG image (*.png)'}, 'Export current preview', [default_name default_ext]);
        if isequal(file_name, 0) || isequal(file_path, 0)
            return;
        end
        copyfile(state.current_file, fullfile(file_path, file_name));
        uialert(app_figure, 'Preview image exported successfully.', 'Export complete', 'Icon', 'info');
    end

    function on_select_mode(~, ~)
        refresh_preview_from_list();
    end

    function refresh_preview_from_list()
        if isempty(state.generated_files)
            return;
        end
        idx = find(strcmp(mode_list.Items, mode_list.Value), 1, 'first');
        if isempty(idx), idx = 1; end
        state.current_file = state.generated_files{idx};
        render_png_preview(preview_axes, state.current_file);
    end

    function on_domain_changed()
        install_boundary_items();
        notes_box.Value = notes_catalog(char(type_dd.Value), char(boundary_dd.Value));
    end

    function on_boundary_changed()
        notes_box.Value = notes_catalog(char(type_dd.Value), char(boundary_dd.Value));
    end

    function install_boundary_items()
        if strcmp(char(type_dd.Value), 'rect')
            items = rect_boundary_options();
            boundary_dd.Items = items;
            if ~any(strcmp(items, boundary_dd.Value))
                boundary_dd.Value = 'ffff';
            end
            boundary_dd.Tooltip = ['rect supports legacy FFFF/CCCC plus analytic SSSS, SSCC, SSFF, SSSC, SSSF, SSCF.'];
        else
            items = {'free', 'simply', 'clamped'};
            boundary_dd.Items = items;
            if ~any(strcmp(items, boundary_dd.Value))
                boundary_dd.Value = 'free';
            end
            boundary_dd.Tooltip = 'circ uses free / simply / clamped edge conditions.';
        end
    end
end

function items = make_mode_items(files)
items = cell(size(files));
for i = 1:numel(files)
    [~, name, ext] = fileparts(files{i});
    items{i} = [name ext];
end
end
