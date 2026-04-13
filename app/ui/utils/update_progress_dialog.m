function update_progress_dialog(dlg, value, message)
%UPDATE_PROGRESS_DIALOG Update a progress dialog if it exists.

if isempty(dlg)
    return;
end
if ~isvalid(dlg)
    return;
end

dlg.Value = max(0, min(1, value));
dlg.Message = message;
drawnow limitrate;
end
