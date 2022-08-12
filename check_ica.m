function EEG = check_ica(ss, visit)
    
    global EEG; % helps with redraw (EEG is a global variable in EEGLAB)

    workspace_prep_rs % Prepares workspace

    % Preallocation ----
    % ss = ss;
    % visit = visit;   % name of the folder for visit number

    % Creating variables ----
    this_ss = num2str(ss);

    visit_name = strcat('av', num2str(visit)); % grabs visit number
    this_ss_path = dir(fullfile(outpath, strcat('rs-', visit_name, '-', this_ss, '-ica.set')));
    this_ss_name = this_ss_path.name;
        
    % Loads in data using EEGLAB ----
    EEG = pop_loadset('filename',this_ss_name,'filepath', this_ss_path.folder);
    
    % Labels ICs for rejection ----
    EEG = pop_iclabel(EEG, 'default');
    EEG = pop_icflag(EEG, ...
        [NaN NaN;...    % brain
        0.8 1;...       % muscle (> 80% probability will reject components)
        0.8 1;...       % eye (> 80% probability will reject components)
        NaN NaN;...     % heart
        NaN NaN;...     % line noise
        NaN NaN;...     % channel noise
        NaN NaN...      % other
        ]);
    
    % Removes artifactual ICs
    this_reject = find(EEG.reject.gcompreject);
    EEG = pop_subcomp(EEG, this_reject, 0);

    eeglab redraw;
    
end