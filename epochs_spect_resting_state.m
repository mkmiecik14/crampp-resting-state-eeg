% EEG ICA cleaning and Spectral Decomposition Step 4 - Resting State
% Matt Kmiecik
% Started 20 FEB 2022

workspace_prep % Prepares workspace

% Preallocation ----
num_iters = size(NUM, 1);       % number of participants in this batch
i=1;                            % for testing purposes
visit = 'assessment-visit-1';   % name of the folder for visit number
csd_switch = 1;                 % 1 == CSD will be computed
plot_switch = 1;                % 1 == PSD plots will be saved

for i = 1:num_iters
    
    % Creating variables ----
    visit_name = strcat('av', visit(end)); % grabs visit number
    this_ss = NUM(i);
    this_ss_path = dir(fullfile(outpath, strcat('rs-', visit_name, '-', num2str(this_ss), '-ica.set')));
    this_ss_name = this_ss_path.name;
        
    % Loads in data using EEGLAB ----
    EEG = pop_loadset('filename',this_ss_name,'filepath', this_ss_path.folder);
    
    % Labels ICs for rejection ----
    EEG = pop_iclabel(EEG, 'default');
    EEG = pop_icflag(EEG, ...
        [NaN NaN;...    % brain
        NaN NaN;...     % muscle
        0.8 1;...       % eye ( > 80% probability will reject components)
        NaN NaN;...     % heart
        NaN NaN;...     % line noise
        NaN NaN;...     % channel noise
        NaN NaN...      % other
        ]);
    
    % Removes artifactual ICs
    this_reject = find(EEG.reject.gcompreject);
    EEG = pop_subcomp(EEG, this_reject, 0);
    
    %eeglab redraw    
   
    if csd_switch == 1
    
        % Compute surface laplacian spatial filter ----
        % Calculating Surface Laplacian via CSD toolbox functions
        % http://psychophysiology.cpmc.columbia.edu/software/CSDtoolbox/tutorial.html#PrepareInput
        chan_mont = cell(EEG.nbchan,1); % initializes cell array
        % Fills cell array with electrode labels
        for j = 1:length(chan_mont)
            chan_mont(j) = cellstr(EEG.chanlocs(j).labels);
        end
        % Derives spherical coordinates via CSD toolbox fx ExtractMontage()
        csd_mont = ExtractMontage('10-5-System_Mastoids_EGI129.csd', chan_mont);
        % To view: MapMontage(csd_mont)
        [G, H] = GetGH(csd_mont); % Calculates G and H matrices

        % Applies surface laplacian to EEG data
        % lambda left at default, 10 = cm head size, so units are microvolt/cm^2
        EEG.data = CSD(EEG.data, G, H, 1.0e-5, 10); 
    
    else
        disp('CSD skipped....');
    end
    
    % Spectral decomposition ----
    % Epoching ----
    % selecting the stimulation blocks (60 second epochs)
    blocks = {'S111' 'S102' 'S103' 'S114' 'S105' 'S116' 'S117' 'S108'};
    blocks_end = {'S211' 'S202' 'S203' 'S214' 'S205' 'S216' 'S217' 'S208'};
    
    % preallocates arrays
    this_spectra = zeros(EEG.nbchan, EEG.srate+1, length(blocks));
    this_freqs = zeros(EEG.srate+1, 1, length(blocks));
    
    % TO DO: In order for this script to accomodate rejected epochs, the
    % epoch start stop must be determined by a 
    % starting trigger (e.g., S111) and its ending e.g., (S211)
    
    % Getting start and stop latencies of the stimulation blocks
    start_times = zeros(1, length(blocks)); % initializes vector
    end_times = start_times; % initializes vector
    for k = 1:length(blocks)
        % Retrieving the latency (in seconds) of each block start
        start_times(1,k) = (eeg_getepochevent(EEG, blocks(k), [], 'latency'))/1000;
        % Retrieving the latency (in seconds) of each block end
        end_times(1,k) = (eeg_getepochevent(EEG, blocks_end(k), [], 'latency'))/1000;
    end
    
    % Calculates the duration of each block (in order)
    block_durations = round(end_times - start_times);
    
    for j = 1:length(blocks)
        
        this_epoch = [0 block_durations(j)]; % Sets the epoch duration
        
        try % This will run if the block exists
            % Selects blocks (in order)
            this_EEG = pop_epoch(EEG,blocks(j),this_epoch,'epochinfo', 'yes');
            % Spectral decomposition here
            [this_spectra(:,:,j), this_freqs(:,:,j)] = spectopo(...
                this_EEG.data(:,:), ... 
                0, ... % frames per epoch (default 0 = data length)
                this_EEG.srate, ... % sampling rate
                'winsize', 2*this_EEG.srate, ... % winsize is 2 seconds
                'overlap', this_EEG.srate, ... % overlap is 1 second
                'plot','off'... % toggles plot
                );
            
            % PROBABLY HERE IS A GOOD PLACE FOR Center of Gravity CALC
            
            % Plots for troubleshooting (if needed)
            if plot_switch == 1
                figure; pop_spectopo(this_EEG,1,[],'EEG','freq',[4 8 12 25 30],'freqrange',[0 75],'electrodes','on');
                saveas(gcf,...
                    fullfile(outpath,...
                    strcat('rs-', visit_name, '-', num2str(this_ss),...
                    '-',blocks{j},'.png')));
                close; % closes figure
            else
                % plots not saved
            end
        
        catch % if the block is missing, then the matrix is filled with NaN
            this_spectra(:,:,j) = NaN; % fills with missing values
            this_freqs(:,:,j) = NaN; % fills with missing values
     
        end
            
    end
    
    % Saving out results ----
    % combines into one variable
    % stimulation results are stored with each index being the
    % blocks in the following order:
    % 'S111' - eyes OPEN
    % 'S102' - eyes CLOSED
    % 'S103' - eyes CLOSED
    % 'S114' - eyes OPEN
    % 'S105' - eyes CLOSED
    % 'S116' - eyes OPEN
    % 'S117' - eyes OPEN
    % 'S108' - eyes CLOSED
    spec_res.spectra   = this_spectra; % 3D mat of spectra
    spec_res.freqs     = this_freqs;   % 3D mat of freqs bins

    % Saving out all data ----
    spec_outname = strcat('rs-', visit_name, '-', num2str(this_ss),...
                    '-spec-res.mat');
    save(fullfile(outpath, spec_outname),'spec_res'); % saves out as matlab struct
end