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

% For spectral decomposition settings
wsize = 8; % FFT window size in seconds
olap = 4; % FFT window overlap in seconds

% For peak alpha frequency (PAF) and center of gravity (COG) calculations
% see: https://github.com/corcorana/restingIAF/blob/master/tutorial/tutorial.m
cmin = 3;           % minimum number of channel estimates required for 
                    % cross-channel averages
fRange = [1 40];    % spectral range (set to filter passband)
w = [7.5 13];       % alpha peak search window (Hz); from McClain et al. 2022
Fw = 23;            % SGF frame width (Fw * freq resoltion = window)
poly = 5;           % SGF polynomial order (same as described in Corcoran et al., 2017)

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
    
    % Epoching
    % selecting the stimulation blocks (60 second epochs)
    blocks = {'S111' 'S102' 'S103' 'S114' 'S105' 'S116' 'S117' 'S108'};
    blocks_end = {'S211' 'S202' 'S203' 'S214' 'S205' 'S216' 'S217' 'S208'};
    
    % preallocates arrays
    N = EEG.srate*wsize; % number of data points in FFT window
    poss_freqs = N/2 + 1; % possible frequencies (freq bins = EEG.srate/N)
    this_spectra = zeros(EEG.nbchan, poss_freqs, length(blocks));
    this_freqs = zeros(poss_freqs, 1, length(blocks));
    this_paf = zeros(EEG.nbchan, length(blocks)); % peak alpha freq
    this_cog = zeros(EEG.nbchan, length(blocks)); % center of gravity
    
    % grand average PAF and COG (rows are blocks, cols are PAF and COG)
    this_iaf = zeros(length(blocks), 2); 
    
    % In order for this script to accomodate rejected epochs, the
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
    
    j=1; % for testing purposes
    
    for j = 1:length(blocks)
        
        this_epoch = [0 block_durations(j)]; % Sets the epoch duration
        
        try % This will run if the block exists
            % Selects blocks (in order)
            this_EEG = pop_epoch(EEG,blocks(j),this_epoch,'epochinfo','yes');
            % Spectral decomposition here
            [this_spectra(:,:,j), this_freqs(:,:,j)] = spectopo(...
                this_EEG.data(:,:), ... 
                0, ... % frames per epoch (default 0 = data length)
                this_EEG.srate, ... % sampling rate
                'winsize', wsize*this_EEG.srate, ... % window size
                'overlap', olap*this_EEG.srate, ... % overlap size
                'plot','off'... % toggles plot
                );
            
            % Peak alpha frequency (PAF) and center of gravity (COG)
            [pSum, pChans, f] = restingIAF(...
                this_EEG.data,... % EEG data
                this_EEG.nbchan,... % number of channels
                cmin,... % min number of channels to grand average PAF
                fRange,... % spectral range
                this_EEG.srate,... % sampling rate
                w,... % alpha peak search window
                Fw,... % window size
                poly,... % polynomial order
                'taper', 'hamming',... % type of taper (default)
                'tlen', N,...% so that window size remains the same as spectral decomp (50% overlap by default)
                'mdiff', .2); % peak must be 20% larger than surrounding (default)  
            
            % fills in preallocated matrices
            this_paf(:,j) = [pChans.peaks]; % peak alpha freq per chan
            this_cog(:,j) = [pChans.gravs]; % center of gravity per chan
            this_iaf(j,:) = [pSum.paf pSum.cog]; %IAF (weighted by Q; see Corcoran et al., 2017)
            
                      
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
    
    % Calculating pink and white noise via PaWNextra.m
        % See Barry & De Blasio (2021)
        this_spectra_t = permute(this_spectra, [2 1 3]); % transposes
        [iterations, PN, PN_slope, WN, n_sols, adjust_value, error_value, FL_chan] = ...
            PaWNextra(this_spectra_t,... % broadband spectra freqs x chans x participants)
            this_freqs(:,:,j),... % column vector of freq bins
            -0.001,... % error threshold (E=-.001 is the default)
            1000); % iterations (I believe this is the default)
    
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
    spec_res.paf       = this_paf;     % 2D mat of PAF per chan per block
    spec_res.cog       = this_cog;     % 2D mat of COG per chan per block
    spec_res.iaf       = this_iaf;     % 2D mat of IAF per block

    % Saving out all data ----
    spec_outname = strcat('rs-', visit_name, '-', num2str(this_ss),...
                    '-spec-res.mat');
    save(fullfile(outpath, spec_outname),'spec_res'); % saves out as matlab struct
end