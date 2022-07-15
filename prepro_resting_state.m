% EEG Preprocessing Pipeline Step 1
% Matt Kmiecik
% Started 15 July 2022

workspace_prep % Prepares workspace (see src/...)

subjs = string({RAW{2:size(RAW,1),1}}); %{'324'}; % Initializes subjects for batch processing (if applicable)
visit = 'assessment-visit-1'; % name of the folder for visit number

i=1; % for testing purposes

% Preprocessing ----
for i = 1:length(subjs)

    % Creating variables ----
    this_ss = NUM(i);
    this_ss_path = dir(fullfile(raw_data_path, visit, num2str(this_ss), 'asy*.vhdr'));
    this_ss_name = this_ss_path.name;

    % Loads in raw data using loadbv() from BrainVision plugin ----
    EEG = pop_loadbv(this_ss_path.folder, this_ss_name, [], []);
    [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'setname',this_ss_name,'gui','off');
    eeglab redraw
    
    % Some participants were accidentally recorded with 64 channels
    % Therefore, these channels need to be deleted for those subjects
    if EEG.nbchan == 65
        disp('65 channels detected...deleting extra channels...');
        % Deleting channels
        EEG = pop_select( EEG, 'nochannel',{'Ch33' 'Ch34' 'Ch35' 'Ch36' 'Ch37' 'Ch38' 'Ch39' 'Ch40' 'Ch41' 'Ch42' 'Ch43' 'Ch44' 'Ch45' 'Ch46' 'Ch47' 'Ch48' 'Ch49' 'Ch50' 'Ch51' 'Ch52' 'Ch53' 'Ch54' 'Ch55' 'Ch56' 'Ch57' 'Ch58' 'Ch59' 'Ch60' 'Ch61' 'Ch62' 'Ch63' 'Ch64'});
        [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'overwrite','on','gui','off'); 
    else
        disp('< 65 channels detected...carry on...');
    end
    
    % Some participants were recorded with a pressure bulb
    % Therefore, this channel will be removed
    if EEG.nbchan == 34
        disp('34 channels detected...deleting pressure bulb...');
        % Deleting pressure bulb channel
        EEG = pop_select( EEG, 'nochannel',{'Pressure Bulb'});
        [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'overwrite','on','gui','off');
    else
        disp('pressure bulb not detected...carry on...');
    end

    % Configuring channel locations ----
    % TP9 is the online reference electrode (channel 10) but is not 
    % included in recording
    chanlocs = loadbvef(chan_loc_path); % with TP9

    % Enters the channel locations for all except TP9
    EEG.chanlocs = chanlocs([1:33]); % All but TP9

    % Adds TP9 as channel, sets it as reference, and sets its location manually
    EEG = pop_chanedit(EEG, 'append',33,...
        'changefield',{34 'labels' 'TP9'},...
        'changefield',{34 'sph_radius' '1'},...
        'changefield',{34 'sph_theta' '108'},...
        'changefield',{34 'sph_phi' '-23'},...
        'changefield',{34 'Z' '-0.3907'},...
        'changefield',{34 'Y' '0.8755'},...
        'changefield',{34 'X' '-0.2845'},...
        'changefield',{34 'radius' '0.6278'},...
        'changefield',{34 'theta' '-108'},...
        'setref',{'1:31' 'TP9'});

    % Re-referencing ----
    % Per: https://sccn.ucsd.edu/wiki/I.4:_Preprocessing_Tools
    % First step is to compute common average reference:
    EEG = pop_reref(EEG, [],'refloc',...
        struct('labels',{'TP9'},...
        'sph_radius',{1},...
        'sph_theta',{108},...
        'sph_phi',{-23},...
        'theta',{-108},...
        'radius',{0.6278},...
        'X',{-0.2845},...
        'Y',{0.8755},...
        'Z',{-0.3907},...
        'type',{''},...
        'ref',{''},...
        'urchan',{[]},...
        'datachan',{0})); 

    % Next compute linked average mastoid
    EEG = pop_reref(EEG, [20 34] ,'keepref','on');

    % Removes Photo channel and EOG ----
    % EOG is removed because it was recorded with a different reference
    % than the cap electrodes and cannot be used in any meaningful analysis
    EEG = pop_select(EEG, 'nochannel',{'Photo', 'EOG'});

    % Downsamples to 256Hz
    EEG = pop_resample( EEG, 256);

    % Removing DC offset by subtracting the mean signal from each electrode
    EEG = pop_rmbase( EEG, [],[]);
    
    % Highpass filter at 1Hz (-6dB @ 1Hz, 425 point highpass, 2Hz transition band width)
    EEG = pop_eegfiltnew(EEG, 'locutoff',2,'plotfreqz',0);
    [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'overwrite','on','gui','off');

    % Cleanline ----
    % Removing electrical line noise @ 60 Hz
    EEG = pop_cleanline(EEG, 'bandwidth',2,'chanlist',[1:32] ,'computepower',1,'linefreqs',60,'normSpectrum',0,'p',0.01,'pad',2,'plotfigures',0,'scanforlines',1,'sigtype','Channels','tau',100,'verb',1,'winsize',4,'winstep',1);
    [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'setname',this_subj.id,'overwrite','on','gui','off'); 

    % Saves out preprocessed data for inspection
    EEG = eeg_checkset( EEG );
    EEG = pop_saveset( EEG, 'filename', outname, 'filepath' ,outpath);
    [ALLEEG, EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

end
