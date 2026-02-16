



% Directory containing subdirectories with .set files
cfg.dir  = '/Users/dhewitt/OneDrive - Nexus365/Data/Touch/SAData24/';

% Path to the corrected channel locations file
corrected_chan_file = '/Users/dhewitt/OneDrive - Nexus365/Analysis/TouchStudy/gtec64new.sph';

% List all .set files in subdirectories
set_files = dir(fullfile(cfg.dir, '**', '*.set'));

% Filter the files to only include those ending with "cleaned.set"
cleaned_set_files = set_files(endsWith({set_files.name}, 'erd_cleaned_2204.set', 'IgnoreCase', true));

% Loop through each .set file
for i = 1:length(cleaned_set_files)
    filename = fullfile(cleaned_set_files(i).folder, cleaned_set_files(i).name);

    % Load .set file
    EEG = pop_loadset(filename);

    % Apply corrected channel locations
    %EEG = pop_chanedit(EEG, 'load', {corrected_chan_file, 'filetype', 'autodetect'}, 'nosedir','+Y','rplurchanloc',1);
    EEG = pop_chanedit(EEG, 'load', {corrected_chan_file, 'filetype', 'autodetect'});


    [file_dir, file_name, ~] = fileparts(filename);
    %new_filename = fullfile([file_name, '_corrected_1604_noy.set']);
    new_filename = fullfile([file_name, '_corrected.set']);
    new_file_dir = strrep(file_dir, 'Library/CloudStorage/OneDrive-Nexus365', 'OneDrive - Nexus365');

    EEG = pop_saveset(EEG, 'filename', new_filename, 'filepath', new_file_dir);

    disp(['Processed: ', filename]);
    disp(['Saved as: ', new_filename]);
end