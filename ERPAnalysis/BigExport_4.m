% Quick way to get all the ERP exports into a big export file for further
% analysis
% Last updated: Danielle Hewitt, 2nd Jan 2024

%% Set your paths and get the data
%=====================================================================
close all; clear all; d = char(datetime('today'));
folderPath = '/Users/dhewitt/Data/Touch/TouchStudy2/ERPAnalysis/Exports/';

files = dir(fullfile(folderPath, 'AllERP_csd_*.xlsx'));

% Initialize an empty table to store the combined data
combinedTable = table();

% Loop through each Excel file
for i = 1:length(files)
    % Extract information from the file name
    fileNameParts = strsplit(files(i).name, '_');
    electrode = fileNameParts{end-2}; % CPz
    timeStart = str2double(fileNameParts{end-1}); % 210
    timeEndsplit = strsplit(fileNameParts{end},'.');
    timeEnd = str2double(timeEndsplit{1}); % 310

    % Load the data from the Excel file
    data = readtable(fullfile(folderPath, files(i).name));

    % Add subheadings to the data
    data.Electrode = repmat({strrep(electrode, '-', '_')}, size(data, 1), 1);
    data.TimeStart = repmat(timeStart, size(data, 1), 1);
    data.TimeEnd = repmat(timeEnd, size(data, 1), 1);

    % Concatenate the data to the combined table
    combinedTable = [combinedTable; data];

end

 combinedTable.Properties.VariableNames(1:5) = {'ID','sOther', 'sSelf', 'fOther', 'fSelf'};

% Optionally, you can save the combined table to a new Excel file
writetable(combinedTable, [folderPath 'AllERPData_UK2_' num2str(d) '.xlsx']);
