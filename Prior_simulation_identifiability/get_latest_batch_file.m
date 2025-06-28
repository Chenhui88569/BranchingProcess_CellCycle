function latestFile = get_latest_batch_file(batch_number, Model_Index)
    % Construct the folder path
    folderPath = sprintf('Batch%d', batch_number);
    
    % Construct the filename pattern
    pattern = sprintf('prior_simulation_identifiability_batch%dModel%d_*', batch_number, Model_Index);

    % Get list of matching files in the folder
    files = dir(fullfile(folderPath, pattern));

    % Filter only files (ignore folders)
    files = files([files.isdir]);

    if isempty(files)
        warning('No matching files found in %s.', folderPath);
        latestFile = '';
        return;
    end

    % Find the file with the latest timestamp
    [~, idx] = max([files.datenum]);

    % Return the filename
    latestFile = files(idx).name;
    
    latestFile = fullfile(folderPath, latestFile)

    % Optionally display the result
    fprintf('Latest file: %s\n', latestFile);
end
