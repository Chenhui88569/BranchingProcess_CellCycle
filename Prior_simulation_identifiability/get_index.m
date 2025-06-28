function   [prioritized_samples, hist_MCMC_samples_folder_coll]   = get_index(history_folder, expected_samples )

hist_MCMC_samples_folder_coll = {};
file_info = [];  % Struct array to store visitation and data size info

rootdir_hist = history_folder;
allFiles = dir(rootdir_hist);
dirFlags = [allFiles.isdir];
subFolders = allFiles(dirFlags);
subFolders = subFolders(3:end);  % Skip '.' and '..'
visited_sample_names = {};

for i = 1:length(subFolders)
    full_subfolder = fullfile(rootdir_hist, subFolders(i).name);
    subSubDirs = dir(full_subfolder);
    subSubDirs = subSubDirs([subSubDirs.isdir]);
    subSubDirs = subSubDirs(3:end);  % Skip '.' and '..'

    for j = 1:length(subSubDirs)
        sample_name = subSubDirs(j).name;
        curr_sub_sub_folder = fullfile(full_subfolder, sample_name);
        mat_files = dir(fullfile(curr_sub_sub_folder, '*.mat'));

        if ~isempty(mat_files)
            full_mat_path = fullfile(curr_sub_sub_folder, mat_files(1).name);
            hist_MCMC_samples_folder_coll{end+1} = full_mat_path;
            visited_sample_names{end+1} = sample_name;

            % Load and inspect file
            try
                s = load(full_mat_path, 'Para_col_accepted');
                if isfield(s, 'Para_col_accepted')
                    nonzeroRows = all(s.Para_col_accepted,2);
                    Para_col_accepted_hist = s.Para_col_accepted(nonzeroRows, :);
                    num_data_points = size(  Para_col_accepted_hist, 1);
                else
                    num_data_points = 0;
                end
                visited_flag = 1;
            catch
                % File couldn't be loaded properly — mark as unvisited
                num_data_points = 0;
                visited_flag = 0;

            end
            file_info(end+1).sample =   sample_name;
            file_info(end).path = full_mat_path;
            file_info(end).visited = 1;
            file_info(end).num_points = num_data_points;

        end
    end
end


% Step 3: Sort by unvisited first, then increasing num_points
info_table = struct2table(file_info);

[unique_samples, ~, group_idx] = unique(info_table.sample, 'stable');

% Initialize logical index to keep the best row per group
keep_idx = false(height(info_table), 1);

% Loop over each group and keep only the row with max num_points
for i = 1:numel(unique_samples)
    group_rows = find(group_idx == i);
    [~, max_idx] = max(info_table.num_points(group_rows));
    keep_idx(group_rows(max_idx)) = true;
end

% Filter table
info_table = info_table(keep_idx, :);


[~, sort_idx] = sortrows([info_table.visited, info_table.num_points], [1 2]);
prioritized_samples = sort_idx;
end