% main.m
% ADNI Data WMH Spatial Heterogeneity Analysis Project
% Note: SPM12 must be installed and added to the MATLAB path before running.

clc; clear; close all;

%% 0. Setup Paths and Validate Environment
if exist('spm', 'file') == 0
    error('SPM12 is not installed or not in the MATLAB path.');
end

% Initialize SPM in fMRI defaults
spm('defaults', 'fmri');
spm_jobman('initcfg');

adni_root = fullfile(pwd, 'images');
output_dir = fullfile(pwd, 'processed');
if ~exist(output_dir, 'dir'), mkdir(output_dir); end

% Get valid subject directories (ignore hidden files/folders)
d = dir(adni_root);
subjects = d([d.isdir] & ~startsWith({d.name}, '.') & ~startsWith({d.name}, '__'));
num_subjects = length(subjects);

if num_subjects == 0
    error('No subject folders found in %s', adni_root);
end

% Locate MNI template for normalization
mni_template = fullfile(pwd, 'MNI152_T1_1mm.nii', 'MNI152_T1_1mm.nii');
if ~exist(mni_template, 'file') || isdir(mni_template)
    mni_template = fullfile(spm('dir'), 'tpm', 'TPM.nii');
end

%% 1. Image Registration (Coregistration & Normalization)
disp('Starting image registration (SPM12)...');
subject_ids = cell(num_subjects, 1);
registered_masks = cell(num_subjects, 1);
ref_vol_size = []; 

for i = 1:num_subjects
    sub_id = subjects(i).name;
    subject_ids{i} = sub_id;
    sub_dir = fullfile(adni_root, sub_id, 'nifti');
    out_sub_dir = fullfile(output_dir, sub_id);
    if ~exist(out_sub_dir, 'dir'), mkdir(out_sub_dir); end
    
    fprintf('Processing Subject %d/%d: %s\n', i, num_subjects, sub_id);
    
    % Locate required NIfTI files
    t1_gz = fullfile(sub_dir, 'T1_brain.nii.gz');
    if ~exist(t1_gz, 'file'), t1_gz = fullfile(sub_dir, 'T1.nii.gz'); end
    t2_gz = fullfile(sub_dir, 'T2_FLAIR.nii.gz');
    mask_gz = fullfile(sub_dir, 'mask_in_T2.nii.gz');
    
    if ~exist(mask_gz, 'file')
        error('WMH mask not found for subject %s', sub_id);
    end
    
    % Unzip to the processed directory for SPM compatibility
    t1_nii = fullfile(out_sub_dir, 'T1.nii');
    t2_nii = fullfile(out_sub_dir, 'T2.nii');
    mask_nii = fullfile(out_sub_dir, 'mask.nii');
    
    unzip_to_nii(t1_gz, out_sub_dir, 'T1.nii');
    unzip_to_nii(t2_gz, out_sub_dir, 'T2.nii');
    unzip_to_nii(mask_gz, out_sub_dir, 'mask.nii');
    
    % Step 1.1: Rigid Coregister T2 -> T1
    % We use Normalized Mutual Information (NMI) as the cost function 
    % because T1 and T2 FLAIR are different modalities.
    if ~exist(fullfile(out_sub_dir, 'rmask.nii'), 'file')
        clear matlabbatch;
        matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {[t1_nii ',1']};
        matlabbatch{1}.spm.spatial.coreg.estwrite.source = {[t2_nii ',1']};
        matlabbatch{1}.spm.spatial.coreg.estwrite.other = {[mask_nii ',1']};
        matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
        matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
        matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
        matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
        % Note: Ideally interpolation for binary mask should be Nearest Neighbor (0)
        % Using 1 (linear) here requires subsequent binarization (threshold > 0.5)
        matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 1;
        matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
        matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
        matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';
        spm_jobman('run', matlabbatch);
    end
    
    % Step 1.2: Non-Rigid Normalize T1 -> MNI Template
    % Applies non-linear deformations to map individual anatomy to standard space
    sn_mat_file = fullfile(out_sub_dir, 'T1_sn.mat');
    if ~exist(fullfile(out_sub_dir, 'wrmask.nii'), 'file')
        clear matlabbatch;
        if exist(sn_mat_file, 'file')
            % Skip estimation if deformation matrix already exists
            matlabbatch{1}.spm.tools.oldnorm.write.subj.matname = {sn_mat_file};
            matlabbatch{1}.spm.tools.oldnorm.write.subj.resample = {[fullfile(out_sub_dir, 'rmask.nii') ',1']};
            matlabbatch{1}.spm.tools.oldnorm.write.roptions.preserve = 0;
            matlabbatch{1}.spm.tools.oldnorm.write.roptions.bb = [-78 -112 -70; 78 76 85];
            matlabbatch{1}.spm.tools.oldnorm.write.roptions.vox = [2 2 2];
            matlabbatch{1}.spm.tools.oldnorm.write.roptions.interp = 1; 
            matlabbatch{1}.spm.tools.oldnorm.write.roptions.wrap = [0 0 0];
            matlabbatch{1}.spm.tools.oldnorm.write.roptions.prefix = 'w';
        else
            matlabbatch{1}.spm.tools.oldnorm.estwrite.subj.source = {[t1_nii ',1']};
            matlabbatch{1}.spm.tools.oldnorm.estwrite.subj.wtsrc = '';
            matlabbatch{1}.spm.tools.oldnorm.estwrite.subj.resample = {[fullfile(out_sub_dir, 'rmask.nii') ',1']};
            matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.template = {[mni_template ',1']};
            matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.weight = '';
            matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.smosrc = 8;
            matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.smoref = 0;
            matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.regtype = 'mni';
            matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.cutoff = 25;
            matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.nits = 16;
            matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.reg = 1;
            matlabbatch{1}.spm.tools.oldnorm.estwrite.roptions.preserve = 0;
            matlabbatch{1}.spm.tools.oldnorm.estwrite.roptions.bb = [-78 -112 -70; 78 76 85];
            matlabbatch{1}.spm.tools.oldnorm.estwrite.roptions.vox = [2 2 2];
            matlabbatch{1}.spm.tools.oldnorm.estwrite.roptions.interp = 1;
            matlabbatch{1}.spm.tools.oldnorm.estwrite.roptions.wrap = [0 0 0];
            matlabbatch{1}.spm.tools.oldnorm.estwrite.roptions.prefix = 'w';
        end
        spm_jobman('run', matlabbatch);
    end
    
    % Step 1.3: Post-processing (Binarize, Smooth and Downsample)
    V_mask = spm_vol(fullfile(out_sub_dir, 'wrmask.nii'));
    mask_vol = spm_read_vols(V_mask) > 0.5; % Fix interpolation artifacts by thresholding
    
    % 3D Gaussian smooth to reduce noise, then downsample by factor of 2
    smoothed_mask = imgaussfilt3(double(mask_vol), 1);
    downsampled_mask = smoothed_mask(1:2:end, 1:2:end, 1:2:end); 
    
    if isempty(ref_vol_size), ref_vol_size = size(downsampled_mask); end
    registered_masks{i} = downsampled_mask(:);
end

% Flatten subject masks into a voxel-by-subject matrix for clustering
valid_subs = ~cellfun(@isempty, registered_masks);
registered_masks = registered_masks(valid_subs);
num_valid = sum(valid_subs);
L_full = cat(2, registered_masks{:}); 

%% 2. Generate Figure 1 (Registration Results)
disp('Generating Figure 1...');
% Pick 3 representative subjects with the highest WMH volume for clear display
total_lesion_per_sub = sum(L_full, 1);
[~, sorted_idx] = sort(total_lesion_per_sub, 'descend');
rep_idx = sorted_idx(1:min(3, num_valid)); 
num_rep = length(rep_idx);

fig1 = figure('Name', 'Registration Results', 'Position', [50, 50, 1400, 900], 'Color', 'k');

for row = 1:num_rep
    si = rep_idx(row);
    sub_id = subject_ids{si};
    sub_proc_dir = fullfile(output_dir, sub_id);
    
    vol_t1 = spm_read_vols(spm_vol(fullfile(sub_proc_dir, 'T1.nii')));
    vol_t2 = spm_read_vols(spm_vol(fullfile(sub_proc_dir, 'T2.nii')));
    vol_rt2 = spm_read_vols(spm_vol(fullfile(sub_proc_dir, 'rT2.nii')));
    vol_rmask = spm_read_vols(spm_vol(fullfile(sub_proc_dir, 'rmask.nii')));
    
    % Automatically select the axial slice containing the most WMH lesions
    mask_sum = squeeze(sum(sum(vol_rmask, 1), 2));
    [~, best_slice] = max(mask_sum);
    if max(mask_sum) == 0, best_slice = round(size(vol_t1, 3) / 2); end
    
    slice_t1 = rot90(vol_t1(:,:,best_slice));
    slice_t2 = rot90(vol_t2(:,:,min(best_slice, size(vol_t2,3))));
    slice_rt2 = rot90(vol_rt2(:,:,best_slice));
    slice_rmask = rot90(vol_rmask(:,:,best_slice));
    
    slice_t1 = slice_t1 / max(slice_t1(:) + eps);
    slice_t2 = slice_t2 / max(slice_t2(:) + eps);
    slice_rt2 = slice_rt2 / max(slice_rt2(:) + eps);
    
    subplot(num_rep, 4, (row-1)*4 + 1);
    imagesc(slice_t1); colormap(gca, 'gray'); axis image off;
    if row == 1, title('T1 (Reference)', 'Color', 'w', 'FontSize', 12); end
    ylabel(sub_id, 'Color', 'w', 'FontSize', 10, 'Interpreter', 'none', 'Rotation', 0, 'HorizontalAlignment', 'right');
    
    subplot(num_rep, 4, (row-1)*4 + 2);
    imagesc(slice_t2); colormap(gca, 'gray'); axis image off;
    if row == 1, title('T2 FLAIR (Original)', 'Color', 'w', 'FontSize', 12); end
    
    subplot(num_rep, 4, (row-1)*4 + 3);
    imagesc(slice_rt2); colormap(gca, 'gray'); axis image off;
    if row == 1, title('T2 FLAIR (Coregistered)', 'Color', 'w', 'FontSize', 12); end
    
    % Overlay the co-registered mask in red channel over the T1 reference
    subplot(num_rep, 4, (row-1)*4 + 4);
    rgb_img = repmat(slice_t1, [1, 1, 3]);
    mask_overlay = slice_rmask > 0.5;
    rgb_img(:,:,1) = rgb_img(:,:,1) + 0.5 * double(mask_overlay);
    rgb_img(:,:,2) = rgb_img(:,:,2) .* (1 - 0.3 * double(mask_overlay));
    rgb_img(:,:,3) = rgb_img(:,:,3) .* (1 - 0.3 * double(mask_overlay));
    imagesc(min(rgb_img, 1)); axis image off;
    if row == 1, title('T1 + WMH Mask', 'Color', 'w', 'FontSize', 12); end
end
sgtitle('Figure 1: Representative Registration Results', 'Color', 'w', 'FontSize', 14, 'FontWeight', 'bold');
saveas(fig1, fullfile(pwd, 'Fig1_Registration_Results.png'));

%% 3. Clustering WMH Voxels (Dimensionality Reduction)
disp('Clustering voxels...');
% Filter out voxels that are less than 1% prevalent in the population
mean_WMH = mean(L_full, 2);
valid_voxels_idx = find(mean_WMH >= 0.01);
L = L_full(valid_voxels_idx, :);

% Use K-Means (r=4) for its efficiency and reasonable performance on 
% continuous probabilistic data resulting from the Gaussian smoothing.
r = 4; 
rng(42); % For reproducibility
idx = kmeans(L, r, 'Distance', 'sqeuclidean', 'Replicates', 10);

% Map 1D cluster indices back to 3D volume space
cluster_vol = zeros(prod(ref_vol_size), 1);
cluster_vol(valid_voxels_idx) = idx;
cluster_vol_3d = reshape(cluster_vol, ref_vol_size);

% Calculate the proportion of WMH within each region for each individual
lesion_percentage = zeros(num_valid, r);
for k = 1:r
    region_mask = (idx == k);
    region_size = sum(region_mask);
    for i = 1:num_valid
        lesion_percentage(i, k) = sum(L(region_mask, i)) / region_size;
    end
end

%% 4. Generate Figure 2 (Clustered Regions Visualization)
disp('Generating Figure 2...');
% Load a reference T1 image to act as background for the cluster overlay
wt1_file = fullfile(output_dir, subject_ids{1}, 'wT1.nii');
vol_wt1 = spm_read_vols(spm_vol(wt1_file));
vol_wt1_ds = vol_wt1(1:2:end, 1:2:end, 1:2:end);

sz = min([size(vol_wt1_ds); ref_vol_size]);
vol_wt1_ds = vol_wt1_ds(1:sz(1), 1:sz(2), 1:sz(3));
cluster_display = cluster_vol_3d(1:sz(1), 1:sz(2), 1:sz(3));
vol_wt1_ds = vol_wt1_ds / max(vol_wt1_ds(:) + eps);

cluster_colors = [1.0, 0.2, 0.2; 0.2, 0.8, 0.2; 0.3, 0.3, 1.0; 1.0, 0.8, 0.0];
active_slices = find(squeeze(sum(sum(cluster_display > 0, 1), 2)) > 0);
if length(active_slices) >= 6
    slice_indices = active_slices(round(linspace(1, length(active_slices), 6)));
else
    slice_indices = active_slices;
    while length(slice_indices) < 6, slice_indices = [slice_indices, round(sz(3)/2)]; end
end

fig2 = figure('Name', 'Clustered WMH Regions', 'Position', [50, 50, 1400, 700], 'Color', 'k');
for s = 1:6
    subplot(2, 3, s);
    sl = slice_indices(s);
    bg = rot90(vol_wt1_ds(:,:,sl));
    rgb = repmat(bg, [1, 1, 3]);
    cl = rot90(cluster_display(:,:,sl));
    
    alpha = 0.6;
    for k = 1:r
        region_mask = (cl == k);
        if any(region_mask(:))
            for ch = 1:3
                rgb(:,:,ch) = rgb(:,:,ch) .* (1 - alpha * double(region_mask)) + cluster_colors(k, ch) * alpha * double(region_mask);
            end
        end
    end
    imagesc(min(max(rgb, 0), 1)); axis image off;
    title(sprintf('Axial Slice %d', sl), 'Color', 'w', 'FontSize', 11);
end
sgtitle('Figure 2: Clustered WMH Dominant Regions in MNI Space', 'Color', 'w', 'FontSize', 14, 'FontWeight', 'bold');
saveas(fig2, fullfile(pwd, 'Fig2_Clustered_WMH_Regions.png'));

%% 5. Phenotypic Data Association (Figure 3)
disp('Correlating with phenotypic data...');
pheno_file = fullfile(pwd, 'all30m.xlsx');
if exist(pheno_file, 'file')
    pheno_table = readtable(pheno_file);
    pheno_vars = {'Age', 'Gender', 'APOE4', 'AV45', 'ADAS11'};
    var_cols = zeros(1, length(pheno_vars));
    
    % Match target variables with actual excel column names
    col_names = pheno_table.Properties.VariableNames;
    for v = 1:length(pheno_vars)
        match_idx = find(strcmpi(col_names, pheno_vars{v}));
        if isempty(match_idx) && strcmpi(pheno_vars{v}, 'Gender'), match_idx = find(strcmpi(col_names, 'PTGENDER')); end
        if isempty(match_idx) && strcmpi(pheno_vars{v}, 'Age'), match_idx = find(strcmpi(col_names, 'AGE')); end
        if ~isempty(match_idx), var_cols(v) = match_idx(1); end
    end
    
    % Extract data matching subject IDs
    X_data = NaN(num_valid, length(pheno_vars));
    for i = 1:num_valid
        ptid_col = find(strcmpi(col_names, 'PTID') | strcmpi(col_names, 'Subject'));
        if ~isempty(ptid_col)
            row_idx = find(strcmp(pheno_table{:, ptid_col(1)}, subject_ids{i}));
            if ~isempty(row_idx)
                row_idx = row_idx(1);
                for v = 1:length(pheno_vars)
                    if var_cols(v) > 0
                        val = pheno_table{row_idx, var_cols(v)};
                        if isnumeric(val), X_data(i, v) = val;
                        elseif iscellstr(val) || isstring(val)
                            if strcmpi(val, 'M') || strcmpi(val, 'Male'), X_data(i, v) = 0;
                            elseif strcmpi(val, 'F') || strcmpi(val, 'Female'), X_data(i, v) = 1;
                            else, X_data(i, v) = str2double(val);
                            end
                        end
                    end
                end
            end
        end
    end
    
    np = sum(var_cols > 0);
    found_vars = pheno_vars(var_cols > 0);
    X_data = X_data(:, var_cols > 0);
    
    % Use modern tiledlayout for automatic margin control
    fig_assoc = figure('Name', 'WMH Regional Associations', 'Position', [50, 50, 1600, 1000], 'Color', 'w');
    t = tiledlayout(np, r, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    for v = 1:np
        for k = 1:r
            ax = nexttile;
            x = X_data(:, v);
            y = lesion_percentage(:, k);
            
            valid = ~isnan(x) & ~isnan(y);
            if sum(valid) > 2
                [R, P] = corrcoef(x(valid), y(valid));
                r_val = R(1,2); p_val = P(1,2);
                
                % Highlight significant correlations (p < 0.05)
                if p_val < 0.05
                    dot_color = [0.85, 0.15, 0.15]; edge_color = [0.7, 0, 0];
                else
                    dot_color = [0.25, 0.55, 0.85]; edge_color = [0.1, 0.3, 0.6];
                end
                
                scatter(x(valid), y(valid), 40, dot_color, 'filled', 'MarkerEdgeColor', edge_color, 'LineWidth', 0.5);
                hold on;
                
                % Add line of best fit
                coeffs = polyfit(x(valid), y(valid), 1);
                x_fit = linspace(min(x(valid)), max(x(valid)), 100);
                plot(x_fit, polyval(coeffs, x_fit), '-', 'Color', [0.3, 0.3, 0.3], 'LineWidth', 1.2);
                hold off;
                
                if p_val < 0.05
                    title(sprintf('%s - Region %d *', found_vars{v}, k), 'FontSize', 11, 'FontWeight', 'bold', 'Color', 'k');
                    ax.XColor = [0.8, 0.15, 0.15]; ax.YColor = [0.8, 0.15, 0.15]; ax.LineWidth = 1.5;
                    stat_str = sprintf('r = %.2f\np = %.3f *', r_val, p_val);
                    text_bg = [1, 0.95, 0.95]; text_edge = [0.8, 0.2, 0.2];
                else
                    title(sprintf('%s - Region %d', found_vars{v}, k), 'FontSize', 11, 'FontWeight', 'bold', 'Color', 'k');
                    ax.XColor = [0.15, 0.15, 0.15]; ax.YColor = [0.15, 0.15, 0.15]; ax.LineWidth = 0.5;
                    stat_str = sprintf('r = %.2f\np = %.3f', r_val, p_val);
                    text_bg = [0.95, 0.95, 1]; text_edge = [0.6, 0.6, 0.8];
                end
                
                xlabel(found_vars{v}, 'FontSize', 10, 'Color', 'k');
                if k == 1, ylabel('Lesion %', 'FontSize', 10, 'Color', 'k'); end
                
                text(0.95, 0.95, stat_str, 'Units', 'normalized', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
                     'FontSize', 9, 'FontWeight', 'bold', 'BackgroundColor', text_bg, 'EdgeColor', text_edge, 'LineWidth', 1, 'Margin', 3);
                 
                % Tight bounds for each subplot ensures data points utilize full vertical space
                y_valid = y(valid);
                y_top = max(y_valid) * 1.20;
                if y_top == 0, y_top = 0.01; end
                ylim([0, y_top]);
                ax.Box = 'on';
            else
                title(sprintf('%s - Region %d\n(Insuff. data)', found_vars{v}, k), 'FontSize', 10, 'FontWeight', 'bold');
            end
            set(ax, 'FontSize', 9); grid on; ax.GridAlpha = 0.15;
        end
    end
    title(t, 'WMH Regional Lesion Proportion vs. Phenotypic Variables', 'FontSize', 16, 'FontWeight', 'bold');
    exportgraphics(fig_assoc, fullfile(pwd, 'WMH_Regional_Associations.png'), 'Resolution', 300);
end

disp('Processing complete.');

%% Helper Function
function unzip_to_nii(gz_file, dest_dir, new_name)
    % Extracts a .gz file to a specific destination with a specific name if it doesn't already exist
    if ~exist(fullfile(dest_dir, new_name), 'file')
        out_files = gunzip(gz_file, dest_dir);
        unzipped_file = out_files{1};
        expected_file = fullfile(dest_dir, new_name);
        if ~strcmp(unzipped_file, expected_file)
            movefile(unzipped_file, expected_file);
        end
    end
end
