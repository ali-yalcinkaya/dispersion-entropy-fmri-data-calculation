% DispEn calculation script for fMRI data
% This function calculates dispersion entropy (DispEn) of a univariate signal
% The code is adapted from:
% H. Azami and J. Escudero, "Amplitude- and Fluctuation-based Dispersion Entropy", Entropy, 2018.
% Available at: https://github.com/HamedAzami/Univariate_Entropy_Methods


% Extract BOLD denoised timeseries from SPM/CONN Toolbox outputs,

clearvars

%% ANALYSIS PARAMETERS
Info.wdir       = 'C:/Users/Ali/Downloads/itf_roiseries/'; 
Info.session    = 0; % 0 for all sessions
Info.nsub       = 22; % Number of subjects
Info.outdir     = pwd;

pre_dir  = fullfile(Info.outdir, 'pre');
post_dir = fullfile(Info.outdir, 'post');

if ~exist(pre_dir, 'dir'); mkdir(pre_dir); end
if ~exist(post_dir, 'dir'); mkdir(post_dir); end

fprintf('---------------------------------------------');
fprintf('\nANALYSIS INFO');
fprintf('\nSession:\t%d', Info.session);
fprintf('\nSubjects:\t%d', Info.nsub);
fprintf('\n---------------------------------------------\n');

%% CREATE CSV FILES FROM CONN OUTPUTS
for i = 1:Info.nsub
    matfile = [ 'ROI_Subject0', num2str(i, '%02i') ,'_Condition000.mat' ];
    load([Info.wdir matfile]);
    
    ROI.names = names;
    ROI.dsess = data_sessions;
    cond = find(ROI.dsess);
    
    max_length = 0;
    for j = 1:length(names)
        ROI_data = cell2mat(data(j));
        max_length = max(max_length, size(ROI_data,1));
    end
    
    subject_data = NaN(max_length, length(names));
    for j = 1:length(names)
        ROI_data = cell2mat(data(j));
        subject_data(1:size(ROI_data,1), j) = ROI_data(:,1);
    end
    
    cleaned_names = regexprep(names, '[\s\(\),]', '_');
    
    pre_data = subject_data(1:min(255, size(subject_data,1)), :);
    post_data = subject_data(max(1, size(subject_data,1)-254):end, :);
    
    csv_header = strjoin(cleaned_names, ',');
    csv_filename_pre = sprintf('%s/Subject_%02d_ROIs_data_pre.csv', pre_dir, i);
    csv_filename_post = sprintf('%s/Subject_%02d_ROIs_data_post.csv', post_dir, i);

    csvwrite_with_headers(csv_filename_pre, pre_data, csv_header);
    csvwrite_with_headers(csv_filename_post, post_data, csv_header);

    fprintf('Subject %02d pre/post CSV files saved.\n', i);
end

%% DispEn parameters (User can choose: 'LM', 'NCDF', 'LOGSIG', 'TANSIG', 'SORT')
m = 2;
nc = 6;
MA = 'NCDF';
tau = 1;

combined_results_pre = table();
combined_results_post = table();

for subject_idx = 1:Info.nsub
    pre_file_path = fullfile(pre_dir, sprintf('Subject_%02d_ROIs_data_pre.csv', subject_idx));
    post_file_path = fullfile(post_dir, sprintf('Subject_%02d_ROIs_data_post.csv', subject_idx));

    for file_path = {pre_file_path, post_file_path}
        file_path = char(file_path);
        data = readtable(file_path, 'VariableNamingRule', 'preserve');
        column_names = matlab.lang.makeValidName(data.Properties.VariableNames);

        Out_DispEn_all = zeros(1, width(data));

        for col_idx = 1:width(data)
            x = table2array(data(:, col_idx))';
            sigma_x = std(x); mu_x = mean(x);
            logsig = @(x) 1 ./ (1 + exp(-x));
            mapminmax_custom = @(x, ymin, ymax) ymin + (ymax - ymin) * (x - min(x)) / (max(x) - min(x));

            switch MA
                case 'LM'
                    y = mapminmax_custom(x, 0, 1);
                case 'NCDF'
                    y = normcdf(x, mu_x, sigma_x);
                case 'LOGSIG'
                    y = logsig((x - mu_x) / sigma_x);
                case 'TANSIG'
                    y = tansig((x - mu_x) / sigma_x) + 1;
                case 'SORT'
                    N = length(x);
                    x = x(1:nc * floor(N / nc));
                    [~, osx] = sort(x);
                    Fl_NC = N / nc;
                    cx = repelem(1:nc, Fl_NC);
                    z = zeros(1, N);
                    z(osx) = cx;
                    goto_entropy = true;
                otherwise
                    error('Unknown mapping approach');
            end

            if ~exist('goto_entropy','var')
                y = mapminmax_custom(y, 0, 1);
                y(y == 1) = 1 - eps; y(y == 0) = eps;
                z = round(y * nc + 0.5);
            else
                clear goto_entropy
            end

            all_patterns = (1:nc)';
            for f = 2:m
                temp = all_patterns;
                all_patterns = [];
                j = 1;
                for w = 1:nc
                    all_patterns(j:j+size(temp,1)-1,:) = [temp,w*ones(size(temp,1),1)];
                    j = j+size(temp,1);
                end
            end
            key = sum(all_patterns .* (100.^(m-1:-1:0)), 2)';
            N = length(z);
            embd2 = zeros(N-(m-1)*tau,1);
            for emb_idx = 1:m
                embd2 = embd2 + z((1+(emb_idx-1)*tau):(N-(m-emb_idx)*tau))' * 100^(m-emb_idx);
            end
            pdf = histcounts(embd2, [key, max(key)+1]);
            npdf = pdf/sum(pdf); p = npdf(npdf>0);
            Out_DispEn_all(col_idx) = -sum(p.*log(p));
        end

        subject_label = sprintf('Subject_%02d', subject_idx);
        subject_results = array2table(Out_DispEn_all, 'VariableNames', column_names, 'RowNames', {subject_label});

        if contains(file_path, 'pre')
            combined_results_pre = [combined_results_pre; subject_results];
            writetable(subject_results, sprintf('%s/Subject_%02d_DispEn_results_pre.csv', pre_dir, subject_idx), 'WriteRowNames',true);
        else
            combined_results_post = [combined_results_post; subject_results];
            writetable(subject_results, sprintf('%s/Subject_%02d_DispEn_results_post.csv', post_dir, subject_idx), 'WriteRowNames',true);
        end

        fprintf('DispEn calculated for Subject %02d (%s).\n', subject_idx, file_path);
    end
end

writetable(combined_results_pre, 'Combined_DispEn_results_pre.csv', 'WriteRowNames',true);
writetable(combined_results_post, 'Combined_DispEn_results_post.csv', 'WriteRowNames',true);

fprintf('All DispEn results saved as combined CSV files.\n');

%% Supporting function: CSV Write with headers
function csvwrite_with_headers(filename, data, headers)
    fid = fopen(filename, 'w');
    fprintf(fid, '%s\n', headers);
    fclose(fid);
    dlmwrite(filename, data, '-append');
end
