% DispEn calculation script

% Specify the number of subjects
num_subjects = 22; % For example, if there are 22 subjects

% Embedding dimension (m), Number of classes (nc), Mapping approach (MA), and Time delay (tau) settings
m = 2;                   % Embedding dimension
nc = 6;                  % Number of classes
MA = 'NCDF';           % Mapping approach
tau = 1;                 % Time delay

% Initialize combined results table
combined_results_pre = table();
combined_results_post = table();

% Folder paths for pre and post files
pre_dir = 'C:/Users/Ali/Downloads/itf_roiseries/pre/';
post_dir = 'C:/Users/Ali/Downloads/itf_roiseries/post/';

% Calculate DispEn for each subject
for subject_idx = 1:num_subjects
	    % Create pre file path
	        pre_file_path = fullfile(pre_dir, sprintf('Subject_%02d_ROIs_data_pre.csv', subject_idx));
		    
		    % Create post file path
		        post_file_path = fullfile(post_dir, sprintf('Subject_%02d_ROIs_data_post.csv', subject_idx));
			    
			    % Perform DispEn calculations for pre and post files
			        for file_path = {pre_file_path, post_file_path}
					        file_path = char(file_path);
						        
						        % Read the CSV file
							        data = readtable(file_path, 'VariableNamingRule', 'preserve');
								        column_names = data.Properties.VariableNames;
									        num_columns = width(data);
										        
										        % Sanitize column names
											        sanitized_column_names = matlab.lang.makeValidName(column_names);
												        
												        % Initialize output arrays
													        Out_DispEn_all = zeros(1, num_columns);
														        
														        % Define custom logsig function
															        logsig = @(x) 1 ./ (1 + exp(-x));
																        
																        % Define custom mapminmax function
																	        mapminmax_custom = @(x, ymin, ymax) ymin + (ymax - ymin) * (x - min(x)) / (max(x) - min(x));
																		        
																		        % Loop through each column to calculate DispEn
																			        for col_idx = 1:num_columns
																					            x = table2array(data(:, col_idx))';
																						                
																						                N = length(x);
																								            sigma_x = std(x);
																									                mu_x = mean(x);
																											        
																											            % Mapping approaches
																												                switch MA
																															                case 'LM'
																																		                    y = mapminmax_custom(x, 0, 1);
																																				                        y(y == 1) = 1 - eps;
																																							                    y(y == 0) = eps;
																																									                        z = round(y * nc + 0.5);
																																												        
																																												                case 'NCDF'
																																															                    y = normcdf(x, mu_x, sigma_x);
																																																	                        y = mapminmax_custom(y, 0, 1);
																																																				                    y(y == 1) = 1 - eps;
																																																						                        y(y == 0) = eps;
																																																									                    z = round(y * nc + 0.5);
																																																											            
																																																											                    case 'LOGSIG'
																																																														                        y = logsig((x - mu_x) / sigma_x);
																																																																	                    y = mapminmax_custom(y, 0, 1);
																																																																			                        y(y == 1) = 1 - eps;
																																																																						                    y(y == 0) = eps;
																																																																								                        z = round(y * nc + 0.5);
																																																																											        
																																																																											                case 'TANSIG'
																																																																														                    y = tansig((x - mu_x) / sigma_x) + 1;
																																																																																                        y = mapminmax_custom(y, 0, 1);
																																																																																			                    y(y == 1) = 1 - eps;
																																																																																					                        y(y == 0) = eps;
																																																																																								                    z = round(y * nc + 0.5);
																																																																																										            
																																																																																										                    case 'SORT'
																																																																																													                        N = length(x);
																																																																																																                    x = x(1:nc * floor(N / nc));
																																																																																																		                        [sx, osx] = sort(x);
																																																																																																					                    Fl_NC = N / nc;
																																																																																																							                        cx = [];
																																																																																																										                    for i = 1:nc
																																																																																																													                            cx = [cx, i * ones(1, Fl_NC)];
																																																																																																																                        end
																																																																																																																			                    for i = 1:N
																																																																																																																						                            z(i) = cx(osx == i);
																																																																																																																									                        end
																																																																																																																												            end
																																																																																																																													            
																																																																																																																													                % Generate all possible patterns
																																																																																																																															            all_patterns = [1:nc]';
																																																																																																																																                for f = 2:m
																																																																																																																																			                temp = all_patterns;
																																																																																																																																					                all_patterns = [];
																																																																																																																																							                j = 1;
																																																																																																																																									                for w = 1:nc
																																																																																																																																												                    [a, b] = size(temp);
																																																																																																																																														                        all_patterns(j:j + a - 1, :) = [temp, w * ones(a, 1)];
																																																																																																																																																	                    j = j + a;
																																																																																																																																																			                    end
																																																																																																																																																					                end
																																																																																																																																																							        
																																																																																																																																																							            % Calculate DispEn
																																																																																																																																																								                for i = 1:nc^m
																																																																																																																																																											                key(i) = 0;
																																																																																																																																																													                for ii = 1:m
																																																																																																																																																																                    key(i) = key(i) * 100 + all_patterns(i, ii);
																																																																																																																																																																		                    end
																																																																																																																																																																				                end
																																																																																																																																																																						        
																																																																																																																																																																						            embd2 = zeros(N - (m - 1) * tau, 1);
																																																																																																																																																																							                for i = 1:m
																																																																																																																																																																										                embd2 = [z(1 + (i - 1) * tau:N - (m - i) * tau)]' * 100^(m - i) + embd2;
																																																																																																																																																																												            end
																																																																																																																																																																													            
																																																																																																																																																																													                pdf = zeros(1, nc^m);
																																																																																																																																																																															        
																																																																																																																																																																															            for id = 1:nc^m
																																																																																																																																																																																	                    [R, C] = find(embd2 == key(id));
																																																																																																																																																																																			                    pdf(id) = length(R);
																																																																																																																																																																																					                end
																																																																																																																																																																																							        
																																																																																																																																																																																							            npdf = pdf / (N - (m - 1) * tau);
																																																																																																																																																																																								                p = npdf(npdf ~= 0);
																																																																																																																																																																																										            Out_DispEn = -sum(p .* log(p));
																																																																																																																																																																																											            
																																																																																																																																																																																											                % Store results in output arrays
																																																																																																																																																																																													            Out_DispEn_all(col_idx) = Out_DispEn;
																																																																																																																																																																																														            end
																																																																																																																																																																																															            
																																																																																																																																																																																															            % Create a table for the current subject's results
																																																																																																																																																																																																            subject_label = sprintf('Subject_%02d', subject_idx);
																																																																																																																																																																																																	            subject_results = array2table(Out_DispEn_all, 'VariableNames', sanitized_column_names, 'RowNames', {subject_label});
																																																																																																																																																																																																		            
																																																																																																																																																																																																		            % Append to the combined results table
																																																																																																																																																																																																			            if contains(file_path, 'pre')
																																																																																																																																																																																																					                combined_results_pre = [combined_results_pre; subject_results];
																																																																																																																																																																																																							            % Create output file name for the subject
																																																																																																																																																																																																								                output_file = sprintf('%s\\Subject_%02d_DispEn_results_pre.csv', pre_dir, subject_idx);
																																																																																																																																																																																																										        else
																																																																																																																																																																																																												            combined_results_post = [combined_results_post; subject_results];
																																																																																																																																																																																																													                % Create output file name for the subject
																																																																																																																																																																																																															            output_file = sprintf('%s\\Subject_%02d_DispEn_results_post.csv', post_dir, subject_idx);
																																																																																																																																																																																																																            end
																																																																																																																																																																																																																	            
																																																																																																																																																																																																																	            % Write results to CSV file
																																																																																																																																																																																																																		            writetable(subject_results, output_file, 'WriteRowNames', true);
																																																																																																																																																																																																																			            
																																																																																																																																																																																																																			            % Inform the user
																																																																																																																																																																																																																				            fprintf('DispEn calculated for Subject %02d and saved to "%s".\n', subject_idx, output_file);
																																																																																																																																																																																																																					        end
																																																																																																																																																																																																																					end

																																																																																																																																																																																																																					% Write combined results to files
																																																																																																																																																																																																																					combined_output_file_pre = 'Combined_DispEn_results_pre.csv';
																																																																																																																																																																																																																					writetable(combined_results_pre, combined_output_file_pre, 'WriteRowNames', true);

																																																																																																																																																																																																																					combined_output_file_post = 'Combined_DispEn_results_post.csv';
																																																																																																																																																																																																																					writetable(combined_results_post, combined_output_file_post, 'WriteRowNames', true);

																																																																																																																																																																																																																					fprintf('DispEn results for all subjects saved to "%s" for pre and "%s" for post.\n', combined_output_file_pre, combined_output_file_post);

