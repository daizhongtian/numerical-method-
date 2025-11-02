function complete_solution_no_simplifications()
    
    
    clear; close all; clc;
    
   
    load('typical_gait.mat');  
   
    
    n = size(lfoo, 1);    
    dt = mean(diff(time));
    fs = 1/dt;           
    
    %% Compute Noise-Free Left Ankle Dorsiflexion Angle α_LA
    LA_alpha_clean_rad = zeros(n,1);
    for i = 1:n
        %  (28)~(30)
        [x_LT, y_LT, z_LT] = construct_LT_coords(ltio(i,:), lank(i,:), lfeo(i,:));
        %  (34)~(36)
        [x_LP, y_LP, z_LP] = construct_LP_coords(lfoo(i,:), lfop(i,:), lfol(i,:));
        % formulas (13) and (15)
        LA_alpha_clean_rad(i) = compute_alpha_LA_strict(x_LP, y_LP, z_LT);
    end
    LA_alpha_clean_deg = rad2deg(LA_alpha_clean_rad);
    
    figure('Name','Noise-Free LA_\alpha (Degrees)');
    plot(time, LA_alpha_clean_deg, 'g', 'LineWidth', 1.5);
    grid on;
    xlabel('Time (s)');
    ylabel('Angle (degrees)');
    title('Left Ankle Dorsiflexion Angle (Noise-Free)');
    
   % Add Noise and Compute Unfiltered Angle 
    noise_std = 15;  
    ltio_noisy = ltio + noise_std * randn(size(ltio));
    lank_noisy = lank + noise_std * randn(size(lank));
    lfeo_noisy = lfeo + noise_std * randn(size(lfeo));
    lfoo_noisy = lfoo + noise_std * randn(size(lfoo));
    lfop_noisy = lfop + noise_std * randn(size(lfop));
    lfol_noisy = lfol + noise_std * randn(size(lfol));
    
    LA_alpha_noisy_deg = zeros(n,1);
    for i = 1:n
        [x_LT, y_LT, z_LT] = construct_LT_coords(ltio_noisy(i,:), lank_noisy(i,:), lfeo_noisy(i,:));
        [x_LP, y_LP, z_LP] = construct_LP_coords(lfoo_noisy(i,:), lfop_noisy(i,:), lfol_noisy(i,:));
        alpha_rad = compute_alpha_LA_strict(x_LP, y_LP, z_LT);
        LA_alpha_noisy_deg(i) = rad2deg(alpha_rad);
    end
    
    figure('Name','Noisy, Unfiltered LA_\alpha (Degrees)');
    plot(time, LA_alpha_noisy_deg, 'r', 'LineWidth', 1.5);
    grid on;
    xlabel('Time (s)');
    ylabel('Angle (degrees)');
    title('Left Ankle Dorsiflexion Angle (Noisy, Unfiltered)');
    
    %%  Single Filtering: Try Different Cutoff Frequencies and Compute MAE
    fc_values = [0.1, 0.5, 1, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5 ,6, 6.5, 10, 15, 20];
    MAE_values = zeros(size(fc_values));
    
    for iCutoff = 1:length(fc_values)
        fc = fc_values(iCutoff);
        
        ltio_filt = zeros(size(ltio_noisy));
        lank_filt = zeros(size(lank_noisy));
        lfeo_filt = zeros(size(lfeo_noisy));
        lfoo_filt = zeros(size(lfoo_noisy));
        lfop_filt = zeros(size(lfop_noisy));
        lfol_filt = zeros(size(lfol_noisy));
        for col = 1:3
            ltio_filt(:,col) = lpfilt(ltio_noisy(:,col), fc, fs);
            lank_filt(:,col) = lpfilt(lank_noisy(:,col), fc, fs);
            lfeo_filt(:,col) = lpfilt(lfeo_noisy(:,col), fc, fs);
            lfoo_filt(:,col) = lpfilt(lfoo_noisy(:,col), fc, fs);
            lfop_filt(:,col) = lpfilt(lfop_noisy(:,col), fc, fs);
            lfol_filt(:,col) = lpfilt(lfol_noisy(:,col), fc, fs);
        end
        
        LA_alpha_filt_deg = zeros(n,1);
        for k = 1:n
            [x_LT, y_LT, z_LT] = construct_LT_coords(ltio_filt(k,:), lank_filt(k,:), lfeo_filt(k,:));
            [x_LP, y_LP, z_LP] = construct_LP_coords(lfoo_filt(k,:), lfop_filt(k,:), lfol_filt(k,:));
            alpha_rad = compute_alpha_LA_strict(x_LP, y_LP, z_LT);
            LA_alpha_filt_deg(k) = rad2deg(alpha_rad);
        end
        
        figure('Name', sprintf('Filtered LA_\\alpha (Degrees) at fc = %.1f Hz', fc));
        plot(time, LA_alpha_filt_deg, 'b', 'LineWidth', 1.5);
        grid on;
        xlabel('Time (s)');
        ylabel('Angle (degrees)');
        title(sprintf('Filtered LA_\\alpha (Degrees) at fc = %.1f Hz', fc));
        
        MAE_values(iCutoff) = mean(abs(LA_alpha_clean_deg - LA_alpha_filt_deg));
        fprintf('fc = %.1f Hz, MAE = %.4f degrees\n', fc, MAE_values(iCutoff));
    end
    
    figure('Name','MAE vs. fc');
    plot(fc_values, MAE_values, 'r-o', 'LineWidth', 1.5, 'MarkerSize', 6);
    grid on;
    xlabel('Cutoff Frequency (Hz)');
    ylabel('Mean Absolute Error (degrees)');
    title('MAE vs. Cutoff Frequency for LA_\alpha');
    
    [bestMAE, bestIdx] = min(MAE_values);
    bestFc = fc_values(bestIdx);
    fprintf('\nBest cutoff frequency: fc = %.1f Hz, MAE = %.4f degrees\n', bestFc, bestMAE);
    
    
    %% Multiple Trials
  
    nTrials = 100;                               
    MAE_trials = zeros(length(fc_values), nTrials);  
    
    for trial = 1:nTrials
        ltio_noisy_trial = ltio + noise_std * randn(size(ltio));
        lank_noisy_trial = lank + noise_std * randn(size(lank));
        lfeo_noisy_trial = lfeo + noise_std * randn(size(lfeo));
        lfoo_noisy_trial = lfoo + noise_std * randn(size(lfoo));
        lfop_noisy_trial = lfop + noise_std * randn(size(lfop));
        lfol_noisy_trial = lfol + noise_std * randn(size(lfol));
        
        for iCutoff = 1:length(fc_values)
            fc = fc_values(iCutoff);
            
            ltio_filt_trial = zeros(size(ltio_noisy_trial));
            lank_filt_trial = zeros(size(lank_noisy_trial));
            lfeo_filt_trial = zeros(size(lfeo_noisy_trial));
            lfoo_filt_trial = zeros(size(lfoo_noisy_trial));
            lfop_filt_trial = zeros(size(lfop_noisy_trial));
            lfol_filt_trial = zeros(size(lfol_noisy_trial));
            for col = 1:3
                ltio_filt_trial(:,col) = lpfilt(ltio_noisy_trial(:,col), fc, fs);
                lank_filt_trial(:,col) = lpfilt(lank_noisy_trial(:,col), fc, fs);
                lfeo_filt_trial(:,col) = lpfilt(lfeo_noisy_trial(:,col), fc, fs);
                lfoo_filt_trial(:,col) = lpfilt(lfoo_noisy_trial(:,col), fc, fs);
                lfop_filt_trial(:,col) = lpfilt(lfop_noisy_trial(:,col), fc, fs);
                lfol_filt_trial(:,col) = lpfilt(lfol_noisy_trial(:,col), fc, fs);
            end
            
            LA_alpha_filt_trial_deg = zeros(n,1);
            for k = 1:n
                [x_LT, y_LT, z_LT] = construct_LT_coords(ltio_filt_trial(k,:), lank_filt_trial(k,:), lfeo_filt_trial(k,:));
                [x_LP, y_LP, z_LP] = construct_LP_coords(lfoo_filt_trial(k,:), lfop_filt_trial(k,:), lfol_filt_trial(k,:));
                alpha_rad = compute_alpha_LA_strict(x_LP, y_LP, z_LT);
                LA_alpha_filt_trial_deg(k) = rad2deg(alpha_rad);
            end
            
            MAE_trials(iCutoff, trial) = mean(abs(LA_alpha_clean_deg - LA_alpha_filt_trial_deg));
        end
    end
    
    avg_MAE = mean(MAE_trials, 2);
    
    fprintf('\n=== Step 5: Average MAE for each cutoff frequency after multiple trials ===\n');
    fprintf('Cutoff Frequency (Hz)    Average MAE (degrees)\n');
    fprintf('-----------------------------------------------\n');
    for iCutoff = 1:length(fc_values)
        fprintf('      %.1f                   %.4f\n', fc_values(iCutoff), avg_MAE(iCutoff));
    end
end


%%  other Function
 %(28)~(30)
function [x_LT, y_LT, z_LT] = construct_LT_coords(c_ltio, c_lank, c_lfeo)
    % (29)
    vec_y = c_lank - c_ltio;
    y_LT = vec_y / norm(vec_y);
    
    % (30) 
    vec_z = c_lfeo - c_ltio;
    z_LT = vec_z / norm(vec_z);
    
    % (28) 
    cross_x = cross(y_LT, z_LT);
    x_LT = cross_x / norm(cross_x);
end

% (34)~(36)
function [x_LP, y_LP, z_LP] = construct_LP_coords(c_lfoo, c_lfop, c_lfol)
    % (34)
    vec_x = c_lfop - c_lfoo;
    x_LP = vec_x / norm(vec_x);
    
    % (35) 
    vec_y = c_lfol - c_lfoo;
    y_LP = vec_y / norm(vec_y);
    
    % (36) 
    cross_z = cross(y_LP, x_LP);
    z_LP = cross_z / norm(cross_z);
end

%% (13) & (15)) ==========
function alpha_LA = compute_alpha_LA_strict(x_LP, y_LP, z_LT)
    % (15) 
    gamma_LA = asin(dot(x_LP, y_LP));
    
    cos_gamma = cos(gamma_LA);
    if abs(cos_gamma) < 1e-8
        cos_gamma = 1e-8;  % Prevent division by zero
    end
    
    % (13) 
    alpha_LA = -asin(dot(x_LP, z_LT) / cos_gamma);
end
