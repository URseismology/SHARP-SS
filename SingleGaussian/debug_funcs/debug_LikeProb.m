% debug: calculate likelihood probability for different models

data = load(['/Users/evan/Documents/Research/THBD_test/' ...
    'H30_vs3p5_Amp0p2_ricker_wide_rayp100.mat']);
Z = data.Z;
R = data.R;
time = data.time;
time = time - 60;
Dt = time(2) - time(1);

prior = define_prior(time, Z, R);

loc_range = linspace(2, 40, 100);
amp_range = linspace(0.02, 0.5, 100);
wid_range = linspace(0.5, 10, 100);

% Define fixed values as variables
fixed_loc = 16.1;
fixed_amp = 0.2;
fixed_wid = 3.0;

% Preallocate arrays to store likelihood probabilities
LikeProb_loc_amp = zeros(length(loc_range), length(amp_range));
LikeProb_loc_wid = zeros(length(loc_range), length(wid_range));
LikeProb_amp_wid = zeros(length(amp_range), length(wid_range));

total_combinations = length(loc_range) * length(amp_range) + length(loc_range) * length(wid_range) + length(amp_range) * length(wid_range);
current_combination = 0;

% Loop over loc and amp
for i = 1:length(loc_range)
    for j = 1:length(amp_range)
        current_combination = current_combination + 1;
        % Use fixed values for amp and wid
        model = create_initial_model_customize(prior, loc_range(i), amp_range(j), fixed_wid);
        [LikeProb_loc_amp(i, j), ~, ~, ~, ~] = calculate_like_prob(Z, R, model, prior, 1);
        percent_finished = 100 * current_combination / total_combinations;
        fprintf('%.2f%% finished\n', percent_finished);
    end
end

% Loop over loc and wid
for i = 1:length(loc_range)
    for j = 1:length(wid_range)
        current_combination = current_combination + 1;
        % Use fixed values for amp
        model = create_initial_model_customize(prior, loc_range(i), fixed_amp, wid_range(j));
        [LikeProb_loc_wid(i, j), ~, ~, ~, ~] = calculate_like_prob(Z, R, model, prior, 1);
        percent_finished = 100 * current_combination / total_combinations;
        fprintf('%.2f%% finished\n', percent_finished);
    end
end

% Loop over amp and wid
for i = 1:length(amp_range)
    for j = 1:length(wid_range)
        current_combination = current_combination + 1;
        % Use fixed values for loc
        model = create_initial_model_customize(prior, fixed_loc, amp_range(i), wid_range(j));
        [LikeProb_amp_wid(i, j), ~, ~, ~, ~] = calculate_like_prob(Z, R, model, prior, 1);
        percent_finished = 100 * current_combination / total_combinations;
        fprintf('%.2f%% finished\n', percent_finished);
    end
end

%% plotting
% Plotting code

figure;

% Plot for loc and amp
subplot(1, 3, 1);
imagesc(amp_range, loc_range, LikeProb_loc_amp);
hold on;
[optimal_loc_amp_i, optimal_loc_amp_j] = find(LikeProb_loc_amp == max(LikeProb_loc_amp(:)));
plot(amp_range(optimal_loc_amp_j), loc_range(optimal_loc_amp_i), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
text(amp_range(optimal_loc_amp_j), loc_range(optimal_loc_amp_i) - 1, ['(' num2str(amp_range(optimal_loc_amp_j)) ', ' num2str(loc_range(optimal_loc_amp_i)) ')'], ...
    'Color', 'red', 'FontSize', 10, 'HorizontalAlignment', 'center');
title('Loc vs. Amp');
xlabel('Amp');
ylabel('Loc');
colorbar;

% Plot for loc and wid
subplot(1, 3, 2);
imagesc(wid_range, loc_range, LikeProb_loc_wid);
hold on;
[optimal_loc_wid_i, optimal_loc_wid_j] = find(LikeProb_loc_wid == max(LikeProb_loc_wid(:)));
plot(wid_range(optimal_loc_wid_j), loc_range(optimal_loc_wid_i), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
text(wid_range(optimal_loc_wid_j), loc_range(optimal_loc_wid_i) - 1, ['(' num2str(wid_range(optimal_loc_wid_j)) ', ' num2str(loc_range(optimal_loc_wid_i)) ')'], ...
    'Color', 'red', 'FontSize', 10, 'HorizontalAlignment', 'center');
title('Loc vs. Wid');
xlabel('Wid');
ylabel('Loc');
colorbar;

% Plot for amp and wid
subplot(1, 3, 3);
imagesc(wid_range, amp_range, LikeProb_amp_wid);
hold on;
[optimal_amp_wid_i, optimal_amp_wid_j] = find(LikeProb_amp_wid == max(LikeProb_amp_wid(:)));
plot(wid_range(optimal_amp_wid_j), amp_range(optimal_amp_wid_i), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
text(wid_range(optimal_amp_wid_j), amp_range(optimal_amp_wid_i) - 0.01, ['(' num2str(wid_range(optimal_amp_wid_j)) ', ' num2str(amp_range(optimal_amp_wid_i)) ')'], ...
    'Color', 'red', 'FontSize', 10, 'HorizontalAlignment', 'center');
title('Amp vs. Wid');
xlabel('Wid');
ylabel('Amp');
colorbar;
