%% COVID-19 Contact Tracing Simulation
% High-risk vs Low-risk Group Model with Contact Tracing
% Modified Gillespie Algorithm

clc; clear;

%% Simulation Settings
num_of_run = 100;       % Number of simulation runs (use 10000 for full analysis)

disp('Starting COVID-19 Contact Tracing Simulation...');
tic;

%% Load Parameters
input = Parameters();

%% Run Simulations
Results = cell(num_of_run, 1);
TotalResults = cell(num_of_run, 1);

parfor i = 1:num_of_run
    if mod(i, 50) == 0
        disp(['Running simulation ', num2str(i), ' of ', num2str(num_of_run)]);
    end
    
    % Run single simulation
    results = ModiGillespie_algorithm(input);
    TotalResults{i} = results;
    
    % Process results for analysis
    interval = 0.1;
    max_time = 34;
    TotalResults{i}.tspan(isinf(TotalResults{i}.tspan)) = max_time;
    TotalResults{i}.tspan(isnan(TotalResults{i}.tspan)) = max_time;
    
    Results{i}.num_reactionlist = zeros(41, ceil(max_time/interval));
    
    for ii = 1:max_time/interval
        t_start = (ii - 1) * interval;
        t_end = ii * interval;
        idxdice = (TotalResults{i}.tspan >= t_start) & (TotalResults{i}.tspan < t_end);
        
        if ii == max_time/interval
            idxdice = (TotalResults{i}.tspan >= t_start) & (TotalResults{i}.tspan <= t_end);
        end
        
        for iii = 1:41
            Results{i}.num_reactionlist(iii, ii) = nnz(TotalResults{i}.reactionlist(idxdice) == iii);
        end
    end
end

simulation_time = toc;
disp(['Simulation completed in ', num2str(simulation_time), ' seconds']);

%% Analyze Results
num_intervals = size(Results{1}.num_reactionlist, 2); 
dayarray = 0.1:0.1:(num_intervals*0.1);  

% Calculate cumulative cases for high-risk group
Hcase = zeros(num_of_run, num_intervals);
for i = 1:num_of_run
    Hmake_cum = zeros(1, size(Results{i}.num_reactionlist, 2));
    for ii = [7 9 14 20 21]  % Confirmed case reactions for high-risk
        Hmake_cum = Hmake_cum + Results{i}.num_reactionlist(ii, :);
    end
    Hcase(i, :) = cumsum(Hmake_cum, 2);
end

% Calculate cumulative cases for low-risk group
Lcase = zeros(num_of_run, num_intervals);
for i = 1:num_of_run
    Lmake_cum = zeros(1, size(Results{i}.num_reactionlist, 2));
    for ii = [8 10 15 22 23]  % Confirmed case reactions for low-risk
        Lmake_cum = Lmake_cum + Results{i}.num_reactionlist(ii, :);
    end
    Lcase(i, :) = cumsum(Lmake_cum, 2);
end

%% Display Summary Statistics
disp('=== Simulation Results Summary ===');
disp(['High-risk group final cases: Mean = ', num2str(mean(Hcase(:, end))), ...
      ', Median = ', num2str(median(Hcase(:, end))), ...
      ', 95% CI = [', num2str(prctile(Hcase(:, end), 2.5)), ', ', ...
      num2str(prctile(Hcase(:, end), 97.5)), ']']);

disp(['Low-risk group final cases: Mean = ', num2str(mean(Lcase(:, end))), ...
      ', Median = ', num2str(median(Lcase(:, end))), ...
      ', 95% CI = [', num2str(prctile(Lcase(:, end), 2.5)), ', ', ...
      num2str(prctile(Lcase(:, end), 97.5)), ']']);

%% Create Simple Histogram Plot
figure('Units', 'centimeters', 'Position', [2, 2, 30, 15]);
t = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% High-risk group histogram
nexttile(1);
h1 = histogram(Hcase(:, end), 'Normalization', 'probability', ...
              'FaceColor', [0.8 0.3 0.3], 'EdgeColor', 'white', ...
              'LineWidth', 0.5, 'NumBins', 30);
xlabel('Cumulative Cases', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Probability', 'FontSize', 12, 'FontWeight', 'bold');
title('High-risk Group Final Case Distribution', 'FontSize', 14, 'FontWeight', 'bold');
grid on;
set(gca, 'GridAlpha', 0.3, 'FontSize', 11);

% Add statistical lines
mean_h = mean(Hcase(:, end));
median_h = median(Hcase(:, end));
ylims = ylim;
line([mean_h mean_h], ylims, 'Color', [0.1 0.1 0.1], 'LineWidth', 2, 'LineStyle', '--');
line([median_h median_h], ylims, 'Color', [0.4 0.4 0.4], 'LineWidth', 2, 'LineStyle', ':');
legend('Distribution', 'Mean', 'Median', 'Location', 'northeast', 'FontSize', 10);

% Low-risk group histogram  
nexttile(2);
h2 = histogram(Lcase(:, end), 'Normalization', 'probability', ...
              'FaceColor', [0 0.4470 0.7410], 'EdgeColor', 'white', ...
              'LineWidth', 0.5, 'NumBins', 30);
xlabel('Cumulative Cases', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Probability', 'FontSize', 12, 'FontWeight', 'bold');
title('Low-risk Group Final Case Distribution', 'FontSize', 14, 'FontWeight', 'bold');
grid on;
set(gca, 'GridAlpha', 0.3, 'FontSize', 11);

% Add statistical lines
mean_l = mean(Lcase(:, end));
median_l = median(Lcase(:, end));
ylims = ylim;
line([mean_l mean_l], ylims, 'Color', [0.1 0.1 0.1], 'LineWidth', 2, 'LineStyle', '--');
line([median_l median_l], ylims, 'Color', [0.4 0.4 0.4], 'LineWidth', 2, 'LineStyle', ':');
legend('Distribution', 'Mean', 'Median', 'Location', 'northeast', 'FontSize', 10);

% Overall title
title(t, 'COVID-19 Simulation: Final Case Distribution Comparison', 'FontSize', 16, 'FontWeight', 'bold');

%% Save Results
save(['Results_', datestr(now, 'yyyymmdd'), '.mat'], 'Results', 'Hcase', 'Lcase', ...
      'input', 'simulation_time', 'dayarray', '-v7.3');

disp('Results saved successfully!');
disp('=== Simulation Complete ===');