%% =========================================================================
%  EXAMPLE: run microEDM_nanoparticle_model.m for several parameter sets
%  and overlay the results on a single pair of plots.
%
%  Edit the `cases` list below to sweep whatever you like -- frequency,
%  voltage, current, duty cycle, or any combination. Each row is one run.
%  =========================================================================
clear; clc; close all;

% ---- Define the set of runs to compare -----------------------------------
% Columns: I_curr [A], V [V], freq [Hz], duty [fraction 0-1]
cases = [
    3,   40,  2400,  0.3;
    3,   40,  4000,  0.3;
    3,   40,  3000, 0.3;
    3,   40,  6000, 0.3;
];

n_cases = size(cases,1);
all_results = cell(n_cases,1);

% ---- Run the model for each case ------------------------------------------
for k = 1:n_cases
    I_curr = cases(k,1);
    V      = cases(k,2);
    freq   = cases(k,3);
    duty   = cases(k,4);

    all_results{k} = microEDM_nanoparticle_model_function(I_curr, V, freq, duty);
end

% ---- Overlay: average diameter vs time ------------------------------------
figure('Name','avg','visible','off');
hold on;
colors = lines(n_cases);
for k = 1:n_cases
    r = all_results{k};
    plot(r.t*1e6, r.d_avg, 'LineWidth', 1.8, 'Color', colors(k,:), ...
        'DisplayName', r.label);
end
hold off;
xlabel('Time (\mus)');
ylabel('Average particle diameter (nm)');
title('Average nanoparticle diameter vs time -- comparison');
legend('Location','best');
grid on;

% ---- Overlay: critical diameter vs time ------------------------------------
figure('Name','crit','visible','off');
hold on;
for k = 1:n_cases
    r = all_results{k};
    plot(r.t*1e9, r.d_star, 'LineWidth', 1.8, 'Color', colors(k,:), ...
        'DisplayName', r.label);
end
hold off;
xlabel('Time (ns)');
ylabel('Critical diameter, d^* (nm)');
title('Critical diameter vs time -- comparison');
legend('Location','best');
grid on;
saveas(figure(1), 'compare_avg.png');
saveas(figure(2), 'compare_crit.png');