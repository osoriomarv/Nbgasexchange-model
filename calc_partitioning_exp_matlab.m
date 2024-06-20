%

% Some important constants/variables

% Arrays of percentiles from which to select error bounds
lower_bnd_sigma_percentiles = [50.0, 15.9, 2.3, 0.15];
upper_bnd_sigma_percentiles = 100.0 - lower_bnd_sigma_percentiles;

% PARAMS
% The names and folder references for the different fits presented in Chidester & Lock
fit_names = {'no_norm_unweighted_all', 'Mg_diss_no_norm_unweighted_all', 'no_norm_unweighted_all_metal', 'no_norm_unweighted_abc_only'};
fit_names_print = {'Full fit with Mg as exchange', 'Full fit with Mg as dissociation', 'Fit with only metal terms', 'Fit with no compositional terms'};

% Select the fit you want to use
ind_fit = 1;
fprintf('Using %s\n', fit_names_print{ind_fit});

% Find the binary fit files associated with the selected fit
binary_fit_file = strcat(fit_names{ind_fit}, '/exLS_fits_ftest_', fit_names{ind_fit}, '.bin');

% What sigma errors to report for errors
sigma = 1;

% Number of MC points for calculating error
NMC_calc = 5000;

% Pressure and temperature at which to calculate the partitioning
pcalc = 10;
Tcalc = 2000;

% Which elements do you want to find the partitioning for ('O' and 'Fe' must be the first two terms)
elem_part_calc = {'O', 'Fe', 'Si', 'Mg', 'S'};

% The composition of the metal and silicate you want to calculate KD for (mol%).
% List in the same order as elements are specified
% O for silicate is calculated self-consistently for charge balance
met_calc_user = [0.02, 0.92, 0.05, 0.01, 0];
sil_calc_user = [0, 0.02, 0.5, 0.48, 0];

% Other examples
% elem_part_calc = {'O', 'Fe', 'Al', 'Ca', 'Mg', 'Si', 'U'};
% met_calc_user = [0.082853, 0.63738565, 0, 0, 0, 0.27976135, 0];
% sil_calc_user = [0, 0.09494723, 0.01505161, 0.03854519, 0.2065641, 0.62451705, 0.02037482];

% elem_part_calc = {'O', 'Fe', 'Al', 'Ca', 'K', 'Mg', 'S', 'Si'};
% met_calc_user = [1.19119320e-01, 4.90837373e-01, 1.11550225e-04, 1.97405083e-04, ...
%                  7.69821521e-05, 1.84339601e-03, 8.27666345e-02, 3.38801874e-03];
% sil_calc_user = [0, 0.15510004, 0.0431348, 0.04117791, 0.00771705, 0.38461031, 0.0046111, 0.3636488];


% Read info from the binary file for the fit
[filename, Nelem, elem, elem_met, elem_sil, valence, ind_elem_disc, ...
    Nelem_part, flag_elem_part, ind_elem_part, flag_term, ind_term_fit, ...
    flag_LS, flag_MC, flag_norm_Nsamp_elem, flag_calc_weights, Ntest, NMC, sigma_omit, ...
    lngammai0_ref, term_default, flag_term_default, ...
    Nsamp, study_names, sample_names, ...
    p, T, met, sil, logK_part, logK_part_err, ...
    ind_active_elem, ind_active_elem_calc, ...
    fits, cov, corr] = read_fit_binary(binary_fit_file);


% Find the indices that relate the elements asked for to those in the fit file
ind_flag_calc_user = zeros(size(met_calc_user));
for i = 1:numel(elem_part_calc)
    temp = find(strcmp(elem, elem_part_calc{i}));
    if isempty(temp)
        error(['Unknown element ', elem_part_calc{i}]);
    end
    ind_flag_calc_user(i) = temp(1);
end

% Define the general composition and flag arrays for all elements
met_calc = zeros(1, Nelem);
sil_calc = zeros(1, Nelem);

met_calc(1, ind_flag_calc_user) = met_calc_user;
sil_calc(1, ind_flag_calc_user) = sil_calc_user;

flag_elem_part_calc = zeros(1, Nelem);
flag_elem_part_calc(ind_flag_calc_user) = 1;


% Calculate logKD using the best fit coefficients
[temp, ind_elem_store, ind_samp_store, logKcalc] = partfit(fits, ...
    pcalc, Tcalc, met_calc, sil_calc, ones(1, Nelem), ones(1, Nelem), ...
    flag_elem_part_calc, ind_flag_calc_user, flag_elem_part_calc, ...
    flag_term, ind_term_fit, Nelem, valence, lngammai0_ref, term_default, ...
    flag_term_default, ind_active_elem_calc, ind_elem_disc, ...
    'flag_output_mode', 1, 'flag_norm_Nsamp_elem', flag_norm_Nsamp_elem);

% Record the relevant values
res = temp{1};
ind_elem_store = temp{2};
ind_samp_store = temp{3};
logKcalc = temp{4};

% Print the 'best fit' values
fprintf('Values calculated with best-fit values\n');
fprintf('Element \t logKD \t KD\n');
for i = 1:numel(ind_elem_store)
    if i ~= 2
        fprintf('%s \t %.6f \t %.6e\n', elem{ind_elem_store(i)}, logKcalc(i), 10^logKcalc(i));
    end
end

% Do a random sampling of the parameters using the covariance matrix to determine
% median and distribution more 'correctly'
param_test = mvnrnd(fits, cov, NMC_calc);
logKcalc_test = zeros(NMC_calc, numel(ind_samp_store));
Kcalc_test = zeros(NMC_calc, numel(ind_samp_store));

% Sometimes the multivariate normal distribution function throws an error, so we need to check
temp = isPSD(cov, corr, tol=1e-8);
if temp(1) && temp(2)
    fprintf('Ignore positive-definite error from mvnrnd, within tolerance\n');
else
    error('Covariance matrix is not positive-definite to within tolerance');
end

for i = 1:NMC_calc
    temp = partfit(param_test(i, :), pcalc, Tcalc, met_calc, sil_calc, ...
        ones(1, Nelem), ones(1, Nelem), flag_elem_part_calc, ind_flag_calc_user, ...
        flag_elem_part_calc, flag_term, ind_term_fit, Nelem, valence, ...
        lngammai0_ref, term_default, flag_term_default, ind_active_elem_calc, ...
        ind_elem_disc, 'flag_output_mode', 1, 'flag_norm_Nsamp_elem', flag_norm_Nsamp_elem);

    logKcalc_test(i, :) = temp{4};
    Kcalc_test(i, :) = 10 .^ temp{4};
end

fprintf('done\n');


% Find the percentiles of logKD
lower_log = prctile(logKcalc_test, lower_bnd_sigma_percentiles(sigma), 1);
upper_log = prctile(logKcalc_test, upper_bnd_sigma_percentiles(sigma), 1);
median_log = prctile(logKcalc_test, 50, 1);

% Report the error
fprintf('Distribution of values calculated with Monte Carlo\n');
for i = 1:numel(ind_elem_store)
    if i ~= 2
        fprintf('logKD_%s = %.6f + %.6f - %.6f\n', elem{ind_elem_store(i)}, ...
            median_log(i), upper_log(i) - median_log(i), median_log(i) - lower_log(i));
    end
end

% Same for KD
fprintf('\n');
lower = prctile(Kcalc_test, lower_bnd_sigma_percentiles(sigma), 1);
upper = prctile(Kcalc_test, upper_bnd_sigma_percentiles(sigma), 1);
median = prctile(Kcalc_test, 50, 1);

% Report the error
for i = 1:numel(ind_elem_store)
    if i ~= 2
        fprintf('KD_%s = %.6e + %.6e - %.6e\n', elem{ind_elem_store(i)}, ...
            median(i), upper(i) - median(i), median(i) - lower(i));
    end
end


% Plot the result
fig = figure;
fig.Units = 'inches';
fig.Position = [0, 0, 3.5, 3 * numel(ind_elem_store)];
gs = gridspec.GridSpec(numel(ind_elem_store), 1);

% Make an array of axes
ax = gobjects(numel(ind_elem_store), 1);

for k = 1:numel(ind_elem_store)
    ax(k) = subplot(gs(k));
end

font = struct('family', 'helvetica', 'weight', 'normal', 'size', 8);
mpl.rcParams('font', font);

for i = 1:numel(ind_elem_store)
    hist = histogram(ax(i), Kcalc_test(:, i), 'BinWidth', 100);
    hold(ax(i), 'on');
    
    plot(ax(i), [median(i), median(i)], [0, max(hist.Values)], 'k-');
    plot(ax(i), [10^logKcalc(i), 10^logKcalc(i)], [0, max(hist.Values)], 'm:');
    plot(ax(i), [lower(i), lower(i)], [0, max(hist.Values)], 'r-');
    plot(ax(i), [upper(i), upper(i)], [0, max(hist.Values)], 'r-');
    
    hold(ax(i), 'off');
end


