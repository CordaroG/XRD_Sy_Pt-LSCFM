% Clear workspace and close all figures
clear
close all

% Specify temperature and load corresponding dataset
scan_temp = '700';
dataset_name = append('Data Pt-LSCFM - ', scan_temp, 'C');
load(dataset_name);
load('zero_shift_mean.mat')                 % Vertical sample displacement calibrated through Pt experiments

% Minimum number of peaks for analysis
min_peaks_Pt = 9;
min_peaks_LSCFM = 2;

%% Constants and initialization
% Experimental setup parameters
energy = 12-0.046;                          % Adjusted beam energy [keV] calculated as Pt_a = Pt_a_cell_par0 @ XY=[0,0]
hc_constant = 1.23984198;                   % Planck constant times speed of light [keV·nm]
wavelength = 10*hc_constant/energy;         % Wavelength [Å] (exp = 1.037031129419311 Å)
omega_angle = 10;                           % Omega Angle for Asymmetric XRD Configuration [°]

% Reference Pt lattice parameter and fitting parameters
Pt_a_cell_par_exp = 3.92217600187318;       % Platinum lattice parameter
LSCFM_a_cell_par_exp = 3.864046561362875;   % LSCFM lattice parameter
zero_shift0 = 0;                            % Initial zero-shift value
min_peak_prominence_Pt = 2;                 % Min peak prominence for findpeaks function
min_peak_prominence_LSCFM = 0.6;            % Min peak prominence for findpeaks function

% Crystal planes (hkl) and squared sum
Pt_hkl_planes = [1 1 1; 2 0 0; 2 2 0; 3 1 1; 
                 2 2 2; 4 0 0; 3 3 1; 4 2 0; 4 2 2;
                 3 3 3; 4 4 0; 5 3 1; 6 0 0; 6 2 0; 5 3 3; 6 2 2];
Pt_hkl_squared_sum = sqrt(sum(Pt_hkl_planes(:,:).^2, 2));
LSCFM_hkl_planes = [2 0 0; 2 1 1; 2 2 0; 3 1 0];
LSCFM_hkl_squared_sum = sqrt(sum(LSCFM_hkl_planes(:,:).^2, 2));

%% Peak Identification and Fitting Initialization
% Define the 2Theta range and corresponding intensity of the first peak
x_range = x_theta(x_theta > 25 & x_theta < 28);
y_range = XRD_spectra(find(sample_coords(:, 1) == 0 & sample_coords(:, 2) == 0), x_theta > 25 & x_theta < 28);
% Identify initial peaks using findpeaks
[~, peak_position] = findpeaks(y_range, x_range, 'MinPeakProminence',10*min_peak_prominence_Pt, 'MinPeakDistance',0.5);
% Calculate experimental Pt lattice parameter from the first identified peak
Pt_a_cell_par0 = wavelength*Pt_hkl_squared_sum(1)/(2*sind((peak_position(1))/2));

% Define the 2Theta range and corresponding intensity of the first peak
x_range = x_theta(x_theta > 29 & x_theta < 32);
y_range = XRD_spectra(find(sample_coords(:, 1) == 6 & sample_coords(:, 2) == -6), x_theta > 29 & x_theta < 32);
% Identify initial peaks using findpeaks
[~, peak_position] = findpeaks(y_range, x_range, 'MinPeakProminence',5*min_peak_prominence_LSCFM, 'MinPeakDistance',0.5);
% Calculate experimental LSCFM lattice parameter from the first identified peak
LSCFM_a_cell_par0 = wavelength*LSCFM_hkl_squared_sum(1)/(2*sind((peak_position(1))/2));

% Initialize 2Theta and fitting parameters
Pt_peak_theta = 2*asind(wavelength./(2*Pt_a_cell_par0./Pt_hkl_squared_sum));
LSCFM_peak_theta = 2*asind(wavelength./(2*LSCFM_a_cell_par0./LSCFM_hkl_squared_sum));
tolerance_theta_Pt = Pt_peak_theta/180./cosd(Pt_peak_theta/2); % Delta2Theta range for find and fit peaks: function of 2Theta
tolerance_theta_Pt(1) = tolerance_theta_Pt(1) + 0.2;
tolerance_theta_LSCFM = LSCFM_peak_theta/180./cosd(LSCFM_peak_theta/2); % Delta2Theta range for find and fit peaks: function of 2Theta
tolerance_theta_LSCFM(1) = tolerance_theta_LSCFM(1) + 0.45;
fitting_width_Pt = 1./cosd(Pt_peak_theta/2);
fitting_width_LSCFM = 1./cosd(LSCFM_peak_theta/2);

% Preallocate arrays for results
num_spectra = length(sample_coords);
num_peak_Pt = length(Pt_peak_theta);
num_peak_LSCFM = length(LSCFM_peak_theta);
Pt_a_cell_par = zeros(num_spectra, 1);
LSCFM_a_cell_par = zeros(num_spectra, 1); 
zero_shift = zeros(num_spectra, 1);
fitted_peaks_index_Pt = zeros(num_spectra, 1);
fitted_peaks_index_LSCFM = zeros(num_spectra, 1);
fitted_peak_stand_deviat_Pt = zeros(num_spectra, 1);
fitted_peak_stand_deviat_LSCFM = zeros(num_spectra, 1);
cell_par_stand_deviat_Pt = zeros(num_spectra, 1);
cell_par_stand_deviat_LSCFM = zeros(num_spectra, 1);
fitted_peak_Pt = zeros(num_spectra, num_peak_Pt);
fitted_peak_LSCFM = zeros(num_spectra, num_peak_LSCFM);
missing_peak_Pt = zeros(num_spectra, num_peak_Pt);
missing_peak_LSCFM = zeros(num_spectra, num_peak_LSCFM);
err_fitted_peak_Pt = zeros(num_spectra, num_peak_Pt);
err_fitted_peak_LSCFM = zeros(num_spectra, num_peak_LSCFM);
err_unfitted_peak_Pt = zeros(num_spectra, num_peak_Pt, 2);
err_unfitted_peak_LSCFM = zeros(num_spectra, num_peak_LSCFM, 2);

% Options for fitting and minimization
options = optimset('MaxFunEvals', 1000, 'MaxIter', 1000, 'TolFun', 1E-12, 'TolX', 1E-12, 'display', 'none');

% Identify gaps in theta range for background adjustments
gap_theta_index = find(diff(x_theta) > 0.1);
gap_theta = [x_theta(gap_theta_index) x_theta(gap_theta_index + 1)];

%% Loop through all spectra for peak fitting and calculation of Pt lattice parameters
for i_spectra = 1 : num_spectra
    % Extract data, calculate background, and find peak parameters
    x_all = x_theta(x_theta > Pt_peak_theta(1)-3 & x_theta < 130);
    y_all = XRD_spectra(i_spectra, x_theta > Pt_peak_theta(1)-3 & x_theta < 130);
    [background, ~, ~] = backcor(x_all, y_all, 10, 0.00001, 'atq');
% % Uncomment the following line to open the GUI and set backcor parameters
% [EST, COEFS, IT] = backcor(x_all, y_all)
    intensity_data_wo_background = y_all - background';
    [peak_height, peak_position, peak_width] = findpeaks(intensity_data_wo_background, x_all, 'MinPeakProminence',min_peak_prominence_Pt, 'MinPeakDistance',0.5);
% % Uncomment the following line to plot and visually verify the fittings
% findpeaks(intensity_data_wo_background, x_all, 'MinPeakProminence',min_peak_prominence_Pt, 'MinPeakDistance',0.5); xline(Pt_peak_theta); hold on 

    key_fitted = zeros(num_peak_Pt, 1);
    fitted_peak_stand_deviat_Pt = zeros(num_peak_Pt, 1);

    [~, index] = max(peak_height);
    if peak_position(index) > 28.8 && peak_position(index) < 29.8
        % figure; findpeaks(intensity_data_wo_background, x_all, 'MinPeakProminence',min_peak_prominence_LSCFM, 'MinPeakDistance',0.5); xline(LSCFM_peak_theta); xline([LSCFM_peak_theta(1)-tolerance_theta_LSCFM(1); LSCFM_peak_theta(1)+tolerance_theta_LSCFM(1)], ':'); xline(Pt_peak_theta, 'b'); hold on
        continue
    end

    % Loop through all peaks to obtain positions and standard deviations
    for i_peak = 1 : num_peak_Pt
        j_peak = find(abs(Pt_peak_theta(i_peak) - peak_position) < tolerance_theta_Pt(i_peak));
        if isempty(j_peak)
            if i_peak == 1
                missing_peak_Pt(i_spectra, :) = 1;
                break
            end
            missing_peak_Pt(i_spectra, i_peak) = 1;
            continue
        elseif length(j_peak) > 1
            [~, jj] = min(abs(peak_position(j_peak) - Pt_peak_theta(i_peak)));
            j_peak = j_peak(jj);
        end

        % Select the range for peak fitting 
        y = intensity_data_wo_background(x_all > peak_position(j_peak) - fitting_width_Pt(i_peak)/2 & x_all < peak_position(j_peak) + fitting_width_Pt(i_peak)/2)';
        x = x_all(x_all > peak_position(j_peak) - fitting_width_Pt(i_peak)/2 & x_all < peak_position(j_peak) + fitting_width_Pt(i_peak)/2);

        % Obtain coefficients for a linear background
        bkg_points = [x(1), mean(mink(y(1:20), 3)); x(end), mean(mink(y(end-20:end), 3))];
        bkg_m = (bkg_points(2,2) - bkg_points(1,2))/((bkg_points(2,1)) - bkg_points(1,1));
        bkg_q = (bkg_points(2,1)*bkg_points(1,2) - bkg_points(1,1)*bkg_points(2,2))/(bkg_points(2,1) - bkg_points(1,1));

        % Fitting parameters as [Center Height FWHM Ratio]
        peak_width(j_peak) = min(peak_width(j_peak), tolerance_theta_Pt(i_peak));
        par0 = [peak_position(j_peak), peak_height(j_peak), peak_width(j_peak)*0.9 0.5];
        lower_bound = [peak_position(j_peak)-peak_width(j_peak)*2, peak_height(j_peak)*0.3, peak_width(j_peak)*0.3, 0];
        upper_bound = [peak_position(j_peak)+peak_width(j_peak)*2, peak_height(j_peak)*2.0, peak_width(j_peak)*2.0, 1];

        % Peak fitting with a Pseudo-Voigt function
        fun_pseudoVoigt = @(p,x) bkg_m*x + bkg_q + p(2)*((p(4)*exp(-((x-p(1))./p(3)*sqrt(4*log(2))).^2))+(1-p(4))./(1+((x-p(1))./p(3)*2).^2));
        [fit_result, fit_err, fit_residual, ~, ~, ~, jacobian] = lsqcurvefit(fun_pseudoVoigt,par0,x,y,lower_bound,upper_bound, options);
        fit_err_perc = 100*norm(fit_residual)/(sqrt(length(y))*max(y));
        r_squared = 1 - (fit_err)/(sum((y - mean(y)).^2));
        stand_deviat = sqrt(diag(inv(jacobian'*jacobian)*(fit_residual'*fit_residual)/(length(x)-length(par0))));
% % Uncomment the following line to plot and visually verify the fittings
% x_dense = linspace(x(1), x(end), 10*length(x)); y_fit = fun_pseudoVoigt(fit_result, x_dense); plot(x_dense,y_fit); scatter(x,y); xline(fit_result(1)); xlim([x(1)-1, x(end)+1]); text(0.2,0.92,[num2str(fit_result(1)), ' +- ', num2str(full(stand_deviat(1)))], 'units','normalized', 'FontSize',16, 'FontWeight','bold');
        % Check fitting results are not equal to boundary conditions
        if  any(abs(1 - fit_result(1:end-1)./upper_bound(1:end-1)) < 0.001) || any(abs(1 - fit_result(1:end-1)./lower_bound(1:end-1)) < 0.001)
            bound_check_Pt = [sample_coords(i_spectra, :), sqrt(sample_coords(i_spectra, 1)^2 + sample_coords(i_spectra, 2)^2), fit_result(1:end-1)./upper_bound(1:end-1), 99, fit_result(1:end-1)./lower_bound(1:end-1)];
        end
        % Check if same peak is fitted twice
        if i_peak > 1 && abs(fitted_peak_Pt(i_spectra, i_peak-1) - fit_result(1)) < 0.1
            fit_err_perc = 100;
        end
        % Check if errors are low
        if fit_err_perc < 10 && r_squared > 0.8 && all([fit_result(1) fit_result(1)-fit_result(3)/4 fit_result(1)+fit_result(3)/4] < gap_theta(:,1) | [fit_result(1) fit_result(1)-fit_result(3)/4 fit_result(1)+fit_result(3)/4] > gap_theta(:,2),'all')
            fitted_peak_Pt(i_spectra, i_peak) = fit_result(1);
            fitted_peak_stand_deviat_Pt(i_peak) = full(stand_deviat(1));
            key_fitted(i_peak) = 1;
            if i_peak == 1
                Pt_peak_theta(2:end) = 2*asind(wavelength./(2*wavelength*Pt_hkl_squared_sum(1)/(2*sind(fit_result(1)/2))./Pt_hkl_squared_sum(2:end)));
            end
        else
            err_unfitted_peak_Pt(i_spectra, i_peak, :) = [fit_err_perc r_squared];
        end
    end

    % Minimization function with peak positions weighted on their standard deviations as input,
    % while Pt lattice parameter, their standard deviation, and zero shift as outputs
    if length(find(key_fitted == 1)) >= min_peaks_Pt
        p1 = fitted_peak_Pt(i_spectra, key_fitted == 1)';
        p2 = Pt_hkl_squared_sum(key_fitted == 1);
        p3 = fitted_peak_stand_deviat_Pt(key_fitted == 1);
        fun_minim = @(X) (p1 + (X(2)*sind(p1)/sind(omega_angle)) - 2*asind(wavelength*p2/(2*X(1))))./p3;
        [minim_result, ~, minim_residual, ~, ~, ~, jacobian] = lsqnonlin(fun_minim, [Pt_a_cell_par0 zero_shift0], [], [], [], [], [], [], [], options);
        stand_deviat = sqrt(diag(inv(jacobian'*jacobian)*(minim_residual'*minim_residual)/(num_peak_Pt - length(minim_result))));

        fitted_peaks_index_Pt(i_spectra) = length(find(key_fitted == 1));
        Pt_a_cell_par(i_spectra) = minim_result(1);
        zero_shift(i_spectra) = minim_result(2);
        cell_par_stand_deviat_Pt(i_spectra) = full(stand_deviat(1));
        err_fitted_peak_Pt(i_spectra, :) = fitted_peak_Pt(i_spectra, :) - (2*asind(wavelength*Pt_hkl_squared_sum'./(2*Pt_a_cell_par(i_spectra))) - zero_shift(i_spectra)*sind(fitted_peak_Pt(i_spectra, :))/sind(omega_angle)).*key_fitted';
    else
        %% LSCFM
        if isnan(zero_shift_mean(i_spectra))
            missing_peak_LSCFM(i_spectra, :) = 2;
            continue
        end
        % Extract data, calculate background, and find peak parameters
        x_all = x_theta(x_theta > LSCFM_peak_theta(1)-5 & x_theta < 90);
        y_all = XRD_spectra(i_spectra, x_theta > LSCFM_peak_theta(1)-5 & x_theta < 90);
        [background, ~, ~] = backcor(x_all, y_all, 10, 0.00001, 'atq');
% % Uncomment the following line to open the GUI and set backcor parameters
% [EST, COEFS, IT] = backcor(x_all, y_all)
        intensity_data_wo_background = y_all - background';
        [peak_height, peak_position, peak_width] = findpeaks(intensity_data_wo_background, x_all, 'MinPeakProminence',min_peak_prominence_LSCFM, 'MinPeakDistance',0.5);
% % Uncomment the following line to plot and visually verify the fittings
% findpeaks(intensity_data_wo_background, x_all, 'MinPeakProminence',min_peak_prominence_LSCFM, 'MinPeakDistance',0.5); xline(LSCFM_peak_theta, ':'); hold on 
    
        key_fitted = zeros(num_peak_LSCFM, 1);
        fitted_peak_stand_deviat_LSCFM = zeros(num_peak_LSCFM, 1);
    
        % Loop through all peaks to obtain positions and standard deviations
        for i_peak = 1 : num_peak_LSCFM
            j_peak = find(abs(LSCFM_peak_theta(i_peak) - peak_position) < tolerance_theta_LSCFM(i_peak));
            if isempty(j_peak)
                if i_peak == 1
                    missing_peak_LSCFM(i_spectra, :) = 1;
                    break
                end
                missing_peak_LSCFM(i_spectra, i_peak) = 1;
                continue
            elseif length(j_peak) > 1
                [~, jj] = min(abs(peak_position(j_peak) - LSCFM_peak_theta(i_peak)));
                j_peak = j_peak(jj);
            end
    
            % Select the range for peak fitting 
            y = intensity_data_wo_background(x_all > peak_position(j_peak) - fitting_width_LSCFM(i_peak)/2 & x_all < peak_position(j_peak) + fitting_width_LSCFM(i_peak)/2)';
            x = x_all(x_all > peak_position(j_peak) - fitting_width_LSCFM(i_peak)/2 & x_all < peak_position(j_peak) + fitting_width_LSCFM(i_peak)/2);
    
            % Obtain coefficients for a linear background
            bkg_points = [x(1), mean(mink(y(1:20), 3)); x(end), mean(mink(y(end-20:end), 3))];
            bkg_m = (bkg_points(2,2) - bkg_points(1,2))/((bkg_points(2,1)) - bkg_points(1,1));
            bkg_q = (bkg_points(2,1)*bkg_points(1,2) - bkg_points(1,1)*bkg_points(2,2))/(bkg_points(2,1) - bkg_points(1,1));
    
            % Fitting parameters as [Center Height FWHM Ratio]
            peak_width(j_peak) = min(peak_width(j_peak), tolerance_theta_LSCFM(i_peak));
            par0 = [peak_position(j_peak), peak_height(j_peak), peak_width(j_peak)*0.9 0.5];
            lower_bound = [peak_position(j_peak)-peak_width(j_peak)*2, peak_height(j_peak)*0.3, peak_width(j_peak)*0.3, 0];
            upper_bound = [peak_position(j_peak)+peak_width(j_peak)*2, peak_height(j_peak)*2.0, peak_width(j_peak)*2.0, 1];
    
            % Peak fitting with a Pseudo-Voigt function
            fun_pseudoVoigt = @(p,x) bkg_m*x + bkg_q + p(2)*((p(4)*exp(-((x-p(1))./p(3)*sqrt(4*log(2))).^2))+(1-p(4))./(1+((x-p(1))./p(3)*2).^2));
            [fit_result, fit_err, fit_residual, ~, ~, ~, jacobian] = lsqcurvefit(fun_pseudoVoigt,par0,x,y,lower_bound,upper_bound, options);
            fit_err_perc = 100*norm(fit_residual)/(sqrt(length(y))*max(y));
            r_squared = 1 - (fit_err)/(sum((y - mean(y)).^2));
            stand_deviat = sqrt(diag(inv(jacobian'*jacobian)*(fit_residual'*fit_residual)/(length(x)-length(par0))));
% % Uncomment the following line to plot and visually verify the fittings
% x_dense = linspace(x(1), x(end), 10*length(x)); y_fit = fun_pseudoVoigt(fit_result, x_dense); plot(x_dense,y_fit); scatter(x,y); xline(fit_result(1)); xlim([x(1)-1, x(end)+1]); text(0.2,0.92,[num2str(fit_result(1)), ' +- ', num2str(full(stand_deviat(1)))], 'units','normalized', 'FontSize',16, 'FontWeight','bold');
            % Check fitting results are not equal to boundary conditions
            if  any(abs(1 - fit_result(1:end-1)./upper_bound(1:end-1)) < 0.001) || any(abs(1 - fit_result(1:end-1)./lower_bound(1:end-1)) < 0.001)
                bound_check_LSCFM = [sample_coords(i_spectra, :), sqrt(sample_coords(i_spectra, 1)^2 + sample_coords(i_spectra, 2)^2), fit_result(1:end-1)./upper_bound(1:end-1), 99, fit_result(1:end-1)./lower_bound(1:end-1)];
            end
            % Check if same peak is fitted twice
            if i_peak > 1 && abs(fitted_peak_LSCFM(i_spectra, i_peak-1) - fit_result(1)) < 0.1
                fit_err_perc = 100;
            end
            % Check if errors are low
            if fit_err_perc < 12 && r_squared > 0.7 && all([fit_result(1) fit_result(1)-fit_result(3)/4 fit_result(1)+fit_result(3)/4] < gap_theta(:,1) | [fit_result(1) fit_result(1)-fit_result(3)/4 fit_result(1)+fit_result(3)/4] > gap_theta(:,2),'all')
                fitted_peak_LSCFM(i_spectra, i_peak) = fit_result(1);
                fitted_peak_stand_deviat_LSCFM(i_peak) = full(stand_deviat(1));
                key_fitted(i_peak) = 1;
                if i_peak == 1
                    LSCFM_peak_theta(2:end) = 2*asind(wavelength./(2*wavelength*LSCFM_hkl_squared_sum(1)/(2*sind(fit_result(1)/2))./LSCFM_hkl_squared_sum(2:end)));
                end
            else
                err_unfitted_peak_LSCFM(i_spectra, i_peak, :) = [fit_err_perc r_squared];
            end
        end
    
        % Minimization function with peak positions weighted on their standard deviations as input,
        % while LSCFM lattice parameter, their standard deviation, and zero shift as outputs
        if length(find(key_fitted == 1)) >= min_peaks_LSCFM
            p1 = fitted_peak_LSCFM(i_spectra, key_fitted == 1)';
            p2 = LSCFM_hkl_squared_sum(key_fitted == 1);
            p3 = fitted_peak_stand_deviat_LSCFM(key_fitted == 1);
            fun_minim = @(X) (p1 + (zero_shift_mean(i_spectra)*sind(p1)/sind(omega_angle)) - 2*asind(wavelength*p2/(2*X(1))))./p3;
            [minim_result, ~, minim_residual, ~, ~, ~, jacobian] = lsqnonlin(fun_minim, LSCFM_a_cell_par0, [], [], [], [], [], [], [], options);
            stand_deviat = sqrt(diag(inv(jacobian'*jacobian)*(minim_residual'*minim_residual)/(num_peak_LSCFM - length(minim_result))));
    
            fitted_peaks_index_LSCFM(i_spectra) = length(find(key_fitted == 1));
            LSCFM_a_cell_par(i_spectra) = minim_result;
            cell_par_stand_deviat_LSCFM(i_spectra) = full(stand_deviat(1));
            err_fitted_peak_LSCFM(i_spectra, :) = fitted_peak_LSCFM(i_spectra, :) - (2*asind(wavelength*LSCFM_hkl_squared_sum'./(2*LSCFM_a_cell_par(i_spectra))) - zero_shift_mean(i_spectra)*sind(fitted_peak_LSCFM(i_spectra, :))/sind(omega_angle)).*key_fitted';
        end
    end
end

%% Results post-treatments to verify and remove outliers
goniom_radius = 645; % mm
sample_displ = zero_shift*pi*goniom_radius/180;

Coeff_a_Pt_into_T = [173404826.220606 -670583934.327650 -2787.51427855167];
fun_a_Pt_into_T = @(a) sqrt(Coeff_a_Pt_into_T(1)*a + Coeff_a_Pt_into_T(2)) + Coeff_a_Pt_into_T(3) - 273.15;

err_sum_table = zeros(max(fitted_peaks_index_Pt), 5);
for iii = min_peaks_Pt : max(fitted_peaks_index_Pt)
    err_sum = 0; count = 0;
    for i_spectra = 1 : length(sample_coords)
        if fitted_peaks_index_Pt(i_spectra) == iii
            err_sum = err_sum + sum(abs(err_fitted_peak_Pt(i_spectra, :)));
            count = count + 1;
        end
    end
    err_sum_table(iii+1, :) = [iii; err_sum; count; err_sum/count; err_sum/count/iii];
end

disp(err_sum_table(min_peaks_Pt+1 : end, :))
min_peaks = input('Select the min number of peaks fitted :   ');

sample_displ_plot = sample_displ;
Pt_a_plot = Pt_a_cell_par;
Pt_stand_deviat_plot = cell_par_stand_deviat_Pt;
LSCFM_a_plot = LSCFM_a_cell_par;
LSCFM_stand_deviat_plot = cell_par_stand_deviat_LSCFM;

sample_displ_plot(fitted_peaks_index_Pt < min_peaks) = NaN;
Pt_a_plot(fitted_peaks_index_Pt < min_peaks) = NaN;
Pt_stand_deviat_plot(fitted_peaks_index_Pt < min_peaks) = NaN;
LSCFM_a_plot(fitted_peaks_index_LSCFM < min_peaks_LSCFM) = NaN;
LSCFM_stand_deviat_plot(fitted_peaks_index_LSCFM < min_peaks_LSCFM) = NaN;

figure; scatter3(sample_coords(:, 1), sample_coords(:, 2), Pt_a_plot, 100, Pt_a_plot, 'square', 'filled'); colormap jet; set(gcf, 'Position', [100, 150, 861, 737]); view(-10, 2)
a_range = [Pt_a_plot - Pt_stand_deviat_plot, Pt_a_plot + Pt_stand_deviat_plot];
t_range = fun_a_Pt_into_T(a_range);
t_diff = (t_range(:, 2) - t_range(:, 1))/2;
figure; scatter3(sample_coords(:, 1), sample_coords(:, 2), t_diff, 100, t_diff, 'square', 'filled'); box on; grid on; colormap jet; c = colorbar; c.Label.String = 'Temperature [°C]'; c.FontSize = 14; clim([2 22]); set(gca, 'fontsize', 12); xlabel('X [mm]'); ylabel('Y [mm]'); zlabel('Temperature Range [°C]'); xlim([-50,50]); ylim([-50,50]); zlim([2 22]); xticks(-50:10:50); yticks(-50:10:50); set(gcf, 'Position', [950, 150, 861, 737]); view(-10, 1)

z_limit0  = input('Select the min & max of aPt values to cut outliers [min,max] :   ');
z_limit_max_t  = input('Select the max T values to cut outliers :   ');
figure; scatter3(sample_coords(:, 1), sample_coords(:, 2), LSCFM_a_plot, 100, LSCFM_a_plot, 's', 'filled'); box on; colormap jet; set(gca, 'fontsize', 20, 'fontname', 'calibri'); c = colorbar; c.Label.String = 'LSCFM Lattice Parameter [Å]'; c.FontSize = 20; c.Label.FontSize = 22; c.Label.FontWeight = 'bold'; xlabel('Horizontal Coordinate X [mm]', 'FontWeight', 'bold'); ylabel('Vertical Coordinate Y [mm]', 'FontWeight', 'bold'); set(gcf, 'Position', [100, 200, 843, 737]); view(-10, 1); xlim([-50,50]); ylim([-50,50]); xticks(-50:10:50); yticks(-50:10:50)
figure; scatter3(sample_coords(:, 1), sample_coords(:, 2), LSCFM_stand_deviat_plot, 100, LSCFM_stand_deviat_plot, 's', 'filled'); box on; colormap jet; set(gca, 'fontsize', 20, 'fontname', 'calibri'); c = colorbar; c.Label.String = 'LSCFM Standard Deviation [Å]'; c.FontSize = 20; c.Label.FontSize = 22; c.Label.FontWeight = 'bold'; xlabel('Horizontal Coordinate X [mm]', 'FontWeight', 'bold'); ylabel('Vertical Coordinate Y [mm]', 'FontWeight', 'bold'); set(gcf, 'Position', [950, 150, 861, 737]); view(-10, 1); xlim([-50,50]); ylim([-50,50]); xticks(-50:10:50); yticks(-50:10:50)

z_limit0_LSCFM  = input('Select the min & max of aLSCFM values to cut outliers [min,max] :   ');
z_limit_max_SD_LSCFM  = input('Select the max SD values to cut outliers :   ');

user_input = [min_peaks z_limit_max_t; z_limit0; min_peaks_LSCFM z_limit_max_SD_LSCFM; z_limit0_LSCFM];

Pt_a_plot(t_diff > z_limit_max_t) = NaN;
a_range(t_diff > z_limit_max_t, :) = NaN;
t_range(t_diff > z_limit_max_t, :) = NaN;
a_stand_deviat = Pt_stand_deviat_plot(t_diff < z_limit_max_t);
t_diff(t_diff > z_limit_max_t) = NaN;
a_limit = [min(Pt_a_plot(Pt_a_plot > z_limit0(1))) max(Pt_a_plot(Pt_a_plot < z_limit0(2)))];
t_limit = fun_a_Pt_into_T(a_limit);
t_plot = fun_a_Pt_into_T(Pt_a_plot);

LSCFM_a_plot(LSCFM_stand_deviat_plot > z_limit_max_SD_LSCFM) = NaN;
LSCFM_stand_deviat_plot(LSCFM_stand_deviat_plot > z_limit_max_SD_LSCFM) = NaN;
LSCFM_stand_deviat_plot(LSCFM_a_plot < z_limit0_LSCFM(1) | LSCFM_a_plot > z_limit0_LSCFM(2)) = NaN;
LSCFM_a_plot(LSCFM_a_plot < z_limit0_LSCFM(1) | LSCFM_a_plot > z_limit0_LSCFM(2)) = NaN;

% Results as [min, mean, max, diff, SD, N° Pt_a] within limits for Pt_a & T
result_table = [a_limit(1) mean(Pt_a_plot, 'omitnan') a_limit(2) diff(a_limit) mean(a_stand_deviat) length(find(Pt_a_plot > 0))
                t_limit(1) fun_a_Pt_into_T(mean(Pt_a_plot, 'omitnan')) t_limit(2) diff(t_limit) mean(t_diff, 'omitnan') length(find(~isnan(Pt_a_plot)))];
disp(result_table)
result_table2 = [sample_coords, Pt_a_plot, t_plot, sample_displ_plot, Pt_stand_deviat_plot];

%% Figures
figure; scatter(sample_coords(:, 1), sample_coords(:, 2), 300, LSCFM_a_plot, 's', 'filled'); box on; colormap jet; set(gca, 'fontsize', 20, 'fontname', 'calibri'); c = colorbar; c.Label.String = 'LSCFM Lattice Parameter [Å]'; c.FontSize = 20; c.Label.FontSize = 22; c.Label.FontWeight = 'bold'; xlabel('Horizontal Coordinate X [mm]', 'FontWeight', 'bold'); ylabel('Vertical Coordinate Y [mm]', 'FontWeight', 'bold'); set(gcf, 'Position', [100, 100, 840, 730]); xlim([-50,50]); ylim([-50,50]); xticks(-50:10:50); yticks(-50:10:50)

figure; scatter(sample_coords(:, 1), sample_coords(:, 2), 300, LSCFM_stand_deviat_plot, 's', 'filled'); box on; colormap jet; set(gca, 'fontsize', 20, 'fontname', 'calibri'); c = colorbar; c.Label.String = 'LSCFM Standard Deviation [Å]'; c.FontSize = 20; c.Label.FontSize = 22; c.Label.FontWeight = 'bold'; xlabel('Horizontal Coordinate X [mm]', 'FontWeight', 'bold'); ylabel('Vertical Coordinate Y [mm]', 'FontWeight', 'bold'); set(gcf, 'Position', [200, 100, 840, 730]); xlim([-50,50]); ylim([-50,50]); xticks(-50:10:50); yticks(-50:10:50)

figure; scatter(sample_coords(:, 1), sample_coords(:, 2), 300, Pt_a_plot, 's', 'filled'); box on; colormap jet; set(gca, 'fontsize', 20, 'fontname', 'calibri'); c = colorbar; c.Label.String = 'Pt Lattice Parameter [Å]'; c.FontSize = 20; c.Label.FontSize = 22; c.Label.FontWeight = 'bold'; xlabel('Horizontal Coordinate X [mm]', 'FontWeight', 'bold'); ylabel('Vertical Coordinate Y [mm]', 'FontWeight', 'bold'); clim(a_limit); set(gcf, 'Position', [300, 100, 840, 730]); xlim([-50,50]); ylim([-50,50]); xticks(-50:10:50); yticks(-50:10:50)

figure; scatter(sample_coords(:, 1), sample_coords(:, 2), 300, t_plot(t_plot ~= 0), 's', 'filled'); box on; colormap jet; set(gca, 'fontsize', 20, 'fontname', 'calibri'); c = colorbar; c.Label.String = 'Temperature [°C]'; c.FontSize = 20; c.Label.FontSize = 22; c.Label.FontWeight = 'bold'; xlabel('Horizontal Coordinate X [mm]', 'FontWeight', 'bold'); ylabel('Vertical Coordinate Y [mm]', 'FontWeight', 'bold'); clim(fun_a_Pt_into_T(a_limit)); set(gcf, 'Position', [400, 100, 840, 730]); xlim([-50,50]); ylim([-50,50]); xticks(-50:10:50); yticks(-50:10:50)

figure; scatter3(sample_coords(:, 1), sample_coords(:, 2), t_diff, 300, t_diff, 'square', 'filled'); box on; grid on; colormap jet; set(gca, 'fontsize', 20, 'fontname', 'calibri'); c = colorbar; c.Label.String = 'Temperature [°C]'; c.FontSize = 20; c.Label.FontSize = 22; c.Label.FontWeight = 'bold'; xlabel('X [mm]', 'FontWeight', 'bold'); ylabel('Y [mm]', 'FontWeight', 'bold'); zlabel('Temperature Range [°C]', 'FontWeight', 'bold'); set(gcf, 'Position', [500, 100, 860, 730]); xlim([-50,50]); ylim([-50,50]); xticks(-50:10:50); yticks(-50:10:50)