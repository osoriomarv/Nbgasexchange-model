% Decay constants (lambda) for Hafnium and Tungsten
lambda_hafnium = log(2) / (9e6);   % Half-life of Hafnium is approximately 9 million years

% Time range for the decay plot
time = linspace(0, 1e8, 1e3);  % Change the range and number of points as desired

% Calculate the decay of Hafnium and Tungsten over time
hafnium = exp(-lambda_hafnium * time);

% Plot the decay curves
plot(time, hafnium, 'K-', 'LineWidth', 1.5);
% Set plot properties
xlabel('Time, Myrs');
ylabel('Mols left');
title('Decay of 182 Hafnium to 182 Tungsten');
legend('182 Hafnium');
