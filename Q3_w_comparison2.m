%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Numerical simulation of the evolution of a wavepacket in a 1D harmonic
%   trap using fast Fourier transform (FFT) with varying omega
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; clc; tic
close all;
set(0, 'DefaultFigureWindowStyle', 'docked');

% Define simulation parameters
a = -20;                       % Left end point
b = +20;                       % Right end point
L = b - a;                     % Width of the space
N = 512;                       % No. of spatial points
X = a + L * (0:N-1) / N;       % Dimensionless coordinates
P = (2 * pi / L) * [0:N/2-1, -N/2:-1]; % Dimensionless momentum

% DVR parameters
T = 200 * pi;                   % Total time duration
M = 10^5;                      % Total number of steps in the evolution
dt = T / M;                    % Time step duration
downsample_factor = 100;       % Factor to reduce the number of data points plotted

% Perturbation parameters
A = 0.01;                      % Fixed perturbation amplitude
w_values = [1.01, 1.05, 1.1]; % Omega values for comparison
w0 = 1.0;                      % Natural frequency
plot_theory_RWA = true; 

% Initialize variables for plots
fs = 22;                       % Font size for plots
transition_probabilities_all = zeros(length(w_values), M); % Preallocate for all omega

% Manually set custom colors for the curves
colors = [
    0.0, 0.5, 1.0; % Light blue
    0.0, 0.0, 0.8; % Dark blue
    0.8, 0.0, 0.0; % Dark red
    1.0, 0.4, 0.4  % Light red
];

% Initial state parameters
X0 = 0.0;                      % Wavepacket / Gaussian center
sigma = 1 / sqrt(w0);          % Width of the wavepacket / Gaussian

% Prepare ground and excited states for comparison
ground_temp = hermiteH(0, X) .* exp(-(X - X0).^2 / (2 * sigma^2));
ground = ground_temp / sqrt(ground_temp * ground_temp');

excited_temp = hermiteH(1, X) .* exp(-(X - X0).^2 / (2 * sigma^2));
excited = excited_temp / sqrt(excited_temp * excited_temp');

% Normalized initial state as ground state
Poly = hermiteH(0, X);
initial_state_temp = Poly .* exp(-(X - X0).^2 / (2 * sigma^2));
initial_state = initial_state_temp / sqrt(initial_state_temp * initial_state_temp');

% Iterate over omega values
for idx = 1:length(w_values)
    w = w_values(idx);         % Current omega
    psi_0 = initial_state;     % Reset initial state
    UT = exp(-1i * (P.^2 / 2) * dt); % Momentum space propagator
    fprintf('Simulating for w = %.2f\n', w);

    % Time evolution loop
    for m = 1:M
        UV = exp(-1i * ((X.^2) / 2 + A * cos(w * dt * (m - 1)) * sin(X)) * dt / 2);
        psi_1 = UV .* psi_0;
        phi_2 = fft(psi_1);
        phi_3 = UT .* phi_2;
        psi_3 = ifft(phi_3);
        psi_4 = UV .* psi_3;
        psi_0 = psi_4; % Update wavefunction

        % Calculate transition probability to the first excited state
        transition_probabilities_all(idx, m) = abs(dot(conj(excited), psi_0))^2;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting the transition probabilities for different omega
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



figure('Name', 'Transition Probabilities vs. Omega', 'Color', 'w');
hold on;

% Downsample the data for plotting
time = (0:M-1) * dt;
% Time evolution using the perturbation theory formula
matrix_element = A / (sqrt(2)*exp(1/4));
fprintf('Vnm/A = %.4f \n', matrix_element/A);
w_diff = w0 - 1.01;
transition_prob_theory_RWA = ( (matrix_element / w_diff) * sin( w_diff * time / 2) ) .^ 2; % remove factor of 4

time_downsampled = time(1:downsample_factor:end) / pi; % Downsample time
for idx = 1:length(w_values)
    probabilities_downsampled = transition_probabilities_all(idx, 1:downsample_factor:end); % Downsample probabilities
    plot(time_downsampled, probabilities_downsampled, 'LineWidth', 2.9, ...
         'Color', colors(idx, :), ...
         'DisplayName', sprintf('\\omega = %.1f', w_values(idx)));
end



    % Plot theoretical results
    if plot_theory_RWA
        theoretical_probabilities_RWA_downsampled = transition_prob_theory_RWA(1:downsample_factor:end);
        plot(time_downsampled, theoretical_probabilities_RWA_downsampled, 'LineWidth', 2.9, ...
             'Color', [1.0, 0.6, 0.2, 1.0], ... % Warm orange for TDPT + RWA
             'LineStyle', '-', ...
             'DisplayName', 'Theoretical (TDPT + RWA)');
    end
    


% Customize axes
xlim([min(time_downsampled), max(time_downsampled)]);

% Ensure valid ylim
max_probability = max(max(transition_probabilities_all));
if max_probability == 0
    max_probability = 1e-6; % Small positive value to avoid errors
end

%ylim([0, 0.08]);

xlabel('Time (multiples of $\pi$)', 'FontSize', fs, 'Interpreter', 'latex');
ylabel('Transition Probability $P_{1 \leftarrow 0}$', 'FontSize', fs, 'Interpreter', 'latex');

% Add grid, frame, and legend
grid on;
box on;
legend('FontSize', fs, 'Location', 'northeast', 'Box', 'on');
%xticks([0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20]);
%yticks([0, 0.02, 0.04, 0.06, 0.08]);

% Final adjustments
set(gca, 'FontSize', fs, 'LineWidth', 1.2);
pbaspect([2 1 1]);

% Save the plot
exportgraphics(gcf, 'Transition_Probabilities_Smoothed.png', 'Resolution', 300);

hold off;

fprintf('Simulation completed.\n');
toc;
