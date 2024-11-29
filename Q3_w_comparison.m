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
T = 40 * pi;                   % Total time duration
M = 10^5;                      % Total number of steps in the evolution
dt = T / M;                    % Time step duration

% Perturbation parameters
A = 0.05;                      % Fixed perturbation amplitude
w_values = [1.05, 1.1, 1.2, 1.3, 2.0]; % 4 omega values for comparison
w0 = 1.0;                      % Natural frequency

% Initialize variables for plots
fs = 22;                       % Font size for plots
transition_probabilities_all = zeros(length(w_values), M); % Preallocate for all omega
cmap = [
    0.0, 0.447, 0.741; % Blue
    0.0, 0.75, 0.75;   % Cyan
    0.466, 0.674, 0.188; % Green
    0.85, 0.325, 0.098; % Orange
    0.635, 0.078, 0.184 % Dark red
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

% Downsample the data for smoother rendering
downsample_factor = 100;
time = (0:M-1) * dt;
time_downsampled = time(1:downsample_factor:end) / pi; % Downsampled time
for idx = 1:length(w_values)
    probabilities_downsampled = transition_probabilities_all(idx, 1:downsample_factor:end);
    plot(time_downsampled, probabilities_downsampled, 'LineWidth', 2.9, ...
         'Color', cmap(idx, :), ...
         'DisplayName', sprintf('\\omega = %.2f', w_values(idx)));
end

% Customize axes
xlim([0, T / pi]);
%ylim([0, 0.08]);

xlabel('Time (multiples of $\pi$)', 'FontSize', fs, 'Interpreter', 'latex');
ylabel('Transition Probability $P_{1 \leftarrow 0}$', 'FontSize', fs, 'Interpreter', 'latex');

% Add grid, frame, and legend
grid on;
box on;
legend('FontSize', fs, 'Location', 'northeast', 'Box', 'on');

% Final adjustments
set(gca, 'FontSize', fs, 'LineWidth', 1.2);
pbaspect([2 1 1]);

% Save the plot
exportgraphics(gcf, 'Transition_Probabilities_Fitted_Style.png', 'Resolution', 300);

hold off;

fprintf('Simulation completed.\n');
toc;
