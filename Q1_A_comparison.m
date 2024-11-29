% Clear workspace and set up
clear all; clc; tic;
close all;
set(0, 'DefaultFigureWindowStyle', 'docked');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Numerical simulation comparing different perturbation amplitudes (A)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define simulation parameters
a = -20;                      % Left end point
b = +20;                      % Right end point
L = b - a;                    % Width of the space
N = 512;                      % No. of spatial points
X = a + L * (0:N-1) / N;      % Dimensionless coordinates
P = (2 * pi / L) * [0:N/2-1, -N/2:-1]; % Dimensionless momentum

% Time and frequency parameters
T = 50*pi;                      % Total time duration
M = 10^5;                     % Fixed number of time steps
dt = T / M;                   % Time step duration
A_values = [0.01, 0.1, 1, 100]; % Array of perturbation amplitudes
w = 1.1;                      % Perturbation frequency
w0 = 1.0;                     % Harmonic oscillator frequency

% Plot parameters
plot_initial_state = true;
plot_animation = true;
plot_final_state = true;
frames = 100;
fs = 22;
x_lim = 20;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Define the initial state psi_0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X0 = 0.0;                     % Wavepacket center
sigma = 1 / sqrt(w0);         % Width of the wavepacket

% Ground state
ground_temp = hermiteH(0, X) .* exp(-(X - X0).^2 / (2 * sigma^2));
ground = ground_temp / sqrt(ground_temp * ground_temp');

% Initial state (normalized ground state)
initial_state_temp = hermiteH(0, X) .* exp(-(X - X0).^2 / (2 * sigma^2));
initial_state = initial_state_temp / sqrt(initial_state_temp * initial_state_temp');

% Store results for each A
psis = cell(1, length(A_values));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Main simulation loop over different A values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for idx = 1:length(A_values)
    A = A_values(idx);        % Current perturbation amplitude
    UT = exp(-1i * (P.^2 / 2) * dt);  % Momentum space propagator

    % Reset psi_0 for each simulation
    psi_0 = initial_state;
    fprintf('Simulating for A = %.2f\n', A);

    for m = 1:M
        % Evolution operator in position space
        UV = exp(-1i * ((X.^2) / 2 + A * cos(w * dt * (m - 1)) * sin(X)) * dt / 2);
        
        % Time evolution steps
        psi_1 = UV .* psi_0;
        phi_2 = fft(psi_1);
        phi_3 = UT .* phi_2;
        psi_3 = ifft(phi_3);
        psi_4 = UV .* psi_3;
        psi_0 = psi_4;            % Update wavefunction

        % Check normalization (optional)
        norm_psi = sum(abs(psi_0).^2);
        if norm_psi > 1 + 1e-5 || norm_psi < 1 - 1e-5
            disp('Normalization error during evolution.');
        end
    end

    % Store the final state for the current A
    psis{idx} = psi_0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Plotting the final state profiles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if plot_final_state
    figure('Name', 'Final State (Comparing A)', 'Color', 'w'); % White background
    hold on;

    % Plot each psi for corresponding A with distinct styles and colors
    colors = {
        [0, 0.447, 0.741],  % Blue
        [0.850, 0.325, 0.098],  % Red
        [0.466, 0.674, 0.188],  % Green
        [0.929, 0.694, 0.125, 0.8]   % Yellow-orange 
    };
    linestyles = {'-', '--', ':', ':'};
    for idx = 1:length(A_values)
        plot(X(1:N), abs(psis{idx}(1:N)).^2, 'LineWidth', 2.9, ...
             'Color', colors{idx}, 'LineStyle', linestyles{idx});
    end

    % Customize axes
    xlim([-x_lim, x_lim]);
    ylim([0, max(cellfun(@(psi) max(abs(psi).^2), psis)) * 1.1]);
    xlabel('Position (x)', 'FontSize', fs, 'Interpreter', 'latex');
    ylabel('$|\psi(x)|^2$', 'FontSize', fs, 'Interpreter', 'latex');
    
    % Reduce the number of y-axis ticks
    yticks([0, 0.01, 0.02, 0.03, 0.04]); % 5 ticks

    % Enable grid and frame (axes box)
    grid on;
    box on; % Add a frame to the plot

    % Add legend
    legend(arrayfun(@(A) sprintf('A = %.2f', A), A_values, 'UniformOutput', false), ...
        'FontSize', fs, 'Location', 'northeast', 'Box', 'on');
    
    annotation('textbox', [0.241, 0.621, 0.1, 0.1], ...
               'String', sprintf('$\\omega = %.2f$, $T = %.0f \\pi, M = 10^5$', w, T / pi), ...
               'Interpreter', 'latex', ...
               'FontSize', fs, ...
               'LineStyle', '-', ... % Add border around the box
               'EdgeColor', 'k', ... % Black border
               'BackgroundColor', [0.9, 0.9, 0.9], ... % Light grey
               'HorizontalAlignment', 'center', ... % Center text horizontally
               'VerticalAlignment', 'middle');      % Center text vertically

    % Final adjustments
    pbaspect([2 1 1]); % Set aspect ratio
    set(gca, 'FontSize', fs, 'LineWidth', 1.2); % Optimize tick marks

    % Save the figure
    exportgraphics(gcf, 'Final_State_Comparison_A.png', 'Resolution', 300);
    hold off;
end

fprintf('Simulation comparing A completed.\n');
toc;
