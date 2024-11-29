% Clear workspace and set up
clear all; clc; tic;
close all;
set(0, 'DefaultFigureWindowStyle', 'docked');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Numerical simulation of the evolution of a wavepacket in a 1D harmonic
%   trap using the FFT method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define simulation parameters
a = -20;                      % Left end point
b = +20;                      % Right end point
L = b - a;                    % Width of the space
N = 512;                      % No. of spatial points
X = a + L * (0:N-1) / N;      % Dimensionless coordinates
P = (2 * pi / L) * [0:N/2-1, -N/2:-1]; % Dimensionless momentum

% Time and frequency parameters
T = 60*pi;                      % Total time duration
myArray = [10, 1000, 10^5];   % Array of M values (time steps)
A = 0.05;                     % Perturbation amplitude
w = 1.1;                      % Perturbation frequency
w0 = 1.0;                     % Harmonic oscillator frequency

% Plot parameters
plot_initial_state = true;
plot_animation = false;
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

% First excited state
excited_temp = hermiteH(1, X) .* exp(-(X - X0).^2 / (2 * sigma^2));
excited = excited_temp / sqrt(excited_temp * excited_temp');

% Initial state (normalized ground state)
initial_state_temp = hermiteH(0, X) .* exp(-(X - X0).^2 / (2 * sigma^2));
initial_state = initial_state_temp / sqrt(initial_state_temp * initial_state_temp');

psi_0 = initial_state;
fprintf('Initial state prepared and normalized.\n');

% Verify normalization
tol = 1e-5;
if sum(abs(initial_state).^2) > 1 + tol || sum(abs(initial_state).^2) < 1 - tol
    error('Normalization error in the initial state.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Main simulation loop over different M values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
psis = cell(1, length(myArray));  % Store final states for each M

for idx = 1:length(myArray)
    M = myArray(idx);            % Number of time steps
    dt = T / M;                  % Time step
    UT = exp(-1i * (P.^2 / 2) * dt);  % Momentum space propagator

    % Reset psi_0 for each simulation
    psi_0 = initial_state;
    fprintf('Simulating for M = %d\n', M);

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

        % Check normalization
        if sum(abs(psi_0).^2) > 1 + tol || sum(abs(psi_0).^2) < 1 - tol
            disp('Normalization error during evolution.');
        end

        % Optional: Plot animation
        if plot_animation && mod(m, M / frames) == 0
            figure(1);
            plot(X(1:N), abs(psi_0(1:N)).^2, 'LineWidth', 2, 'Color', 'b');
            xlim([-x_lim, x_lim]);
            ylim([0, max(abs(psi_0).^2) * 1.1]);
            drawnow;
        end
    end

    % Store the final state for the current M
    psis{idx} = psi_0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Plotting the final state profiles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting the final state profiles
if plot_final_state
    figure('Name', 'Final State', 'Color', 'w');
    hold on;

    % Plot each psi for corresponding M with distinct styles
    colors = {[0.466, 0.674, 0.188], [0, 0.447, 0.741], [0.850, 0.325, 0.098]};
    linestyles = {'-', '-', '--'};
    for idx = 1:length(myArray)
        plot(X(1:N), abs(psis{idx}(1:N)).^2, 'LineWidth', 2.9, ...
             'Color', colors{idx}, 'LineStyle', linestyles{idx});
    end

    % Customize axes
    xlim([-x_lim, x_lim]);
    ylim([0, max(cellfun(@(psi) max(abs(psi).^2), psis)) * 1.1]);
    xlabel('Position (x)', 'FontSize', fs, 'Interpreter', 'latex');
    ylabel('$|\psi(x)|^2$', 'FontSize', fs, 'Interpreter', 'latex');
    grid on;
    box on; % Add a frame to the plot
    yticks([0, 0.01, 0.02, 0.03, 0.04]);


    % Add legend
    legend(arrayfun(@(M) sprintf('M = %d', M), myArray, 'UniformOutput', false), ...
        'FontSize', fs, 'Location', 'northeast', 'Box', 'on');
    annotation('textbox', [0.241, 0.621, 0.1, 0.1], ...
               'String', sprintf('$A = %.2f, \\omega = %.2f$, $T = %.0f \\pi$', A, w, T / pi), ...
               'Interpreter', 'latex', ...
               'FontSize', fs, ...
               'LineStyle', '-', ... % Add border around the box
               'EdgeColor', 'k', ... % Black border
               'BackgroundColor', [0.9, 0.9, 0.9], ... % Light grey
               'HorizontalAlignment', 'center', ... % Center text horizontally
               'VerticalAlignment', 'middle');      % Center text vertically
    % Final adjustments
    pbaspect([2 1 1]);
    set(gca, 'FontSize', fs, 'LineWidth', 1.2);

    % Save the figure
    exportgraphics(gcf, 'Final_State_Comparison_M.png', 'Resolution', 300);
    hold off;
end
