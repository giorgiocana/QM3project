%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Numerical simulation of the evolution of a wavepacket in a 1D harmonic
%   trap using fast fourier transport (fft) method
%   In matlab, there is a internal function fft(...) which can do the job
%   For more detials on fft, type "help fft" in matlab command window and
%   press enter button to check the matlab help file  

% Unit of energy: hbar*omega, where h_bar is the Planck constant and
%   omega is the frequency of the trap

%   Unit of length: l=sqrt(h_bar/(m*omega)), where sqrt(...) is the square
%  root function and m is the mass of the particle

%   Unit of momentum: hbar/l
%    energy unit: hbar\omega,  Hamiltonian --> dimensionless

%%   time dimensionless: omega*t    i d/dt | >= dimension H |>
%    dimensionless time = 2pi. one classical period
clear all; clc; tic
close all;
set(0, 'DefaultFigureWindowStyle', 'docked')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Define initial parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------
a = -20;                       % Left end point 
b = +20;                       % Right end point 
L = b-a;                        % Width of the space
N = 512;                       % No. of cells
X = a+L*(0:N-1)/N;                % Dimensionless coordinates
P = (2*pi/L)*[0:N/2-1,-N/2:-1]; % Dimensionless momentum

% Set DVR parameters 
T = 200; %pi*pi*20; %pi*pi ; %5;                  % Time duration of the evolution
M = 10^3; %10^3;                % Total No. of steps in the evolution
dt = T/M;                       % Time step

% Set system and perturbation parameters
n_eigenstate = 0; % 0 is ground state, 1 is first excited state, etc.

A = 0.05;
w = 1.1; %1.2;
w0 = 1.0;

% Set and display plots parameters
plot_initial_state = true;
plot_animation = true;
plot_final_state = true;
plot_transition_probabilities = true;
plot_theory = true;

frames = 100;
fs = 22;
x_lim = 20;

fprintf('\n Initial state: %d', n_eigenstate);
fprintf('\n A = %d', A);
fprintf('\n w = %d', w);
fprintf('\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Define vectors to store split step propagators in position and
%   momentum space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NOTE: hbar=1 in our dimensionless units
UT = exp(-1i*(P.^2/2)*dt);  % One-setp propagator in momentum space

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Define the initial state psi_0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Define the initial state
X0 = 0.0; % wavepacket / Gaussian center
sigma=1/sqrt(w0);  % sigma is the width of the initial wavepacket / Gaussian

% Prepare ground and excited state for comparison
ground_temp = (hermiteH(0,X)).*exp(-(X-X0).^2/(2*sigma^2));
ground = ground_temp /  sqrt(ground_temp * ground_temp');

excited_temp = (hermiteH(1,X)).*exp(-(X-X0).^2/(2*sigma^2));
excited = excited_temp / sqrt(excited_temp * excited_temp');

%   Normalized initial state as ground state or excited of harmonic oscillators
Poly=hermiteH(n_eigenstate,X);  % this is to create hermite polymonial)
initial_state_temp = Poly.*exp(-(X-X0).^2/(2*sigma^2));%1i means "i"
initial_state = initial_state_temp/sqrt(initial_state_temp * initial_state_temp'); %normalization

psi_0 = initial_state; %/sqrt(sum(abs(psiprep).^2));%normalized state
psi = psi_0;
fprintf('State prep done\n');

% Check normalisation
tol = 1e-15;
if sum(abs(initial_state).^2) < 1 + tol && sum(abs(initial_state).^2) > 1 - tol 
    disp('State is properly normalised')
end
%Demofftp=fft(VE_INI);

%   Verify that the initial state is indeed normalized
fprintf('Abs_val_squared = %.10f\n', sum(abs(initial_state).^2) );


%   Check the probability distribution of the initial state
%fprintf('Area = %.10f\n', area(X,abs(initial_state).^2) );

%   Check dimensions
if length(psi_0) ~= length(excited)
    error('Dimension mismatch: psi0 and excited must have the same length.');
    fprintf('\nlength(psi0) = %d\n', length(psi_0));
    fprintf('\nlength(excited) = %d\n', length(excited));
end

% Plot initial state
if plot_initial_state
    % Create a new figure with a specific name
    figure('Name', 'Initial State', 'Color', 'w'); % White background for clarity

    % Plot the initial state with a blue line
    plot(X(1:N), abs(psi(1:N)).^2, 'LineWidth', 2, 'Color', 'b'); 

    % Customize axes limits
    xlim([-x_lim, x_lim]);  % Fix the X-axis limits
    ylim([0, max(abs(psi).^2) * 1.1]); % Add padding to the Y-axis limits

    % Add title and axis labels with consistent font properties
    title('Initial State', 'FontSize', 16, 'FontWeight', 'normal'); % Match axis font
    xlabel('Position (x)', 'FontSize', 16, 'Interpreter', 'latex');
    ylabel('$|\psi(x)|^2$', 'FontSize', 16, 'Interpreter', 'latex');

    % Add a grid for better visual guidance
    grid on;

    % Optimize tick marks and their size
    set(gca, 'FontSize', 16, 'LineWidth', 1.2);

    % Refresh the display
    drawnow;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Evolve the state 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define transition probabilities array
transition_probabilities = zeros(1, M);

% Initialize plot
if plot_animation
    % Create a new figure with a specific name
    figure('Name', 'Wavepacket Evolution', 'Color', 'w'); % White background for clarity

    % Initialize the plot handle with a blue line
    h = plot(X(1:N), abs(psi(1:N)).^2, 'LineWidth', 2, 'Color', 'b'); 

    % Customize axes limits
    xlim([-x_lim, x_lim]);  % Fix the X-axis limits
    ylim([0, max(abs(psi).^2) * 1.1]); % Add padding to the Y-axis limits

    % Add title and axis labels with consistent font properties
    title('Wavepacket Evolution', 'FontSize', fs, 'FontWeight', 'normal'); % Match axis font
    xlabel('Position (x)', 'FontSize', fs, 'Interpreter', 'latex');
    ylabel('$|\psi(x)|^2$', 'FontSize', fs, 'Interpreter', 'latex');

    % Add a grid for better visual guidance
    grid on;

    % Optimize tick marks and their size
    set(gca, 'FontSize', 16, 'LineWidth', 1.2);

    % Refresh the display
    drawnow;
end


for m = 1:M
        UV = exp(-1i*( (X.^2)/2 + A*cos(w*dt*(m-1))*sin(X) )*dt/2);

        psi_1 = UV.*psi_0;
        phi_2 = fft(psi_1);   %wavefunction in momentum space
        phi_3 = UT.*phi_2;
        psi_3 = ifft(phi_3);
        psi_4 = UV.*psi_3;
        psi_0 = psi_4; %prepare a new cycle 
        
        % Re-normalize psi_0 to avoid cumulative drift
        %psi_0 = psi_0_temp / sqrt(psi_0_temp * psi_0_temp');
        if sum(abs(psi_0).^2) > 1 + tol && sum(abs(psi_0).^2) < 1 - tol 
            disp('NORMALISATION ERROR')
        end

        % Calculate transition probability to the first excited state
        transition_probabilities(m) = abs(dot(conj(excited), psi_0))^2;
    
        if plot_animation
            if mod(m, M/frames) == 0
                set(h, 'YData', abs(psi_0(1:N)).^2);  % Update the plot data
                drawnow;                              % Force the plot update
            end
        end
end

psi=psi_0; %final state updated 

%plotting the final state profile
if plot_final_state
    % Create a new figure with a specific name
    figure('Name', 'Final State', 'Color', 'w'); % White background for clarity

    % Plot the final state with a blue line
    plot(X(1:N), abs(psi(1:N)).^2, 'LineWidth', 2, 'Color', 'b'); 

    % Customize axes limits
    xlim([-x_lim, x_lim]);  % Fix the X-axis limits
    ylim([0, 0.045]); % Add padding to the Y-axis limits
    %ylim([0, max(abs(psi).^2) * 1.1]); % Add padding to the Y-axis limits

    % Add title and axis labels with the same font properties
    %title('Final State', 'FontSize', 16, 'FontWeight', 'normal'); % Match axis font
    xlabel('Position (x)', 'FontSize', fs, 'Interpreter', 'latex');
    ylabel('$|\psi(x)|^2$', 'FontSize', fs, 'Interpreter', 'latex');

    % Add a grid for better visual guidance
    grid on;

    % Add a legend inside the plot area
    legend(sprintf('A = %.2f, w = %.2f, T = %d, M = %d', A, w, T, M), 'FontSize', fs, ...
           'Location', 'northeast', 'Box', 'off'); % Position it in the top right corner

    % Adjust figure aspect ratio (optional)
    pbaspect([3 1 1]); % Uncomment if you want a specific aspect ratio

    % Optimize tick marks and their size
    set(gca, 'FontSize', fs, 'LineWidth', 1.2);

    % Save the figure as a high-quality PNG
    exportgraphics(gcf, sprintf('Final_state_A=%.3f_w=%.3f_M=%d.png', A, w, M), 'Resolution', 300);

    % Refresh the display
    drawnow;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  COMPARE SIMULATION WITH TIME-DEPENDENT PERTURBATION THEORY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

time = (0:M-1) * dt;

% Display maximum transition probability (optional)
[max_prob, max_idx] = max(transition_probabilities);
fprintf('Maximum transition probability: %.7f at time %.4f pi \n', max_prob, time(max_idx)/pi);

% Time evolution using the perturbation theory formula
matrix_element = A / (sqrt(2)*exp(1/4));
fprintf('Vnm/A = %.4f \n', matrix_element/A);
w_diff = w - w0;
transition_prob_theory = ( (matrix_element / w_diff) * sin( w_diff * time / 2) ) .^ 2; % remove factor of 4

if length(transition_prob_theory) ~= length(transition_probabilities)
    error('Dimension mismatch: transition_prob_theory and transition_prob must have the same length.');
end

% Plot the transition probability from code and theoretical calculations% Plot the transition probability from code and theoretical calculations
if plot_transition_probabilities
    % Create a new figure with a specific name
    figure('Name', 'Transition Probabilities', 'Color', 'w'); % White background for clarity

    % Hold the figure to overlay multiple plots
    hold on;

    % Plot numerical transition probabilities
    plot(time / pi, transition_probabilities, 'LineWidth', 2, 'Color', 'b', ...
         'DisplayName', 'Numerical (Code)');

    % Plot theoretical transition probabilities if enabled
    if plot_theory
        plot(time / pi, transition_prob_theory, 'LineWidth', 2, 'Color', 'r', ...
             'DisplayName', 'Theoretical (TDPT)');
    end

    % Customize legend inside the plot area
    legend('FontSize', 14, 'Location', 'northeast', 'Box', 'on'); % Place in top-right corner

    % Add title and axis labels with consistent font properties
    title('Transition Probability to First Excited State', 'FontSize', fs, 'FontWeight', 'normal');
    xlabel('Time (multiples of $\pi$)', 'FontSize', fs, 'Interpreter', 'latex');
    ylabel('Transition Probability', 'FontSize', fs, 'Interpreter', 'latex');

    % Add a grid for better visual guidance
    grid on;

    % Optimize tick marks and their size
    set(gca, 'FontSize', fs, 'LineWidth', 1.2);

    exportgraphics(gcf, 'Transition_Probabilities.png', 'Resolution', 300);

    % Release the hold on the figure
    hold off;
end

