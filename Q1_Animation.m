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
T = 40*pi; %pi*pi*20; %pi*pi ; %5;                  % Time duration of the evolution
M = 10^5; %10^3;                % Total No. of steps in the evolution
dt = T/M;                       % Time step

% Set system and perturbation parameters
n_eigenstate = 0; % 0 is ground state, 1 is first excited state, etc.

A = 0.1;
w = 1.1; %1.2;
w0 = 1.0;

% Set and display plots parameters
plot_initial_state = false;
plot_animation = true;
plot_animation_gif = false;

plot_final_state = false;
plot_transition_probabilities = false;
plot_theory_RWA = false;
plot_theory_noRWA = false;


frames = 200;
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
ground_state = ground_temp /  sqrt(ground_temp * ground_temp');

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
    plot(X(1:N), abs(psi(1:N)).^2, 'LineWidth', 2.9, 'Color', 'b'); 

    % Customize axes limits
    xlim([-x_lim, x_lim]);  % Fix the X-axis limits
    ylim([0, max(abs(psi).^2) * 1.1]); % Add padding to the Y-axis limits

    % Add title and axis labels with consistent font properties
    %title('Initial State', 'FontSize', 16, 'FontWeight', 'normal'); % Match axis font
    xlabel('Position (x)', 'FontSize', 16, 'Interpreter', 'latex');
    ylabel('$|\psi(x)|^2$', 'FontSize', 16, 'Interpreter', 'latex');
    
    annotation('textbox', [0.2, 0.615, 0.1, 0.1], ...
               'String', sprintf('$\\sigma = %.1f$, $x_0 = %.1f$', sigma, X0), ...
               'Interpreter', 'latex', ...
               'FontSize', fs, ...
               'LineStyle', '-', ... % Add border around the box
               'EdgeColor', 'k', ... % Black border
               'BackgroundColor', [0.9, 0.9, 0.9], ... % Light grey
               'HorizontalAlignment', 'center', ... % Center text horizontally
               'VerticalAlignment', 'middle');      % Center text vertically


    % Add a grid for better visual guidance
    grid on;

    % Final adjustments for style
    set(gca, 'FontSize', fs, 'LineWidth', 1.2);
    pbaspect([2 1 1]); % Maintain consistent aspect ratio
    
    % Save the plot
    exportgraphics(gcf, sprintf('Q2_initial_state.png', A), 'Resolution', 300);

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
    h = plot(X(1:N), abs(psi(1:N)).^2, 'LineWidth', 2.9, 'Color', 'b'); 

    % Customize axes limits
    xlim([-x_lim, x_lim]);  % Fix the X-axis limits
    ylim([0, max(abs(psi).^2) * 1.1]); % Add padding to the Y-axis limits

    % Add title and axis labels with consistent font properties
    title('Wavepacket Evolution', 'FontSize', fs, 'FontWeight', 'normal'); % Match axis font
    xlabel('Position (x)', 'FontSize', fs, 'Interpreter', 'latex');
    ylabel('$|\psi(x)|^2$', 'FontSize', fs, 'Interpreter', 'latex');

    % Add a grid for better visual guidance
    grid on;
    pbaspect([2 1 1]); % Maintain consistent aspect ratio

    % Optimize tick marks and their size
    set(gca, 'FontSize', 16, 'LineWidth', 1.2);

    % Refresh the display
    drawnow;
end

% Initialize plot
if plot_animation_gif
    % Create a new figure
    fig = figure('Name', 'Wavepacket Evolution', 'Color', 'w'); % White background for clarity
    filename = 'Wavepacket_Evolution_Optimized.gif'; % File name for the GIF

    % Downsample data
    downsample_factor = 100; % Reduce data size for faster plotting
    X_downsampled = X(1:downsample_factor:end);

    % Initialize the plot
    h = plot(X_downsampled, abs(psi(1:downsample_factor:end)).^2, 'LineWidth', 2.9, 'Color', 'b');
    xlim([-x_lim, x_lim]);
    ylim([0, max(abs(psi).^2) * 1.1]);
    xlabel('Position (x)', 'FontSize', fs, 'Interpreter', 'latex');
    ylabel('$|\psi(x)|^2$', 'FontSize', fs, 'Interpreter', 'latex');
    title('Wavepacket Evolution', 'FontSize', fs, 'FontWeight', 'normal');
    grid on;
    set(gca, 'FontSize', fs, 'LineWidth', 1.2);

    % Time evolution loop with frame skipping
    frame_step = 10; % Update every 10 iterations
    for m = 1:frame_step:M
        % Update psi
        UV = exp(-1i * ((X.^2) / 2 + A * cos(w * dt * (m - 1)) * sin(X)) * dt / 2);
        psi_1 = UV .* psi;
        phi_2 = fft(psi_1);
        phi_3 = UT .* phi_2;
        psi_3 = ifft(phi_3);
        psi_4 = UV .* psi_3;
        psi = psi_4; % Update wavefunction

        % Update plot
        set(h, 'YData', abs(psi(1:downsample_factor:end)).^2);
        drawnow;

        % Capture frame
        frame = getframe(fig);
        img = frame2im(frame);
        [imind, cm] = rgb2ind(img, 256);

        % Write to GIF
        if m == 1
            imwrite(imind, cm, filename, 'gif', 'Loopcount', inf, 'DelayTime', 0.1);
        else
            imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.1);
        end
    end
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
w_diff = w0 - w;
transition_prob_theory_RWA = ( (matrix_element / w_diff) * sin( w_diff * time / 2) ) .^ 2; % remove factor of 4

first = ( (matrix_element / (w+w0) ) * sin( (w+w0) * time / 2) ) .^ 2;
second = ((matrix_element / (w0 - w)) * sin( (w0 - w) * time / 2) ) .^ 2;
third = ((matrix_element).^2 / (w0 - w^2)) * ( (cos(w*time)).^2 - cos(time).*cos(w*time) ); 
transition_prob_theory_noRWA = first + second + third;

%transition_prob_theory = transition_prob_theory_noRWA;

if length(transition_prob_theory_RWA) ~= length(transition_probabilities)
    error('Dimension mismatch: transition_prob_theory and transition_prob must have the same length.');
end

% Plot the transition probability from code and theoretical calculations
if plot_transition_probabilities
    figure('Name', 'Transition Probabilities: Computational vs. Theory', 'Color', 'w');
    hold on;
    
    % Downsample the data for smoother rendering
    downsample_factor = 100; % Adjust this for smoother curves
    time_downsampled = time(1:downsample_factor:end) / pi; % Downsampled time
    
    % Plot computational results
    probabilities_downsampled = transition_probabilities(1:downsample_factor:end);
    plot(time_downsampled, probabilities_downsampled, 'LineWidth', 2.9, ...
         'Color', [0.0, 0.447, 0.741], ... % Blue for computational
         'LineStyle', '-', ...
         'DisplayName', 'Computational (Split-Operator)');

    % Plot theoretical results
    if plot_theory_RWA
        theoretical_probabilities_RWA_downsampled = transition_prob_theory_RWA(1:downsample_factor:end);
        plot(time_downsampled, theoretical_probabilities_RWA_downsampled, 'LineWidth', 2.9, ...
             'Color', [1.0, 0.6, 0.2, 1.0], ... % Warm orange for TDPT + RWA
             'LineStyle', '-', ...
             'DisplayName', 'Theoretical (TDPT + RWA)');
    end
    
    if plot_theory_noRWA
        theoretical_probabilities_noRWA_downsampled = transition_prob_theory_noRWA(1:downsample_factor:end);
        plot(time_downsampled, theoretical_probabilities_noRWA_downsampled, 'LineWidth', 2.9, ...
             'Color', [0.8, 0.2, 0.0, 1.0], ... % Dark red for TDPT, no RWA
             'LineStyle', '--', ...
             'DisplayName', 'Theoretical (TDPT, no RWA)');
    end


    % Customize axes
    xlim([0, T / pi]);
    %ylim([0, max([max(probabilities_downsampled), max(theoretical_probabilities_downsampled)]) * 1.1]);
    %ylim([0,0.0035])

    % Axis labels with LaTeX formatting
    xlabel('Time (multiples of $\pi$)', 'FontSize', fs, 'Interpreter', 'latex');
    ylabel('Transition Probability $P_{1 \leftarrow 0}$', 'FontSize', fs, 'Interpreter', 'latex');
    
    % Add grid, frame, and legend
    grid on;
    box on;
    legend('FontSize', fs, 'Location', 'northeast', 'Box', 'on', 'Interpreter', 'latex');
    annotation('textbox', [0.24, 0.625, 0.1, 0.1], ...
               'String', sprintf('$A = %.2f$, $\\omega = %.2f$, $T = %.0f \\pi$', A, w, T / pi), ...
               'Interpreter', 'latex', ...
               'FontSize', fs, ...
               'LineStyle', '-', ... % Add border around the box
               'EdgeColor', 'k', ... % Black border
               'BackgroundColor', [0.9, 0.9, 0.9], ... % Light grey
               'HorizontalAlignment', 'center', ... % Center text horizontally
               'VerticalAlignment', 'middle');      % Center text vertically


    % Final adjustments for style
    set(gca, 'FontSize', fs, 'LineWidth', 1.2);
    pbaspect([2 1 1]); % Maintain consistent aspect ratio
    
    % Save the plot
    exportgraphics(gcf, sprintf('P10_Comp_v_Theory_RWA+noRWA_A=%.2f.png', A), 'Resolution', 300);
    
    hold off;
    
    fprintf('Simulation completed.\n');
    
end
toc;


