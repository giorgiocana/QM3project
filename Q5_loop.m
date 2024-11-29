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
T = pi*pi*20; %pi*pi ; %5;                  % Time duration of the evolution
M = 10^4; %10^3;                % Total No. of steps in the evolution
dt = T/M;                       % Time step

% Set system and perturbation parameters
n_eigenstate = 0; % 0 is ground state, 1 is first excited state, etc.

number_of_trials = 50;
A_vector = (1:number_of_trials).*0.005;
w = 0.75; %or Q6 2.1;
w0 = 1.0;

% Set and display plots parameters
plot_initial_state = false;
plot_animation = false;
plot_final_state = false;
plot_transition_probabilities = false;
plot_theory = false;

frames = 100;
x_lim = 20;

fprintf('\n Initial state: %d', n_eigenstate);
%fprintf('\n A = %d', A);
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

%   Plot initial state
if plot_initial_state
    figure('Name','Initial State');
    disp(length(psi))
    plot(X(1:N),abs(psi(1:N)).^2);   % plotting initial state
    xlim([-x_lim,x_lim]);           % Fix the X-axis limits
    ylim([-0.001, max(abs(psi).^2)]);           % Fix the Y-axis limits
    title('Initial State |psi|^2');
    xlabel('Position');
    ylabel('|psi|^2');
    drawnow;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Evolve the state 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define transition probabilities array
transition_probabilities = zeros(1, M);
max_probabilities = zeros(1, number_of_trials);

% Initialize plot
if plot_animation
    figure('Name','Animation');
    h = plot(X(1:N), abs(psi(1:N)).^2);    % Initialize the plot handle
    xlim([-x_lim,x_lim]);           % Fix the X-axis limits
    ylim([-0.001, max(abs(psi).^2)]);           % Fix the Y-axis limits
    title('Wavepacket Evolution');
    xlabel('Position');
    ylabel('|psi|^2');
    drawnow;
end


% Loop over different values of A
for i = (1:number_of_trials) 
    A = A_vector(i);
    psi_0 = initial_state;
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
    max_probabilities(i) = max(transition_probabilities);
    disp([A,max_probabilities(i)]);
end

psi=psi_0; %final state updated 

%plotting the final state profile
if plot_final_state
    figure('Name','Final State');
    plot(X(1:N), abs(psi(1:N)).^2);    % Initialize the plot handle
    xlim([-x_lim,x_lim]);           % Fix the X-axis limits
    ylim([-0.001, max(abs(psi).^2)]);           % Fix the Y-axis limits
    title('Final state');
    xlabel('Position');
    ylabel('|psi|^2');
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

% Plot the transition probability from code and theoretical calculations
if plot_transition_probabilities
    figure('Name','Transition Probabilities');
    hold on;
        plot(time/pi, transition_probabilities, 'LineWidth', 1.5, 'DisplayName', 'Numerical (Code)');
        if plot_theory
            plot(time/pi, transition_prob_theory, 'LineWidth', 1.5, 'DisplayName', 'Theoretical (TDPT)');
        end
    legend
    title('Transition Probability to First Excited State');
    xlabel('Time, multiples of PI');
    ylabel('Transition Probability');
    grid on;
    hold off;
end




% Giorgio's code
fs = 22; % this is the fontsize I used, seems to be working pretty well on latex

% Plot the transition amplitude against A
figure('Name','Transition Probability against A');
hold on;
scatter(log(A_vector), log(max_probabilities), 'LineWidth', 2.2, 'Color', 'b', 'DisplayName', 'Numerical (Split-Operator)');

% Calculate best-fit line coefficients: slope (m) and intercept (b)
coeffs = polyfit(log(A_vector), log(max_probabilities), 1);
m = coeffs(1); % Slope
b = coeffs(2); % Intercept

% Generate the best-fit line
y_fit = m * log(A_vector) + b;
plot(log(A_vector), y_fit, 'LineWidth', 2.2, 'Color', 'r', 'DisplayName', 'Best Fit')

xlabel('$\ln A$', 'FontSize', fs, 'Interpreter', 'latex');
ylabel('$\ln(\mathrm{Max} ~P_{0\to 1})$', 'FontSize', fs, 'Interpreter', 'latex');
grid on;
pbaspect([2 1 1]); % Uncomment if you want a specific aspect ratio
% Optimize tick marks and their size
set(gca, 'FontSize', fs, 'LineWidth', 1.2);

r=corr(log(A_vector)', log(max_probabilities)');
disp(r)
fprintf('r-value=% \n', r);

legend('Location', 'southeast')
%exportgraphics(gcf, sprintf('Q5_Transition_Growth_w=%.3f_m=%.5f_r=%.5f.png', w, m, r), 'Resolution', 300);
hold off;

figure('Name','Transition growth with A');
plot(A_vector, max_probabilities);    % Initialize the plot handle
%xlim([-x_lim,x_lim]);           % Fix the X-axis limits
%ylim([-0.001, max(abs(psi).^2)]);           % Fix the Y-axis limits
title('Transition growth with A');
xlabel('A');
ylabel('Maximum transition probability');
drawnow;