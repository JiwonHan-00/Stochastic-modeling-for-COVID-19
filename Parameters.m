function input = Parameters()
% COVID-19 Contact Tracing Simulation Parameters
% High-risk vs Low-risk Group Model

%% Initial Population Settings
input.SH0 = 8833;        % Initial susceptible high-risk population
input.SL0 = 2490490;     % Initial susceptible low-risk population

% Initial infected compartments
input.E1H0 = 394;        % Initial exposed high-risk (will self-report)
input.E2H0 = 44;         % Initial exposed high-risk (will not self-report)
input.I1H0 = 57;         % Initial infectious high-risk (will self-report)
input.I2H0 = 6;          % Initial infectious high-risk (will not self-report)
input.E1L0 = 144;        % Initial exposed low-risk (will self-report)
input.E2L0 = 16;         % Initial exposed low-risk (will not self-report)
input.I1L0 = 14;         % Initial infectious low-risk (will self-report)
input.I2L0 = 2;          % Initial infectious low-risk (will not self-report)

%% Disease Parameters
% Latent period: Gamma distribution (5.48, 2.72)
input.l1 = (5.48^2)/2.72;
input.l2 = 2.72/5.48;

% Infectious period for non-reporting: Uniform distribution (7.39, 17.39)
input.i1 = 7.39;
input.i2 = 17.39;

% Self-reporting delay: Lognormal distribution
input.aH1 = log(4.04^2/sqrt(27.24 + 4.04^2));  % High-risk group
input.aH2 = sqrt(log(27.24/(4.04^2) + 1));
input.aL1 = log(4.04^2/sqrt(27.24 + 4.04^2));  % Low-risk group
input.aL2 = sqrt(log(27.24/(4.04^2) + 1));

%% Contact and Transmission Parameters
% Pre-social distancing contact rates
input.cHH_pre = 435.8849 * 2.4;    % High-risk to High-risk
input.cHL_pre = 0.4466 * 2.2;      % High-risk to Low-risk
input.cLH_pre = 0.0731 * 2.4;      % Low-risk to High-risk
input.cLL_pre = 8.1729 * 2.2;      % Low-risk to Low-risk

% Post-social distancing contact rates
input.cHH_post = 0.1404 * 2.4;
input.cHL_post = 0.1288 * 2.2;
input.cLH_post = 0.0211 * 2.4;
input.cLL_post = 0.1263 * 2.2;

% Transmission probabilities
input.pHH = 0.721;      % High-risk transmission probability
input.pLL = 0.118;      % Low-risk transmission probability

%% Self-reporting Rates
input.rhoH = 0.4;       % High-risk self-reporting rate
input.rhoL = 0.1;       % Low-risk self-reporting rate

%% Contact Tracing Parameters
% Contact tracing delay: Uniform distribution
input.t1 = 1;           % Minimum delay (days)
input.t2 = 4;           % Maximum delay (days)

% Contact tracing rates (time-dependent)
% Pre-social distancing
input.CTrateH_pre = 0.0008;  % High-risk contact tracing rate
input.CTrateL_pre = 0;       % Low-risk contact tracing rate

% Post-social distancing  
input.CTrateH_post = 1;      % High-risk contact tracing rate
input.CTrateL_post = 0.4;    % Low-risk contact tracing rate

%% Time Parameters
input.SDtime = 11;      % Start time of social distancing (day 11)
input.t_iso = 14;       % Isolation period for infected individuals
input.t_isoC = 14;      % Isolation period for contacts
input.t_suv = 7;        % Contact surveillance period

%% Simulation Parameters
input.Iter_max = round(input.SL0 * 1);  % Maximum iterations
end
