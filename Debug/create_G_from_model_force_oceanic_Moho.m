function G_model = create_G_from_model_force_oceanic_Moho(model, prior)
% Create G in the time domain (as a time series)
% Debug version: force an oceanic Moho pair
% H = 7.0 km, vs = 3.5 km/s --> location at 4.0 s
% Assume amplitude of 0.2
%
% model contains:
% location, amplitude, width, Sigma, NoiseCorr, NoiseCorr2
%
% prior contains:
% AmpChange, WidthChange, LocChange, SigChange, CorrChange, CorrChange_2,
% tlen (half), dt, std_P,
% minWid, maxWid,
% 

G_model = zeros(1, prior.tlen * 2 - 1); % because tlen is half of the full length

center_amplitude = 1.0; % unit amplitude
center_width     = 3; % in sample points

Moho_location  = round(4.0 / prior.dt);
Moho_amplitude = 0.3;
Moho_width     = 10; % in sample points

% using gausswin to generate gaussians
gauss_0 = center_amplitude * gausswin(center_width);
gauss_1 = Moho_amplitude * gausswin(Moho_width);
gauss_2 = model.amplitude * gausswin(round(model.width));

% center
G_model(prior.tlen - (center_width - 1) / 2 : prior.tlen + (center_width - 1) / 2) = gauss_0;

% Moho
G_model(round(prior.tlen + Moho_location - Moho_width / 2) : ...
    round(prior.tlen + Moho_location - Moho_width / 2) + round(Moho_width) - 1) = -gauss_1;
G_model(round(prior.tlen - Moho_location - Moho_width / 2) : ...
    round(prior.tlen - Moho_location - Moho_width / 2) + round(Moho_width) - 1) = gauss_1;

% modeled upper mantle
G_model(round(prior.tlen + model.location - model.width / 2) : ...
    round(prior.tlen + model.location - model.width / 2) + round(model.width) - 1) = gauss_2;
G_model(round(prior.tlen - model.location - model.width / 2) : ...
    round(prior.tlen - model.location - model.width / 2) + round(model.width) - 1) = -gauss_2;