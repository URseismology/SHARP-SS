function G_model = create_G_from_model(model, prior)
% Create G in the time domain (as a time series)
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

% using gausswin to generate gaussians
gauss_0 = center_amplitude * gausswin(center_width);
gauss_1 = model.amplitude * gausswin(round(model.width));

G_model(prior.tlen - (center_width - 1) / 2 : prior.tlen + (center_width - 1) / 2) = gauss_0;
G_model(round(prior.tlen + model.location - model.width / 2) : ...
    round(prior.tlen + model.location - model.width / 2) + round(model.width) - 1) = gauss_1;
G_model(round(prior.tlen - model.location - model.width / 2) : ...
    round(prior.tlen - model.location - model.width / 2) + round(model.width) - 1) = -gauss_1;