function G_model = create_G_from_model_transdimensional(model, prior)
% Create G in the time domain (as a time series)
%
% model contains:
% location, amplitude, width: 2-element vectors
% Sigma, NoiseCorr, NoiseCorr2: scalars
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
gauss_1 = model.amplitude(1) * gausswin(round(model.width(1)));

G_model(prior.tlen - (center_width - 1) / 2 : prior.tlen + (center_width - 1) / 2) = gauss_0;
G_model(round(prior.tlen + model.location(1) - model.width(1) / 2) : ...
    round(prior.tlen + model.location(1) - model.width(1) / 2) + round(model.width(1)) - 1) = gauss_1;
G_model(round(prior.tlen - model.location(1) - model.width(1) / 2) : ...
    round(prior.tlen - model.location(1) - model.width(1) / 2) + round(model.width(1)) - 1) = -gauss_1;

% if gauss 2 exist, build it
if ~isinf(model.location(2))
    gauss_2 = model.amplitude(2) * gausswin(round(model.width(2)));
    G_model(round(prior.tlen + model.location(2) - model.width(2) / 2) : ...
        round(prior.tlen + model.location(2) - model.width(2) / 2) + round(model.width(2)) - 1) = gauss_2;
    G_model(round(prior.tlen - model.location(2) - model.width(2) / 2) : ...
        round(prior.tlen - model.location(2) - model.width(2) / 2) + round(model.width(2)) - 1) = -gauss_2;
end