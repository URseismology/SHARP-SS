function model = create_initial_model_transdimensional(prior)
% model contains:
% location, amplitude, width: 2-element vectors
% Sigma, NoiseCorr, NoiseCorr2: scalars
%
% prior contains:
% AmpChange, WidthChange, LocChange, SigChange, CorrChange, CorrChange_2,
% tlen (half), dt, std_P,
% minWid, maxWid,

% gaussian parameters
model.location = [round(rand * (prior.tlen - prior.minLoc / prior.dt) + prior.minLoc / prior.dt) inf];
model.amplitude = [(prior.maxAmp - prior.minAmp) * rand + prior.minAmp inf];
model.width = [round(rand * (prior.maxWid - prior.minWid)) + prior.minWid inf];

% make sure the Gaussian does not fall off the edges
if model.location(1) < (model.width(1) + 1) / 2
    model.location(1) = ceil((model.width(1) + 1) / 2);
elseif model.location(1) > prior.tlen - (model.width(1) + 1) / 2
    model.location(1) = prior.tlen - ceil((model.width(1) + 1) / 2);
end

% noise parameters
model.Sigma = prior.std_P;
model.NoiseCorr  = 0.25;
model.NoiseCorr2 = 1.40;