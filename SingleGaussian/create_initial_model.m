function model = create_initial_model(prior)
% model contains:
% location, amplitude, width, Sigma, NoiseCorr, NoiseCorr2
%
% prior contains:
% AmpChange, WidthChange, LocChange, SigChange, CorrChange, CorrChange_2,
% tlen (half), dt, std_P,
% minWid, maxWid,

% gaussian parameters
model.location = round(rand * (prior.tlen - prior.minLoc / prior.dt) + prior.minLoc / prior.dt);
model.amplitude = (prior.maxAmp - prior.minAmp) * rand + prior.minAmp;
model.width = round(rand * (prior.maxWid - prior.minWid)) + prior.minWid;

% make sure the Gaussian does not fall off the edges
if model.location < model.width / 2
    model.location = ceil(model.width / 2);
elseif model.location > prior.tlen - model.width / 2
    model.location = prior.tlen - ceil(model.width / 2);
end

% noise parameters
model.Sigma = prior.std_P;
model.NoiseCorr  = 0.25;
model.NoiseCorr2 = 1.40;