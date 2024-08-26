function model_new = create_new_model(model, action, prior)
% model_new = create_new_model(model, action, prior)
%
% Create a new model, given the current model, prior, and the chosen
% action.
%
% model and model_new has the same structure:
% location, amplitude, and width of the Gaussian pair (three scalars since
% they're symmetrical), noise standard deviation (sigma), noise correlation
% parameters (two of them)
%
% model (and model_new) contains:
% location, amplitude, width, Sigma, NoiseCorr, NoiseCorr2
%
% list of actions:
% (1) Change the location of the Gaussian pair: slightly
% (2) Change the location of the Gaussian pair: randomly to another location
% (3) Change the amplitude of the Gaussian pair
% (4) Change the width of the Gaussian pair
% (5) Change the noise standard deviation (Sigma)
% (6) Change the noise correlation parameters (two params, but changing 
% together)
%
% prior contains:
% AmpChange, WidthChange, LocChange, SigChange, CorrChange, CorrChange_2,
% tlen (half), dt, std_P, minWid, maxWid
%
% Author: Ziqi Zhang (ziqi.zhang@rochester.edu)

model_new = model;

switch action

    case 1 % change location slightly

        model_new.location = model_new.location + round(prior.LocChange * randn/prior.dt);

        % make sure the Gaussian does not fall off the edges
        while model_new.location < ((model_new.width + 1) / 2 + prior.minLoc / prior.dt) || ...
                model_new.location > prior.tlen - (model_new.width + 1) / 2
            model_new.location = model_new.location + round(prior.LocChange * randn/prior.dt);
        end

    case 2 % change location entirely

        model_new.location = round(rand * (prior.tlen - prior.minLoc / prior.dt) + prior.minLoc / prior.dt);

        % make sure the Gaussian does not fall off the edges
        while model_new.location < ((model_new.width + 1) / 2 + prior.minLoc / prior.dt) || ...
                model_new.location > prior.tlen - (model_new.width + 1) / 2
            model_new.location = round(rand * (prior.tlen - prior.minLoc / prior.dt) + prior.minLoc / prior.dt);
        end

    case 3 % change amplitude

        if rand < 0.8
            new_amp = model_new.amplitude + prior.AmpChange * randn;
        else
            new_amp = -model_new.amplitude;
        end

        % make sure amplitude is within range
        model_new.amplitude = min(max(new_amp, prior.minAmp), prior.maxAmp);

    case 4 % change width

        new_wid = model_new.width + round(prior.WidthChange * randn/prior.dt);

        % make sure the Gaussian does not fall off the edges
        while model_new.location < (new_wid + 1) / 2 || ...
                model_new.location > prior.tlen - (new_wid + 1) / 2 || ...
                new_wid < prior.minWid || new_wid > prior.maxWid
            new_wid = model_new.width + round(prior.WidthChange * randn/prior.dt);
        end

        model_new.width = new_wid;

    case 5 % change sigma

        model_new.Sigma = ...
            max(eps, model_new.Sigma + prior.SigChange * randn * prior.std_P);

    case 6 % change noise correlation parameters

        model_new.NoiseCorr = ...
            max(eps, model_new.NoiseCorr + prior.CorrChange * randn);
        model_new.NoiseCorr2 = ...
            min(1, model_new.NoiseCorr2 + prior.CorrChange_2 * randn);

end


% %%%%%%%%%%%%%%%%%%%%%%%% debug: visualize
% time = - (prior.tlen - 1) * prior.dt : prior.dt : (prior.tlen - 1) * prior.dt;
% G_model = create_G_from_model(model_new, prior);
% figure(1);
% clf;
% plot(time, G_model);
% xlim([min(time) max(time)]);
% ylim([-0.5 1]);
