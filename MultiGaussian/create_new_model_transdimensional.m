function model_new = create_new_model_transdimensional(model, action, prior)
% model_new = create_new_model(model, action, prior)
%
% Create a new model, given the current model, prior, and the chosen
% action.
%
% model and model_new has the same structure:
% location, amplitude, and width of the Gaussian pair (three 2-element
% vectors), noise standard deviation (sigma), noise correlation
% parameters (two of them)
%
% model (and model_new) contains:
% location, amplitude, width: 2-element vectors
% Sigma, NoiseCorr, NoiseCorr2: scalars
%
% list of actions:
% (1) Add (if currently 1) or remove (if currently 2) one Gaussian pair
% (2) Change the location of one Gaussian pair
% (3) Change the amplitude of one Gaussian pair
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

    case 1 % add or remove a gaussian

        if ~isinf(model.location(2)) % 2 gaussians in current model
            % decide which gaussian pair to remove
            gauss_ind = randi([1, 2], 1);
            % remove
            model_new.location(gauss_ind)  = inf;
            model_new.amplitude(gauss_ind) = inf;
            model_new.width(gauss_ind)     = inf;
            % put the only gaussian in the first element
            model_new.location  = sort(model_new.location);
            model_new.amplitude = sort(model_new.amplitude);
            model_new.width     = sort(model_new.width);
        else % 1 gaussian in current model
            % randomly choose location and width for new gaussian
            model_new.location(2)  = round(rand * (prior.tlen - prior.minLoc / prior.dt) + prior.minLoc / prior.dt);
            model_new.width(2)     = round(rand * (prior.maxWid - prior.minWid)) + prior.minWid;
            model_new.amplitude(2) = (prior.maxAmp - prior.minAmp) * rand + prior.minAmp;
            % make sure the Gaussian does not fall off the edges
            % and does not overlap with the other Gaussian
            [~, g1_ind] = min(model_new.location); % find the first gaussian (returns 1 if only one because the other one is inf)
            if (model_new.location(2) < (model_new.width(2) + 1) / 2 + prior.minLoc / prior.dt) || ...
                    model_new.location(2) > prior.tlen - (model_new.width(2) + 1) / 2 || ...
                    ((model_new.location(g1_ind) + (model_new.width(g1_ind) + 1) / 2) ...
                    >= (model_new.location(3 - g1_ind) - (model_new.width(3 - g1_ind) - 1) / 2))
                model_new = model;
            end
        end
    
    case 2 % change location

        % decide which gaussian pair to make change
        if ~isinf(model.location(2))
            gauss_ind = randi([1, 2], 1);
        else
            gauss_ind = 1;
        end

        % decide change slightly or randomly: 50% chance
        rand_change = 1;
        if rand <= 0.5
            rand_change = 0;
        end

        if rand_change
            model_new.location(gauss_ind) = round(rand * (prior.tlen - prior.minLoc / prior.dt) + prior.minLoc / prior.dt);
        else
            model_new.location(gauss_ind) = ...
                model_new.location(gauss_ind) + round(prior.LocChange * randn/prior.dt);
        end

        % make sure the Gaussian does not fall off the edges
        % and does not overlap with the other Gaussian (if exist)
        [~, g1_ind] = min(model_new.location); % find the first gaussian (returns 1 if only one because the other one is inf)
        if (model_new.location(gauss_ind) < (model_new.width(gauss_ind) + 1) / 2 + prior.minLoc / prior.dt) || ...
                model_new.location(gauss_ind) > prior.tlen - (model_new.width(gauss_ind) + 1) / 2 || ...
                (~isinf(model_new.location(2)) ...
                && ((model_new.location(g1_ind) + (model_new.width(g1_ind) + 1) / 2) ...
                >= (model_new.location(3 - g1_ind) - (model_new.width(3 - g1_ind) - 1) / 2)))
            model_new = model;
        end
    
    case 3 % change amplitude

        % decide which gaussian pair to make change
        if ~isinf(model_new.location(2))
            gauss_ind = randi([1, 2], 1);
        else
            gauss_ind = 1;
        end
        % generate new amplitude
        new_amp = model_new.amplitude(gauss_ind) + prior.AmpChange * randn;
        % make sure amplitude is within range
        model_new.amplitude(gauss_ind) = min(max(new_amp, prior.minAmp), prior.maxAmp);

    case 4 % change width

        % decide which gaussian pair to make change
        if ~isinf(model_new.location(2))
            gauss_ind = randi([1, 2], 1);
        else
            gauss_ind = 1;
        end
        % generate new width
        model_new.width(gauss_ind) = model_new.width(gauss_ind) + round(prior.WidthChange * randn/prior.dt);
        % make sure the Gaussian does not fall off the edges
        % and does not overlap with the other Gaussian if exist
        [~, g1_ind] = min(model_new.location); % find the first gaussian (returns 1 if only one because the other one is inf)
        if (model_new.location(gauss_ind) < (model_new.width(gauss_ind) + 1) / 2 + prior.minLoc / prior.dt) || ...
                model_new.location(gauss_ind) > prior.tlen - (model_new.width(gauss_ind) + 1) / 2 || ...
                (~isinf(model.location(2)) ...
                && ((model_new.location(g1_ind) + (model_new.width(g1_ind) + 1) / 2) ...
                >= (model_new.location(3 - g1_ind) - (model_new.width(3 - g1_ind) - 1) / 2))) || ...
                model_new.width(gauss_ind) < prior.minWid || model_new.width(gauss_ind) > prior.maxWid
            model_new = model;
        end

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
