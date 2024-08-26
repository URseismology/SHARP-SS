function [LocSave, AmpSave, WidSave, SigSave, NC1Save, NC2Save, LikSave, prior] = ...
    THBD_SS(time, P, D, ManArgs, OptArgs)

% read in mandatory and optional arguments
totalSteps = ManArgs.totalSteps;
saveSteps = ManArgs.saveSteps;
output_time_len = ManArgs.output_time;

dt = time(2) - time(1);
output_time = linspace(-output_time_len, output_time_len, ...
    round(2 * output_time_len / dt) + 1);

isBP = OptArgs.isBP;
lf = OptArgs.lf;
hf = OptArgs.hf;
isPlotTrue = OptArgs.isPlotTrue;
trueG_filename = OptArgs.trueG_filename;

totalActions = 3; % actions per step
ChangeCorrPercent = 0.1; % chance of doing action 6 (noise correlation)

prior = define_prior(output_time, P);

% optional zero-phase bandpass filtering
if isBP
    [P, D] = preprocess_data(P, D, lf, hf, prior);
end

% initialize saving vectors
totalSavedSteps = ceil(totalSteps / saveSteps);
LocSave = zeros(totalSavedSteps, 1);
AmpSave = zeros(totalSavedSteps, 1);
WidSave = zeros(totalSavedSteps, 1);
SigSave = zeros(totalSavedSteps, 1);
NC1Save = zeros(totalSavedSteps, 1);
NC2Save = zeros(totalSavedSteps, 1);
LikSave = zeros(totalSavedSteps, 1);
iSaveStep = 1;

% initialize starting model and calculate its likelihood probability
model = create_initial_model(prior);
model.LikeProb = calculate_like_prob(P, D, model, prior, 1);

% start MCMC loop
for iStep = 1:totalSteps

    changedCorr = 0; % flag for CdInvMatrix

    % do multiple actions per step
    for iAction = 1:totalActions

        action = randperm(5, 1); % action 6 decided separately

        % decide if go to action 6
        if rand < ChangeCorrPercent
            action = 6;
            changedCorr = 1; % flag for CdInvMatrix
        end

        % create new model based on the action chosen
        model_new = create_new_model(model, action, prior);

    end

    % calculate likelihood probability
    % only do CdInvMatrix part if noise correlation changes or step 1
    if changedCorr || iStep == 1
        [model_new.LikeProb, D_model, R_LT ,R_UT ,R_P, LogDetR] = ...
            calculate_like_prob(P, D, model_new, prior, 1);
        
        if iStep == 1
            D_model_toplot = D_model;
        end

    else % in other cases skip CdInvMatrix part
        [model_new.LikeProb, D_model, ~, ~, ~, ~] = ...
            calculate_like_prob(P, D, model_new, prior, 0, R_LT ,R_UT ,R_P, LogDetR);

    end

    % calculate probability of acceptance
    ProbAccept = exp(model_new.LikeProb - model.LikeProb);

    % accept model based on the comparison of likelihood probability
    if isfinite(model.LikeProb) && isfinite(model_new.LikeProb)
        if ProbAccept > 1 % || ProbAccept > rand % does this make sense?
            model = model_new;
            D_model_toplot = D_model;
        end
    end

    % if opt for align, pad D_model_toplot
    if prior.align
        D_model_toplot_tmp = D_model_toplot;
        D_model_toplot = NaN(1, length(D));
        D_model_toplot(round(prior.align / prior.dt):round(prior.align / prior.dt) + length(D_model) - 1) =...
            D_model_toplot_tmp(~isnan(D_model_toplot_tmp));
    end

    % save results
    if mod((iStep - 1), saveSteps) == 0

        LocSave(iSaveStep) = model.location;
        AmpSave(iSaveStep) = model.amplitude;
        WidSave(iSaveStep) = model.width;
        SigSave(iSaveStep) = model.Sigma;
        NC1Save(iSaveStep) = model.NoiseCorr;
        NC2Save(iSaveStep) = model.NoiseCorr2;
        LikSave(iSaveStep) = model.LikeProb;

        %%%%%%%% debug - visualize accepted model %%%%%%%%%
        if prior.FOM
            G_model = create_G_from_model_force_oceanic_Moho(model, prior);
        else
            G_model = create_G_from_model(model, prior);
        end

        if isPlotTrue && ~exist('trueG', 'var')
            trueG = load(trueG_filename);
        end
        figure(1);
        clf;
        subplot(3,1,1);
        if isPlotTrue
            plot(time, trueG.G_model, 'b-', 'LineWidth', 1, 'DisplayName', 'True G');
            hold on;
        end
        plot(output_time, G_model, 'r-', 'LineWidth', 2, 'DisplayName', 'Predicted G');
        legend;
        xlim([min(output_time) max(output_time)]);
        ylim([-0.5 1]);
        text(-20, 0.5, sprintf('Loc = %3.1f, Amp = %3.2f', model.location * prior.dt, model.amplitude),'FontSize', 16);
        xlabel('Time (s)');
        ylabel('Amplitude');
        subplot(3,1,2);
        plot(time, D_model_toplot, 'r--', 'LineWidth', 1, 'DisplayName', 'Predicted D');
        hold on;
        if prior.negOnly
            xline(0, 'r-', 'DisplayName', '');
        end
        hold on;
        plot(time, D, 'k-', 'LineWidth', 2, 'DisplayName', 'Observed D');
        legend;
        xlim([min(time) max(time)]);
        ylim([min(D) max(D)] * 1.1);
        xlabel('Time (s)');
        subplot(3,1,3);
        plot(1:iSaveStep, LikSave(1:iSaveStep), 'k-', 'LineWidth', 2);
        xlim([0 totalSavedSteps]);
        xticklabels(num2str(str2double(xticklabels()) * saveSteps));
        xlabel('Step Count');
        ylabel('Likelihood Probability');
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        iSaveStep = iSaveStep + 1;

    end

end