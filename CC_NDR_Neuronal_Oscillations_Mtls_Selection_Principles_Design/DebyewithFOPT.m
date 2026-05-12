function output = DebyewithFOPT(FitOrEval, varargin)
% Function "DebyewithFOPT.m" v3.0.0, tested 18 September 2024
% Written by T.D. Brown Fall 2023
%
% This function uses many nested functions in order to perform computations
% related to a Piecewise Debye model for temperature-dependent
% thermal capacitance (specific heat capacity)
%
% Option to either fit data or evaluate a model using FitOrEval
% Function syntax depends on whether fit or evaluation is desired
% For evaluation: econd01 = DebyewithFOPT(0, modelpars, temps, order)
% For fitting: modelpars = DebyewithFOPT(1, econddata, phaseflag,images)
%
% See documentation for lower functions for more details
%
% Input: FitOrEval = Boolean, (0) for evaluating, (1) for fitting model
% FitOrEval == 0:
% Input: modelpars = struct of model parameters for heat capacity
% Input: temps = vector of temperatures [K] at which to evaluate model
% Output: output = vector of temperature dependent specific heat [J/m3K]
% FitOrEval == 1:
% Input: tcapdata = array of data in columns [temperatures, specific heat]
% Input: phaseflag = Boolean, (0) single-phase, (1) two-phase Debye
% Input: guess = vector of guessed temperature thresholds for 2-phase [0-1]
% Input: images = Boolean, (0) suppress plots, (1) show plots

% Use function to evaluate model described by modelpars
if FitOrEval == 0
    modelpars = varargin{1}; tempers = varargin{2};
    output = DebyewithFOPTEval(modelpars, tempers);
% Use function to fit model to heat capacity data
else
    tcapdata = varargin{1}; phaseflag = varargin{2}; ...
        guess = varargin{3}; images = varargin{4};
    output = DebyewithFOPTFitter(tcapdata, phaseflag, guess, images);
end
end

%% Function Fitting Block

function modelpars = DebyewithFOPTFitter(tcapdata, phaseflag, guess, images)
% Create best-fit Piecewise Debye model to specific heat capacity data. 
% Option to create single-phase Debye only. Option to create
% plot of raw data with best fit model
% NewInput: guess = vector [min, max] of guess temperature thresholds
%           thresholds are as a fraction of full-scale temperature [0-1]
% Output: modelpars = vector of model parameters fit to data

    % Option to create Debye model without FOPT
    if phaseflag == 0
        modelpars.lowTpars = DebyeFit(tcapdata);
        modelpars.highTpars = modelpars.lowTpars;
        modelpars.FOPTpars = [0,0,0];
    
    % Option to create Debye model with peak for FOPT
    else
        
        % Use fmincon to find best lowT and highT cutoffs
        temps = tcapdata(:,1);
        init = [temps(floor(guess(1)*end)), temps(ceil(guess(2)*end))];
        opts = optimoptions(@fmincon,'StepTolerance',0.05);
        warning('off','all');

        bestcutoffs = fmincon(@(cutoffs)objfun(cutoffs, tcapdata),...
            init,[1 -1],0,[0,0],0,...
            [temps(2),temps(2)], [temps(end-2),temps(end-2)],...
            [],opts);
    
        % Use best cutoffs to evaluate the rest of the model
        lowTdata = tcapdata(tcapdata(:,1)<bestcutoffs(1),:);
        highTdata = tcapdata(tcapdata(:,1)>bestcutoffs(2),:);
        modelpars.lowTpars = DebyeFit(lowTdata);
        modelpars.highTpars = DebyeFit(highTdata);
        
        % Fit first order phase transition heat capacity peak
        bestgauss = tcapdata(:,2) - ComputeFullDebye(bestcutoffs,tcapdata);
        modelpars.FOPTpars = GaussFit([temps, bestgauss]);    
    end

    % Option to calculate model mse and plot images
    if images 
        modeltcap = DebyewithFOPTEval(modelpars, tcapdata(:,1));
        mae = mean( abs(tcapdata(:,2)-modeltcap) );

        figure(1); scatter(tcapdata(:,1), tcapdata(:,2), 10, 'ok', 'filled');
        hold on; plot(tcapdata(:,1), modeltcap, '-k');
        set(gcf,'color','w'); ylim([0,max(tcapdata(:,2))])
        title(['MAE = ',num2str(mae)]);
        xlabel('Temperature (K)'); ylabel('Heat Capacity (J m-3 K-1)');
    end
end

function val = objfun(cutoffs, tcapdata)
% Computes the loss for determining whether the cutoffs are chosen "well",
% using a simple overall MSE loss
%
% NewInput: cutoffs = vector of threshold temperatures for 2-phase model
% Output: val = scalar, loss function for fitting sigmoid
    FDebye = ComputeFullDebye(cutoffs, tcapdata);
    temps = tcapdata(:,1);

    FDebye(or( temps<cutoffs(1), temps>cutoffs(2)  )) = [];
    tcapdata(or( temps<cutoffs(1), temps>cutoffs(2)  ),:) = []; 

    % Objective function is squared sum of errors
    val = sum((FDebye-tcapdata(:,2)).^2);
end

function fullDebye = ComputeFullDebye(cutoffs,tcapdata)
% Computes the best-fit 2-phase Debye background assuming single phase
% thresholds at cutoffs for data in tcapdata
%
% Output: fullDebye = vector of computed specific heat as 2-phase Debye

    % Use cutoffs to specify lowT and highT data for Debye fit
    lowtdata = tcapdata(tcapdata(:,1)<=cutoffs(1),:);
    hightdata = tcapdata(tcapdata(:,1)>=cutoffs(2),:);
    Tdiff = cutoffs(2) - cutoffs(1);
    Tmean = 0.5*(cutoffs(1) + cutoffs(2));

    % Determine best fit Debye model through lowT and highT data
    lowTpars = DebyeFit(lowtdata);
    highTpars = DebyeFit(hightdata);

    % Interpolate between low and high T heat capacity
    sig = (1 + exp(-5/Tdiff*(tcapdata(:,1) - Tmean))).^-1;
    fullDebye = (1-sig).*DebyeEval(lowTpars,tcapdata(:,1)) + sig.*DebyeEval(highTpars,tcapdata(:,1));
end

function gausspars = GaussFit(gaussdata)
% Computes the best-fit sigmoid parameters (height, mu, sigma) to the data
% Gauss: y = height*exp( (-(T-mu)/sigma))^2)
%
% Input: gaussdata = array of data in columns [temperatures, gauss]
% Output: gausspars = vector of best-fit Gauss parameters (height, mu, sigma) 
    temps = gaussdata(:,1);
    gauss = gaussdata(:,2);
    guess = [max(gauss),mean(temps), 0.1*mean(temps)];
    gausspars = nlinfit(temps,gauss, @GaussianEval,guess);
end

function debyepars = DebyeFit(tcapdata)
% Computes the best-fit Debye model parameters (molvol, T_D) to the data in
% input array tcapdata
%
% Output: debyepars = vector of bestfit Debye parameters (molvol, T_D) 
    temps = tcapdata(:,1);
    debyes = tcapdata(:,2);
    guess = [mean(debyes)/(8.314), 3*max(temps)];
    debyepars = nlinfit(temps,debyes, @DebyeEval,guess);
end

%% Function Evaluation Block

function tcap0 = DebyewithFOPTEval(modelpars, temps)
% Evaluates a DebyewithFOPT model described by parameters modelpars at
% the specified inputs. Only returns function, not derivative 
%
% Output: tcap0 = vector of specific heats [J/m3K] corresponding to input temps

    % Unpackage modelpars struc
    lowTpars = modelpars.lowTpars;
    highTpars = modelpars.highTpars;
    FOPTpars = modelpars.FOPTpars;

    % Evaluate just the Debye part of the model at the input temperatures   
    % Evaluate the FOPT peak part of the model at the input temperatures
    if FOPTpars(1) ~= 0
        lowTtcap = DebyeEval(lowTpars,temps);
        highTtcap = DebyeEval(highTpars,temps);
        sig = (1+exp(-1/abs(FOPTpars(3))*(temps-abs(FOPTpars(2))))).^-1;
        backgroundtcap = (1-sig).*lowTtcap + sig.*highTtcap;
        FOPTtcap = GaussianEval(FOPTpars,temps);
    elseif FOPTpars(1) == 0
        backgroundtcap = DebyeEval(lowTpars,temps);
        FOPTtcap = zeros(size(temps));
    end

    % Final model is the sum of the Debye background, plus the FOPT peak
    tcap0 = backgroundtcap + FOPTtcap;

    % No need to compute first derivative, e.g., for gamma parameters
end

function gaussvals = GaussianEval(gausspars,temps)
% Evaluates a Gaussian model described by parameters gausspars to
% generate the gauss values at input temperatures
% 
% NewInput: gausspars = vector of Gauss parameters (height, mu, sigma)
% Output: gaussvals = vector of Gauss values computed for input temps
    normT = (temps - gausspars(2)) / gausspars(3);
    gaussvals = abs(gausspars(1))*exp(-0.5*normT.^2);
end

function debyevals = DebyeEval(debyepars,temps)
% Evaluates a Debye model described by parameters debyepars to generate
% the background Debye Cth values at input temperatures
% 
% NewInput: debyepars = vector of Debye parameters (molvol, T_D)
% Output: debyevals = vector of Debye values computed for input temps
    
    pref = debyepars(1); % Reciprocal Molar volume [mol / m3]
    theta = max(10,abs(debyepars(2)) ); % Debye temperature TD [K]

    normT = theta./temps; % TD/T
    
    %Slowww naive implementation 
    %{ 
    debyevals = pref*9*8.3145./normT.^3.*...
        arrayfun(@(y)DebyeInt(y),normT); % Vectorized
    %}
    
    normd = zeros(size(normT));
    for ii = 1:numel(normT)
        normd(ii) = DebyeInt(normT(ii));
    end
    debyevals = pref*9*8.3145./normT.^3.*normd;
end

function out = DebyeInt(in)
    N = 200;
    x = linspace(1e-6,in,N);
    out = sum(kernf(x)).*(in/(N-1));
end

function out = kernf(in)
% Evaluates the kernel function (integrand) of the Debye function, ie 
% f = x^4*exp(-x)/(1-exp(-x))^2
    u = exp(-in);
    v = in.^2 ./ (1-u);
    out = v.^2.*u;
end