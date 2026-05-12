function output = PiecewiseArrhenian(FitOrEval, varargin)
% Function "PiecewiseArrhenian.m" v3.0.0, tested 18 September 2024
% Written by T.D. Brown Fall 2023
%
% This function uses many nested functions in order to perform computations
% related to a Piecewise Arrhenian model for temperature-dependent
% electrical conductivity.
%
% Option to either fit data or evaluate a model using boolean FitOrEval
% Function syntax depends on whether fit or evaluation is desired
% For evaluation: econd01 = PiecewiseArrhenian(0, modelpars, temps, order)
% For fitting: modelpars = PiecewiseArrhenian(1, econddata, phaseflag,guess,images)
%
% See documentation for lower functions for more details
%
% Input: FitOrEval = Boolean, (0) for evaluating, (1) for fitting model
% FitOrEval == 0:
% Input: modelpars = struct of model parameters for electrical conductivity
% Input: temps = vector of temperatures [K] at which to evaluate modeel
% Input: order = Boolean, (0) E conductivity, (1) log derivative by T
% Output: output = vector of temperature dependent model output
% FitOrEval == 1:
% Input: econddata = array of data in columns [temperatures, conductivity]
% Input: phaseflag = Boolean, (0) single-phase, (1) two-phase Arrhenian
% Input: guess = vector of guessed temperature thresholds for 2-phase [0-1]
% Input: images = Boolean, (0) suppress plots, (1) show plots

% Use function to evaluate model described by modelpars
if FitOrEval == 0
    modelpars = varargin{1}; tempers = varargin{2}; order = varargin{3};
    output = PiecewiseArrhenianEval(modelpars, tempers, order);
% Use function to fit model to electrical conductivity data
else
    econddata = varargin{1}; phaseflag = varargin{2}; ...
        guess = varargin{3}; images = varargin{4};
    output = PiecewiseArrhenianFitter(econddata, phaseflag, guess, images);
end
end

%% Function Fitting Block

function modelpars = PiecewiseArrhenianFitter(econddata, phaseflag, guess,images)
% Create best-fit Piecewise Arrhenian model to electrical conductivity
% data. Option to create single-phase Arrhenian only. Option to create
% plot of raw data with best fit model
% NewInput: guess = vector [min, max] of guess temperature thresholds
%           thresholds are as a fraction of full-scale temperature [0-1]
% Output: modelpars = struct of model parameters fit to data

    % Option to create single-phase Arrhenian model (semiconductors)
    if phaseflag == 0
        modelpars.lowTpars = ArrhenianFit(econddata);
        modelpars.highTpars = modelpars.lowTpars;
        modelpars.sigpars = [0,0.1];
    
    % Option to create two-phase Piecewise Arrhenian model
    else
        % Use fmincon to find best lowT and highT cutoffs
        temps = econddata(:,1);
        init = [temps(floor(guess(1)*end)), temps(ceil(guess(2)*end))];
        optimoptions('fmincon','OptimalityTolerance',1e-10);

        bestcutoffs = fmincon(@(cutoffs)objfun(cutoffs, econddata),...
            init,[1 -1],0,[0,0],0,...
            [temps(2),temps(2)], [temps(end-2),temps(end-2)]);
    
        % Use best cutoffs to evaluate the rest of the model
        lowTdata = econddata(econddata(:,1)<bestcutoffs(1),:);
        highTdata = econddata(econddata(:,1)>bestcutoffs(2),:);
        modelpars.lowTpars = ArrhenianFit(lowTdata);
        modelpars.highTpars = ArrhenianFit(highTdata);    
        bestsig = ComputeSigmoid(bestcutoffs, econddata);
        modelpars.sigpars = LogisticFit([temps bestsig]);
    end

    % Option to calculate model mse and plot images
    if images 
        modelecond = PiecewiseArrhenianEval(modelpars, econddata(:,1), 0);
        mae = mean( abs(econddata(:,2)-modelecond) );

        figure(1); scatter(1000./econddata(:,1), econddata(:,2), 10, 'ok', 'filled');
        hold on; plot(1000./econddata(:,1), modelecond, '-k');
        set(gcf,'color','w'); set(gca, 'YScale', 'log'); 
        title(['MAE = ',num2str(mae)]);
        xlabel('1000/T (K-1)'); ylabel('E Conductivity (S m-1)');
    end
end

function val = objfun(cutoffs, econddata)
% Computes the loss for determining whether the cutoffs are chosen "well",
% related to properties of the resultant sigmoid and its derivative
%
% NewInput: cutoffs = vector of threshold temperatures for 2-phase model
% Output: val = scalar, loss function for fitting sigmoid
    sigmoid = ComputeSigmoid(cutoffs, econddata);
    
    % Reinterpolate data to assist with derivative approx.
    temps = econddata(:,1);
    moretemps = linspace(temps(1),temps(end),5000);
    sigmoid2 = interp1(temps,sigmoid,moretemps,'spline');
    
    % When sigmoid is "good" these should be proportional
    f1 = sigmoid2.*(1-sigmoid2);
    f2 = numericaldydx(moretemps,sigmoid2);

    % Find best proportionality constant, and evaluate SSE
    %bestscale = sum(abs(sqrt(f1)).*f1./f2) / sum(abs(sqrt(f1)));
    bestscale = sum(f1.*f1./f2) / sum(f1);
    val = sum( (f1-bestscale*f2).^2 );
end

function sigmoid = ComputeSigmoid(cutoffs,econddata)
% Computes the sigmoidal interpolation function for a data set, with the
% lowT and highT Arrhenian fits determined by input cutoffs
%
% Output: sigmoid = vector of computed sigmoid at temperatures in econddata

    % Use cutoffs to specify lowT and highT data
    lowtdata = econddata(econddata(:,1)<=cutoffs(1),:);
    hightdata = econddata(econddata(:,1)>=cutoffs(2),:);

    % Determine best fit Arrhenian for lowT and highT data
    lowtpars = ArrhenianFit(lowtdata);
    hightpars = ArrhenianFit(hightdata);
    lowtlecond = ArrhenianEval(lowtpars,econddata(:,1));
    hightlecond = ArrhenianEval(hightpars,econddata(:,1));

    % Compute sigmoid as rule of mixtures of log conductivities
    sigmoid = (log(econddata(:,2)) - lowtlecond) ./ (hightlecond-lowtlecond);
end

function sigpars = LogisticFit(sigmoiddata)
% Computes the best-fit sigmoid parameters (mu, sigma) to the data
% Logistic: y = (1+exp(-(T-mu)/sigma))^-1
%
% Input: sigmoiddata = array of data in columns [temperatures, sigmoid]
% Output: sigpars = vector of best-fit logistic parameters (mu, sigma) 
    temps = sigmoiddata(:,1);
    sigmoid = sigmoiddata(:,2);
    guess = [mean(temps), 0.1*mean(temps)];
    sigpars = nlinfit(temps,sigmoid, @LogisticEval,guess);
end

function Arrhpars = ArrhenianFit(econddata)
% Computes the best-fit Arrhenian parameters (-Ea/kB, log(sig0)) to the data
% Arrhenian: log(sigma) = -Ea/kB*1/T + log(sig0)
%
% Output: Arrhpars = vector of best-fit Arrhenian parameters (-Ea/kB,log(sig0))
    temps = econddata(:,1);
    lecond = log(econddata(:,2));
    Arrhpars = polyfit(1./temps, lecond,1);
end

%% Function Evaluation Block

function econd01 = PiecewiseArrhenianEval(modelpars, temps, order)
% Evaluates a PiecewiseArrhenian model described by parameters modelpars at
% the specified inputs. Option to return either the function itself (order
% = 0) or its log-derivative for gamma calculations (order=1); 
%
% Output: econd01 = vector of conductivities [S/m] or log-derivatives [K-1]
%                   corresponding to input temps

    % Unpackage modelpars struc
    lowTpars = modelpars.lowTpars;
    highTpars = modelpars.highTpars;
    sigmoidpars = modelpars.sigpars;

    % Evaluate both Arrhenian models at the input temperatures
    lowTlecond = ArrhenianEval(lowTpars,temps);
    highTlecond = ArrhenianEval(highTpars,temps);
    
    % Evaluate the sigmoid model at the input temperatures
    sigmoid = LogisticEval(sigmoidpars,temps);

    % Zeroth-order derivative is simple function evaluation
    if order == 0

        % Multiply according to rule of mixtures and exponentiate
        lecond = (1-sigmoid).*lowTlecond + sigmoid.*highTlecond;
        econd01 = exp(lecond);

    % First-order derivative can be computed analytically
    % This is actually the first-order log-derivative, used for gammas
    elseif order == 1
        
        % Log-derivatives of Arrhenian fits
        lowTdlecond = -lowTpars(1)./temps.^2;
        highTdlecond = -highTpars(1)./temps.^2;

        % First two terms in product rule
        t12 = (1-sigmoid).*lowTdlecond + sigmoid.*highTdlecond;

        % Second two terms in product rule
        dnormT = 1/sigmoidpars(2);
        t34 = dnormT*sigmoid.*(1-sigmoid).*(highTlecond-lowTlecond);
        econd01 = t12 + t34;
    end
end

function sigvals = LogisticEval(sigpars,temps)
% Evaluates a logistic sigmoid model described by parameters sigpars to
% generate the sigmoid values at input temperatures
% 
% NewInput: sigpars = vector of sigmoid parameters (mu, sigma)
% Output: sigvals = vector of sigmoid values computed for input temps
    normT = (temps - sigpars(1)) / sigpars(2);
    sigvals = (1 + exp(-normT)).^-1;
end

function lecond = ArrhenianEval(Arrhpars,temps)
% Evaluates an Arrhenian model described by parameters Arrhpars to generate
% the logarithmic conductivities at input temperatures
%
% NewInput: Arrhpars = vector of Arrhenian parameters (-Ea/kB,log(sig0))
% Output: lecond = vector of log-conductivities corresponding to temps
    % New implementation to improve speed. Notice Arrhpars from polyfit is
    % left right flipped so evaluation is from right to left of Arrhpars
    lecond = Arrhpars(2) + Arrhpars(1)./temps;
end


function dy = numericaldydx(x,y)
% Computes a central difference (finite difference) for vector data
% Uses unbiased estimate, ie (f(x+h)-f(x-h))/2h
%
% Input: x = vector of abscissa values
% Input: y = vector of ordinate values
% Output: dy = vector of estimated dy/dx

    dy = diff(y);dx = diff(x);
    newx = (x([1:numel(x)-1])+x([2:numel(x)]))./2;

    % Calculate derivative at midpoints
    y = (dy(1:numel(dy)-1)+dy(2:numel(dy)))./(dx(1:numel(dx)-1)+dx(2:numel(dx)));
    y1 = y(2) - (y(3)-y(2))./(x(3)-x(2)).*(x(2)-x(1));
    y = flip(y); y(end+1) = y1; y = flip(y);

    % Reinterpolate onto original gridpoints
    y = interp1(newx,y,x,'linear');

    % Linear Extrapolation on endpoints
    y(1) = y(2) - (y(3)-y(2))./(x(3)-x(2)).*(x(2)-x(1));
    y(end) = y(end-1) + (y(end-1)-y(end-2))./(x(end-1)-x(end-2)).*(x(end)-x(end-1));

    dy=y;
end