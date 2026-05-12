function output = PiecewiseLaurent(FitOrEval, varargin)
% Function "PiecewiseLaurent.m" v3.0.0, tested 18 September 2024
% Written by T.D. Brown Fall 2023
%
% This function uses many nested functions in order to perform computations
% related to a Piecewise Arrhenian model for temperature-dependent
% electrical conductivity
%
% Option to either fit data or evaluate a model using boolean FitOrEval
% Function syntax depends on whether fit or evaluation is desired
% For evaluation: econd01 = PiecewiseLaurent(0, modelpars, temps, order)
% For fitting: modelpars = PiecewiseLaurent(1, econddata, phaseflag,images)
%
% See documentation for lower functions for more details
%
% Input: FitOrEval = Boolean, (0) for evaluating, (1) for fitting model
% FitOrEval == 0:
% Input: modelpars = struct of model parameters for thermal conductivity
% Input: temps = vector of temperatures [K] at which to evaluate modeel
% Input: order = Boolean, (0) T conductivity, (1) log derivative by T
% Output: output = vector of temperature dependent model output
% FitOrEval == 1:
% Input: tconddata = array of data in columns [temperatures, conductivity]
% Input: phaseflag = Boolean, (0) single-phase, (1) two-phase Laurent
% Input: guess = vector of guessed temperature thresholds for 2-phase [0-1]
% Input: images = Boolean, (0) suppress plots, (1) show plots

% Use function to evaluate model described by modelpars
if FitOrEval == 0
    modelpars = varargin{1}; tempers = varargin{2}; order = varargin{3};
    output = PiecewiseLaurentEval(modelpars, tempers, order);
% Use function to fit model to electrical conductivity data
else
    tconddata = varargin{1}; phaseflag = varargin{2}; ...
        guess = varargin{3}; images = varargin{4};
    output = PiecewiseLaurentFitter(tconddata, phaseflag, guess, images);
end
end

%% Function Fitting Block

function modelpars = PiecewiseLaurentFitter(tconddata, phaseflag, guess,images)
% Create best-fit Piecewise Laurent model to thermal conductivity
% data. Option to create single-phase Laurent only. Option to create
% plot of raw data with best fit model
% NewInput: guess = vector [min, max] of guess temperature thresholds
%           thresholds are as a fraction of full-scale temperature [0-1]
% Output: modelpars = struct of model parameters fit to data

    % Option to create single-phase Laurent model 
    if phaseflag == 0
        modelpars.lowTpars = LaurentFit(tconddata);
        modelpars.highTpars = modelpars.lowTpars;
        modelpars.sigpars = [0,0.1];
    
    % Option to create two-phase Piecewise Laurent model
    else
        % Use fmincon to find best lowT and highT cutoffs
        temps = tconddata(:,1);
        init = [temps(floor(guess(1)*end)), temps(ceil(guess(2)*end))];
        optimoptions('fmincon','OptimalityTolerance',1e-10);

        bestcutoffs = fmincon(@(cutoffs)objfun(cutoffs, tconddata),...
            init,[1 -1],0,[0,0],0,...
            [temps(2),temps(2)], [temps(end-2),temps(end-2)]);
    
        % Use best cutoffs to evaluate the rest of the model
        lowTdata = tconddata(tconddata(:,1)<bestcutoffs(1),:);
        highTdata = tconddata(tconddata(:,1)>bestcutoffs(2),:);
        modelpars.lowTpars = LaurentFit(lowTdata);
        modelpars.highTpars = LaurentFit(highTdata);    
        bestsig = ComputeSigmoid(bestcutoffs, tconddata);
        modelpars.sigpars = LogisticFit([temps bestsig]);

    end
    modelpars.endpoints = [min(tconddata(:,2)),max(tconddata(:,2))]; %!!!

    % Option to calculate model mse and plot images
    if images 
        modeltcond = PiecewiseLaurentEval(modelpars, tconddata(:,1), 0);
        lowTmodel = LaurentEval(modelpars.lowTpars, tconddata(:,1));
        highTmodel = LaurentEval(modelpars.highTpars, tconddata(:,1)); 
        mae = mean( abs(tconddata(:,2)-modeltcond) );

        figure(1); scatter(tconddata(:,1), tconddata(:,2), 10, 'ok', 'filled');
        hold on; plot(tconddata(:,1), modeltcond, '-k');
        set(gcf,'color','w');
        title(['MAE = ',num2str(mae)]);
        xlabel('Temperature (K)'); ylabel('T Conductivity (W m-1K-1)');
    end
end

function val = objfun(cutoffs, tconddata)
% Computes the loss for determining whether the cutoffs are chosen "well",
% related to properties of the resultant sigmoid and its derivative
%
% NewInput: cutoffs = vector of threshold temperatures for 2-phase model
% Output: val = scalar, loss function for fitting sigmoid
    sigmoid = ComputeSigmoid(cutoffs, tconddata);
    
    % Reinterpolate data to assist with derivative approx.
    temps = tconddata(:,1);
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

function sigmoid = ComputeSigmoid(cutoffs,tconddata)
% Computes the sigmoidal interpolation function for a data set, with the
% lowT and highT Laurent fits determined by input cutoffs
%
% Output: sigmoid = vector of computed sigmoid at temperatures in tconddata

    % Use cutoffs to specify lowT and highT data
    lowtdata = tconddata(tconddata(:,1)<=cutoffs(1),:);
    hightdata = tconddata(tconddata(:,1)>=cutoffs(2),:);

    % Determine best fit Laurent for lowT and highT data
    lowtpars = LaurentFit(lowtdata);
    hightpars = LaurentFit(hightdata);
    lowttcond = LaurentEval(lowtpars,tconddata(:,1));
    highttcond = LaurentEval(hightpars,tconddata(:,1));

    % Compute sigmoid as rule of mixtures of log conductivities
    sigmoid = (tconddata(:,2) - lowttcond) ./ (highttcond-lowttcond);
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
    sigpars(1) = max(sigpars(1),temps(2));
    sigpars(2) = max(sigpars(2),5);
end

function Laurentpars = LaurentFit(tconddata)
% Computes the best-fit third-order Laurent parameters 
% (b{-2}, b{-1}, b{0}, b{1}) to the data
% Laurent: kappa = b{-2}*T^-2 + b{-1}*T^-1 + b{0} + b{1}*T^1
%
% Output: Laurentpars = vector of best-fit Laurent parameters (b{-2}...b{1})
    temps = tconddata(:,1);
    tcond = tconddata(:,2);
    inarr = [temps.^-2, temps.^-1, ones(numel(temps),1),temps];

    % If positive slope, only use positive exponents
    if tcond(1)<=tcond(end)
        inarr(:,1:2) = [];
    % If negative slope, only use negative exponents
    else
        inarr(:,4) = [];
    end

    bb = (inarr'*inarr)\eye(size(inarr'*inarr));
    B_mat = bb*inarr';
    Laurentpars = B_mat*tcond;

    if tcond(1)<=tcond(end)
        Laurentpars = [0; 0; Laurentpars];
    else
        Laurentpars = [Laurentpars; 0;];
    end

end

%% Function Evaluation Block

function tcond01 = PiecewiseLaurentEval(modelpars, temps, order)
% Evaluates a PiecewiseLaurent model described by parameters modelpars at
% the specified inputs. Option to return either the function itself (order
% = 0) or its log-derivative for gamma calculations (order=1);
%
% Output: econd01 = vector of conductivities [W/mK] or log-derivatives [K-1]
%                   corresponding to input temps

    % Unpackage modelpars struc
    lowTpars = modelpars.lowTpars;
    highTpars = modelpars.highTpars;
    sigmoidpars = modelpars.sigpars;

    endpoints = modelpars.endpoints; %!!!

    % Evaluate both Arrhenian models at the input temperatures
    lowTtcond = LaurentEval(lowTpars,temps);
    highTtcond = LaurentEval(highTpars,temps);
    
    % Evaluate the sigmoid model at the input temperatures
    sigmoid = LogisticEval(sigmoidpars,temps);

    % Zeroth-order derivative is simple function evaluation
    if order == 0

        % Multiply according to rule of mixtures and exponentiate
        tcond01 = (1-sigmoid).*lowTtcond + sigmoid.*highTtcond;
        tcond01 = min(tcond01, 10*endpoints(2)); %!!! Extrapolation bounds
        tcond01 = max(tcond01, 0.1*endpoints(1)); %!!! Extrapolation bounds

    % First-order derivative can be computed analytically
    % This is actually the first-order log-derivative, used for gammas
    elseif order == 1
        
        % Log-derivatives of Laurent fits
        lowTdmat = [-2*lowTpars(1), -lowTpars(2), ...
            0, lowTpars(4),];
        lowTdtcond = [temps.^-3, temps.^-2, zeros(size(temps)), ...
            ones(size(temps))]*lowTdmat';
        lowTtcond = LaurentEval(lowTpars,temps);

        highTdmat = [-2*highTpars(1), -highTpars(2), ...
            0, highTpars(4),];
        highTdtcond = [temps.^-3, temps.^-2, zeros(size(temps)), ...
            ones(size(temps))]*highTdmat';
        highTtcond = LaurentEval(highTpars,temps);

        tcond = (1-sigmoid).*lowTtcond + sigmoid.*highTtcond;

        % First two terms in product rule
        t12 = (1-sigmoid).*lowTdtcond + sigmoid.*highTdtcond;

        % Second two terms in product rule
        dnormT = 1/sigmoidpars(2);
        t34 = dnormT*sigmoid.*(1-sigmoid).*(highTtcond-lowTtcond);
        tcond01 = (t12 + t34)./tcond;
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

function tcond = LaurentEval(Laurentpars,temps)
% Evaluates a Laurent model described by parameters Laurentpars to generate
% the thermal conductivities at input temperatures
%
% NewInput: Laurentpars = vector of Laurent parameters (b{-2},...,b{1})
% Output: tcond = vector of T conductivities corresponding to temps
    
    % A lot of heartache with temps being a row vector! Avoid
    tarr = [temps.^-2, temps.^-1, ones(size(temps)),temps];
    tcond = tarr*Laurentpars;
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
