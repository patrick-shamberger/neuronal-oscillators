function [econdpars, tcondpars, tcappars] = modelfitter(econdfile, tcondfile, tcapfile, emodelopts, tmodelopts, tcapmodelopts, trange)
% Function "modelfitter.m" v3.0.0, tested 18 September 2024
% Written by T.D. Brown Fall 2023
%
% This function loads user data stored in excel files and fits appropriate
% models to each data set. Data and models are: PiecewiseArrhenian fit to
% T-dependent electrical conductivity, PiecewiseLaurent fit to T-dependent
% thermal conductivity, DebyewithFOPT fit to T-dependent specific heat.
%
% Option to keep trying different guesses for model parameters until user
% is satisfied with fit model.
%
% See "modelfitter.m", "PiecewiseArrhenian.m", "PiecewiseLaurent.m",
% and "DebyewithFOPT.m"
%
% Input: econdfile = string, name of *.xlsx file with electrical conductivity data
% Input: tcondfile = string, name of *.xlsx file with thermal conductivity data
% Input: tcapfile = string, name of *.xlsx file with specific heat data
%   Data in columns [temperatures, electrothermal data] 
% Input: emodelopts = vector of fitting options for electrical conductivity
% Input: tmodelopts = vector of fitting options for thermal conductivity
% Input: tcapmodelopts = vector of fitting options for specific heat
%   Variables emodelopts and tmodel need to be in a specific format
%   They are all a 1x3 vector
%   The first element is 0 or 1, depending on if you want a single-phase or
%   two-phase mode lwith sigmoid interpolant
%   The second two elements only apply for a sigmoidal 2-phase model;
%   they give an initial guess for what the temperature cutoffs should be
%   Example: emodel = [1,0.4,0.6] fits a piecewise Arrhenian model, with
%   the initial guess being at 40% and 60% of the temperature range
% Input: trange = vector of temperatures [K] over which to fit data
%
% Output: econdpars = struct of model parameters for electrical
%   conductivity. Right syntax for "PiecewiseArrhenian.m"
% Output: tcondpars = struct of model parameters for thermal
%   conductivity. Right syntax for "PiecewiseLaurent.m"
% Output: tcappars = struct of model parameters for specific heat.
%   Right syntax for "DebyewithFOPT.m"

    close all

    %% First load data and fit models to it
    econddata = xlsread(econdfile);
    tconddata = xlsread(tcondfile);
    tcapdata = xlsread(tcapfile);

    % Make sure trange is a column vector
    if size(trange,1) < size(trange,2)
        trange = trange';
    end
    
    %% Loop until model parameters are well fit
    
    % Fit PiecewiseArrhenian model to electrical conductivity data
    fprintf('Fitting Electrical Conductivity Model\n')

    for ii = 1:50
        % Attempt Fit
        econdpars = PiecewiseArrhenian(1, econddata, emodelopts(1), emodelopts(2:3), 0);
        econdmodel = [trange, PiecewiseArrhenian(0, econdpars, trange, 0)];
    
        % Plot data and model together
        figure(1) % Electrical conductivity + fit 
        scatter(1000./econddata(:,1),econddata(:,2),10,'ok','filled');
        hold on
        plot(1000./econdmodel(:,1), econdmodel(:,2),'-k');
        set(gcf,'color','w'); set(gca,'YScale','log'); title('E Conductivity');
        xlabel('1000/T'); ylabel('Electrical conductivity (S m-1)');
        hold off

        % Loop if refitting is desired
        reply = input('Re-fit Model? (Enter ''Y''/''N''):\n');

        if reply == 'Y'
            emodelopts = input(['Enter new vector of options for fitting model:\n' ...
            'First element: Use sigmoid (1) or no (0)\n' ...
            'Second element: If sigmoid, guess for first cutoff T (%)\n' ...
            'Third element: If sigmoid, guess for second cutoff T (%)\n']);
        else
            break
        end
    end

    % Fit PiecewiseLaurent model to thermal conductivity data
    fprintf('Fitting Thermal Conductivity Model\n')

    for jj = 1:50
        tcondpars = PiecewiseLaurent(1, tconddata, tmodelopts(1), tmodelopts(2:3), 0);
        tcondmodel = [trange, PiecewiseLaurent(0, tcondpars, trange, 0)];
    
        figure(1) % Thermal conductivity + fit 
        scatter(tconddata(:,1),tconddata(:,2),10,'ok','filled');
        hold on
        plot(tcondmodel(:,1), tcondmodel(:,2),'-k');
        set(gcf,'color','w'); title('T Conductivity');
        xlabel('Temperature (K)'); ylabel('Thermal conductivity (W m-1K-1)');
        hold off

        reply = input('Re-fit Model? (Enter ''Y''/''N''):\n');

        if reply == 'Y'
            tmodelopts = input(['Enter new vector of options for fitting model:\n' ...
            'First element: Use sigmoid (1) or no (0)\n' ...
            'Second element: If sigmoid, guess for first cutoff T (%)\n' ...
            'Third element: If sigmoid, guess for second cutoff T (%)\n']);
        else
            break
        end
    end

    % Fit DebyewithFOPT model to specific heat data
    fprintf('Fitting Thermal Capacity Model\n')

    for kk = 1:50
        tcappars = DebyewithFOPT(1, tcapdata, tcapmodelopts(1), tcapmodelopts(2:3), 0);
        tcapmodel = [trange, DebyewithFOPT(0, tcappars, trange)];
    
        figure(1) % Thermal conductivity + fit 
        scatter(tcapdata(:,1),tcapdata(:,2),10,'ok','filled');
        hold on
        plot(tcapmodel(:,1), tcapmodel(:,2),'-k');
        set(gcf,'color','w'); title('T Capacitance');
        ylim([0,max( max(tcapmodel(:,2)), max(tcapdata(:,2)) )]);
        xlabel('Temperature (K)'); ylabel('Thermal capacitance (J m-3K-1)');
        hold off

        reply = input('Re-fit Model? (Enter ''Y''/''N''):\n');

        if reply == 'Y'
            tcapmodelopts = input(['Enter new vector of options for fitting model:\n' ...
            'First element: Use FOPT peak (1) or no (0)\n' ...
            'Second element: If FOPT, guess for first cutoff T (%)\n' ...
            'Third element: If FOPT, guess for second cutoff T (%)\n']);
        else
            break
        end
    end

    fprintf('If re-fit, overwrite initial emodelopts, tmodelopts, and tcapopts with these values\n');
    fprintf(['Emodelopts = ',num2str(emodelopts),'\n']);
    fprintf(['Tmodelopts = ',num2str(tmodelopts),'\n']);
    fprintf(['Tcapmodelopts = ',num2str(tcapmodelopts),'\n']);

 end