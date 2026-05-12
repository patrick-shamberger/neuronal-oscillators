function [tcondpars,tcappars] = modelfitter(tcondfile,tcapfile,tmodelopts,tcapmodelopts,trange)
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


    close all

    %% First load data and fit models to it
    tconddata = xlsread(tcondfile);
    tcapdata = xlsread(tcapfile);

    % Make sure trange is a column vector
    if size(trange,1) < size(trange,2)
        trange = trange';
    end
    
    %% 
   
    % Fit PiecewiseLaurent model to thermal conductivity data
    fprintf('Fitting Thermal Conductivity Model\n')


        tcondpars = PiecewiseLaurent(1, tconddata, tmodelopts(1), tmodelopts(2:3), 0);
        tcondmodel = [trange, PiecewiseLaurent(0, tcondpars, trange, 0)];
    
        figure(1) % Thermal conductivity + fit 
        scatter(tconddata(:,1),tconddata(:,2),10,'ok','filled');
        hold on
        plot(tcondmodel(:,1), tcondmodel(:,2),'-k');
        set(gcf,'color','w'); title('T Conductivity');
        xlabel('Temperature (K)'); ylabel('Thermal conductivity (W m-1K-1)');
        hold off

    % Fit DebyewithFOPT model to specific heat data
    fprintf('Fitting Thermal Capacity Model\n')

        tcappars = DebyewithFOPT(1, tcapdata, tcapmodelopts(1), tcapmodelopts(2:3), 0);
        tcapmodel = [trange, DebyewithFOPT(0, tcappars, trange)];
    
        figure(2) % Thermal conductivity + fit 
        scatter(tcapdata(:,1),tcapdata(:,2),10,'ok','filled');
        hold on
        plot(tcapmodel(:,1), tcapmodel(:,2),'-k');
        set(gcf,'color','w'); title('Specific heat');
        ylim([0,max( max(tcapmodel(:,2)), max(tcapdata(:,2)) )]);
        xlabel('Temperature (K)'); ylabel('Specific heat (J m-3K-1)');
        hold off

  
 end