% main.m

% Script "main.m" v3.0.0, tested 02 Feb 2026
% Originally Written by T.D. Brown Fall 2019 and updated to VC-NDR
% materials by F Jardali Fall 2025

% This script is the main program that calls RK4Nsolver.m and derivfunc.m
% in order to implement Runge-Kutta 4th order solving of the differential
% equation system programmed into derivfunc.m.
%
% It is mostly a convenient compilation of user-defined simulation
% metaparameters, as well as parameters specific to whatever system of
% differential equations is being solved at the time. See RK4Nsolver.m for
% additional details.
%
% The program is currently set-up to perform simulations of steady-state IV
% curves and nonlinear neuronal oscillations for electrothermal neuron
% materials. The program operates in two steps. In the
% first step, the electrothermal models for the given material are used to
% compute the predicted steady state IV curve. In the second step, 
% the prediccted time series of coupled 
% thermal and electrical oscillations are computed, as a function of the 
% user-input circuit capacitance. 
% 
% Input [control]: diffsolve = Boolean; (0) plot IV, (1) simulate oscillations
% Input [control]: lims = [initial_time, final_time], [ROW vector]
% Input [control]: numsteps = number of integrator steps between lims
%
% Input [model]: dims = struct containing device dimensions
% 
% Input [model]: params = struct to pass all external variables to derivfunc
%   params.CelFactor = circuit capacitance as ratio of critical capacitance [1]
%   params.T0 = ambient temperature reservoir [K]
%   params.rhoSS = Limiter series resistance as fraction of R(Tss) [1]
%       *not currently used, set to 0*
%
% Input [model]: material select [Int] = flag to select given material model
% List of flags and material models:
% 1-LCMO
% 
% When diffsolve == 0
% Output: SSTemp = vector of steady state temperatures [K]
% Output: SSCurrent = vector of steady state currents [A]
% Output: SSVolt = vector of steady state voltages [V]
% 
% When diffsolve == 1
% Output: sol = solution array including timesteps and all integrated system variables
% Output: result = vector of oscillation properties, extrema and frequency
%
    close all; clear all
%% User-Modified Block / Parameters
% Boolean variable switches between steady state IV modeling and associated
% figures, or simulating oscillations, Choose =0 or =1. Must run 0 step
% before 1 step in order to determine standard bias point
diffsolve =1; % [Bool] Either plots I-V (0) or runs RK4 integrator (1)

if diffsolve == 0
    clear all;
    diffsolve = 0;
end

% Simulation metaparameters
numsteps = 200000;
lims = [0,5e-3]; % [s], integration domain
initcond = [262, 0.0045]'; % [K, A] initial temperature, current state 

% Define "external" parameters, i.e., those that appear explicitly in
% differential equation system but not in Gel(T) or Gth(T)
params.L = 40e-6; % Series electrical inductance [H]
params.T0 = 100; % Ambient temperature [K]
params.Tss = 262; % Quasisteady state device temperature [K] (=T at bias V)

% Define "internal parameters", those that appear in Gel(T) and Gth(T)
% Internal Parameters: Device dimensions
dims.devtype = 1; % 0 for lateral, 1 for cylindrical
dims.length = 2.5e-6; % Dimension parallel to current flow [m]
dims.width = 5e-6; % Dimension orthogonal to current, parallel to heat [m]
dims.thick = 2.5e-6; % Dimension out-of-plane [m]
dims.radius = 5.0e-6; % Dimension of radius for cylinder

material_select = 1; % [Int] Use to determine which material to model

%% Don't Change Anything Below This!
% Initialize simulations and fit material models

% Determine dimension scaling factors from user input
if dims.devtype == 0
    escale = dims.thick*dims.width / dims.length;
    tscale = dims.width*dims.length / dims.thick;
    tcscale = dims.length*dims.width*dims.thick;
elseif dims.devtype == 1
    escale = pi*dims.radius^2 / dims.thick;
    tscale = 2*pi*dims.thick;
    tcscale = pi*dims.radius^2*dims.thick;
end
funpars.dims = [escale, tscale, tcscale];

% Fit models and repackage in structs for electrothermalmodels.m
% if diffsolve == 0
modelopts = selectMaterial(material_select);

[tcondpars, tcappars] = modelfitter(modelopts.tcondfile,modelopts.tcapfile, ...
    modelopts.tmodelopts,modelopts.tcapmodelopts,modelopts.trange);     

%Constants needed for fitting electrical conductivity data according to
%Alexandrov et al.,PRL 96, 117003 (2006)
econdpars.lowT = [2.8E-3,17E-8,0.05e-12];     %rho02+a*T.^2+b.*T.^4.5;Low temp Ferromagnetic phase
econdpars.highT = [3.44E-6,90.5E-3,8.617e-5]; % rho1= rho0.*T.*exp(Ep./(kB*T));High Temp Paramagnetic phase
econdpars.sigmoid = [267, 8];        % s =(1/2).*erfc((T-TC)/Gamma);

funpars.Epars = econdpars;
funpars.Tpars = tcondpars;
funpars.Tcappars = tcappars;
params.funpars = funpars;

% end
%%
% Initialize RK4 state variables as empty arrays. (Custom RK4 solver is
% capable of integrating with hysteretic state variable, comparing two
% simulations with different timesteps for convergence studies, and
% creating its own figures; none of these features used here).
state.blank = []; prev = []; imgflag = [];

switch diffsolve

 %% diffsolve = 0 performs first step, steady state modeling 
    case 0
     % Calculate steady state IV curve from formulas
     SSTemp = linspace(params.T0,900,5000)';
     ETcond = electrothermalmodels(funpars,SSTemp);
     Econd = ETcond.Econd; Tcond = ETcond.Tcond; Tcap = ETcond.Tcap;

     % Current, voltage as f(temperature) 
     SSCurrent = sqrt((SSTemp-params.T0).*Tcond.*Econd);
     SSVolt = sqrt((SSTemp-params.T0).*Tcond./Econd);

     % Make Figures for steady-state IVT properties
     figure() % Temperature-dependent electrical conductance
     plot(SSTemp, Econd); 
     set(gcf,'color','w'); xlabel('Temperature (K)'); ylabel('E conductance(S)');

     figure() % Temperature-dependent thermal conductance
     plot(SSTemp, Tcond); 
     set(gcf,'color','w'); xlabel('Temperature (K)'); ylabel('T conductance(W/K)');

     figure() % Temperature-dependent thermal capacitance
     plot(SSTemp, Tcap); title('Thermal Capacitance vs. T');
     set(gcf,'color','w'); xlabel('Temperature (K)'); ylabel('T Capacitance (J/K)');
           
      figure() % Steady-state current-voltage curve
      plot(SSVolt,SSCurrent ); title('Steady State I vs V');
      set(gcf,'color','w'); %set(gca,'XScale','log');
      xlabel('SS voltage (V)'); ylabel('SS Current (A)');

      figure() % Steady state temperature-current curve
      plot(SSTemp, SSCurrent); 
      title('Steady State I vs T');
      set(gcf,'color','w'); xlabel('Temperature (K)'); ylabel('SS Current (A)');
case 1
        %% diffsolve = 1 performs second step, oscillation simulations 
       % This block is part of renormalizing time and temperature to help
       % with convergence
% Need to calculate Gth(T0), but this is a constant, so compute it just 
% once, outside of for loops

    outT0 = electrothermalmodels(funpars,params.T0);
    params.GthofT0 = outT0.Tcond;
    params.GelofT0 = outT0.Econd;
    params.CthofT0 = outT0.Tcap;
     
% Also compute at bias point 
    outTss = electrothermalmodels(funpars,params.Tss);
    params.GthofTss = outTss.Tcond;
    params.GelofTss = outTss.Econd;
    params.CthofTss = outTss.Tcap;
 
    params.V0 = sqrt((params.Tss-params.T0).*params.GthofTss./params.GelofTss);
 
    % Normalize by thermal time constant, tth = 2*pi*Cth*Rth(T0)
    thermalt = 2*pi*params.CthofT0/params.GthofT0;


     % Time integrator method with tic / toc
        tic
        
        % This block is part of renormalizing time and temperature to help
        % with convergence of diffeq solver
        lims = lims/thermalt;
        params.Ith = sqrt(params.T0*params.GthofT0*params.GelofT0);
        initcond = [initcond(1)/params.T0; initcond(2)/params.Ith];

        % Function call to RK4Nsolver to implement Runge-Kutta integrator
        [sol,~,err] = RK4Nsolver(lims,1/numsteps,initcond,prev,imgflag,[],params);        
        toc

        % Overwrite diffeq solution with rescaled values
        sol(1,:) = thermalt*sol(1,:);
        sol(2,:) = params.T0*sol(2,:);
        sol(3,:) = params.Ith*sol(3,:);
%% Make Figures
        figure()% Temp oscillations
        plot(sol(1,:),sol(2,:),'Color', [0.8125 0 0] );
        set(gca,'XScale','linear')
        set(gca,'YScale','log')
        ylabel('T (K)');
        xlabel('Time (s)');
        title('Temperature Trajectory')
   
        figure()% Current oscillations
        plot(sol(1,:),sol(3,:))
        set(gca,'XScale','linear')
        ylabel('Current (A)')
        xlabel('Time (s)')
        title('Current Trajectory')

% Compute minT, maxT, minI, maxI, and freq
oscResults = oscAnalyzer(sol);

% Calculate power [W] at bias point 
P = params.GthofTss*(params.Tss-params.T0);

% Calculate deltaT and deltaI
deltaT = oscResults(2)-oscResults(1);
deltaI = oscResults(4)-oscResults(3);

%  Cap, minT, maxT, minI, maxI, freq, power, deltaT, deltaI to results 
result = [params.L, oscResults, P,deltaT,deltaI];
end
%% 

function result = oscAnalyzer(in)

    % Remove all but the last 10% of the data
    cutoff = round(0.9*size(in,2));
    dat = in(:,cutoff:size(in,2));

    % Maxima analysis
    [pks,locs] = findpeaks(dat(2,:),dat(1,:));
    perd = mean(diff(locs));
    freq = 1/perd;
    tmax = mean(pks);
    pks = findpeaks(dat(3,:),dat(1,:));
    vmax = mean(pks);

    % Minima analysis
    dat2 = dat;
    dat2(2,:) = 1.5*max(dat(2,:)) - dat2(2,:);
    dat2(3,:) = 1.5*max(dat(3,:)) - dat2(3,:);

    pks = findpeaks(dat2(2,:),dat2(1,:));
    tmin = 1.5*max(dat(2,:))-mean(pks);
    pks = findpeaks(dat2(3,:),dat2(1,:));
    vmin = 1.5*max(dat(3,:))-mean(pks);

    % Give results
    result = [tmin, tmax, vmin, vmax, freq];

end

function modelopts = selectMaterial(material_select)
% Function "selectMaterial.m" v2.0.0, tested 18 September 2024

    switch material_select         
     
        case 1 % LCMO
        trange = linspace(10,900, 5000); name = 'CMR';
        tcondfile = 'ThermCond.xlsx'; 
        econdfile = ''; 
        tcapfile = 'HeatCapa.xlsx';
%         emodelopts = [1,0.6,0.8]; 
        tmodelopts = [1, 0.25, 0.33]; 
        tcapmodelopts = [1,0.1,0.4];

    end

    % Wrap up output in struct modelopts
    modelopts.trange = trange; 
%     modelopts.econdfile = econdfile; 
    modelopts.tcondfile = tcondfile; 
    modelopts.tcapfile = tcapfile;
%     modelopts.emodelopts = emodelopts; 
    modelopts.tmodelopts = tmodelopts; 
    modelopts.tcapmodelopts = tcapmodelopts; 
    modelopts.name = name;


end