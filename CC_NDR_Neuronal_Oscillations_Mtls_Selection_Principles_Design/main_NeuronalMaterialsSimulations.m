% main_NeuronalMaterialsSimulations.m

% Script "main_NeuronalMaterialsSimulations.m" v3.0.0, tested 29 June 2021
% Written by T.D. Brown + J. Chong Spring 2024
%
% THIS PROGRAM DEPENDS ON RK4NSOLVER.M, DERIVFUNC.M!!
% 
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
% materials, as in the main text. The program operates in two steps. In the
% first step, the electrothermal models for the given material are used to
% compute the predicted steady state IV curve, and the standardized bias
% point with maximum dV/dT is found, along with the critical capacitance at
% that point. In the second step, the prediccted time series of coupled 
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
%   params.T0 = ambient temperature reservoir, usually 300 K
%   params.rhoSS = Limiter series resistance as fraction of R(Tss) [1]
%       *not currently used, set to 0*
%
% Input [model]: material select [Int] = flag to select given material model
% List of flags and material models:
% 1-AlN, 2-Ga2O3(beta), 3-Bi2Te3, 4-CdS, 5-CdSe, 6-Co3O4, 7-CuO, 8-Fe3O4,
% 9-FeO, 10-GaAs, 11-GaN, 12-GaP, 13-GaSb, 14-GaSe, 15-Ge, 16-InSb,
% 17-InSe, 18-LaCoO3, 19-MoO3, 20-MoS2, 21-MoSe2, 22-NbO2, 23-NdNiO3, 24-PbTe,
% 25-Sb2Te3, 26-Si, 27-SiC, 28-SnO2, 29-SnSe, 30-SrTiO3, 31-TaS2, 32-TbMnO3,
% 33-TiO2, 34-V3O5, 35-VO2, 36-WO3, 37-WSe2, 38-ZnO, 39-ZnSe, 40-ZnTe
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

close all; 

%% User-Modified Block / Parameters

% Boolean variable switches between steady state IV modeling and associated
% figures, or simulating oscillations, Choose =0 or =1. Must run 0 step
% before 1 step in order to determine standard bias point
diffsolve =0; % [Bool] Either plots I-V (0) or runs RK4 integrator (1)

if diffsolve == 0
    clear all;
    diffsolve = 0;
end

% Simulation metaparameters
numsteps = 400000; % [Int] Number of timesteps for RK4 integrator
lims = [0,300E-6]; % [s], integration domain

% Define "external" parameters, i.e., those that appear explicitly in
% differential equation system but not in Gel(T) or Gth(T)
params.CelFactor = 1.002; % Electrical capacitance multiplication factor
params.T0 = 300; % Ambient temperature [K]
params.rhoSS = 0.00; % What fraction of total resistance is Rlim? [0-1]

% Define "internal parameters", those that appear in Gel(T) and Gth(T)
% Internal Parameters: Device dimensions
dims.devtype = 1; % 0 for lateral, 1 for cylindrical
dims.length = 1e-6; % Dimension parallel to current flow [m]
dims.width = 1e-6; % Dimension orthogonal to current, parallel to heat [m]
dims.thick = 1e-7; % Dimension out-of-plane [m]
dims.radius = 1e-6; % Dimension of radius for cylinder

% Internal Parameters: Intrinsic material properties
% List of flags and material models:
% 1-AlN, 2-Ga2O3(beta), 3-Bi2Te3, 4-CdS, 5-CdSe, 6-Co3O4, 7-CuO, 8-Fe3O4,
% 9-FeO, 10-GaAs, 11-GaN, 12-GaP, 13-GaSb, 14-GaSe, 15-Ge, 16-InSb,
% 17-InSe, 18-LaCoO3, 19-MoO3, 20-MoS2, 21-MoSe2, 22-NbO2, 23-NdNiO3, 24-PbTe,
% 25-Sb2Te3, 26-Si, 27-SiC, 28-SnO2, 29-SnSe, 30-SrTiO3, 31-TaS2, 32-TbMnO3,
% 33-TiO2, 34-V3O5, 35-VO2, 36-WO3, 37-WSe2, 38-ZnO, 39-ZnSe, 40-ZnTe
material_select = 18; % [Int] Use to determine which material to model

%% Don't Change Anything Below This!

%% Initialize simulations and fit material models

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
if diffsolve == 0
modelopts = selectMaterial(material_select);
[econdpars, tcondpars, tcappars] = modelfitter(modelopts.econdfile, ...
    modelopts.tcondfile, modelopts.tcapfile,...
    modelopts.emodelopts, modelopts.tmodelopts, modelopts.tcapmodelopts,...
    modelopts.trange);     
funpars.Epars = econdpars; funpars.Tpars = tcondpars; funpars.Tcappars = tcappars;
params.funpars = funpars;
end

% Initialize RK4 state variables as empty arrays. (Custom RK4 solver is
% capable of integrating with hysteretic state variable, comparing two
% simulations with different timesteps for convergence studies, and
% creating its own figures; none of these features used here).
state.blank = []; prev = []; imgflag = [];

switch diffsolve
    
    %% diffsolve = 0 performs first step, steady state modeling 
    case 0
        clear('biasPointData.mat'); % Guarantee over-write

        % Calculate steady state IV curve from formulas
        SSTemp = linspace(params.T0,900,5000)';
        ETcond = electrothermalmodels(funpars,SSTemp);
        Econd = ETcond.Econd; Tcond = ETcond.Tcond; Tcap = ETcond.Tcap;
        
        % Current, voltage as f(temperature) 
        SSCurrent = sqrt((SSTemp-params.T0).*Tcond.*Econd);
        SSVolt = sqrt((SSTemp-params.T0).*Tcond./Econd);

        % Determine bias point with steepest dV/dT for oscillations:
    
        % Intermediate computation using "gamma parameters", a rescaling of
        % the slopes m_I and m_V. See DOI 10.1063/5.0070558
        gmminus = gammacalculator(econdpars, tcondpars, params.T0, SSTemp);
        gmminus(gmminus<0)=0; %Only care about NDR states

        % Calculate dV/dT of steady state V-T curve (analytic formula)
        dvdT = -1/2.*sqrt(Tcond./((SSTemp-params.T0).*Econd)).*gmminus;

        % Find where dV/dT is the steepest within NDR and get I,V,T 
        [~, steepestVindex] = min(dvdT);
        steepestT = SSTemp(steepestVindex);

        % We standardize simulations by selecting max dV/dT as bias point
        Tss = steepestT; % Bias so that steady state T = steepest T
        
        % Evaluate models at Tss to determine bias current and critical C 
        outTss = electrothermalmodels(funpars,Tss);
        GelofTss = outTss.Econd;
        GthofTss = outTss.Tcond; CthofTss = outTss.Tcap;
        [gminus,gplus] = gammacalculator(econdpars, tcondpars, params.T0, Tss);

        % Calculate current and critical capacitance at this bias
        % Recalculate bias current and voltage analytically        
        steepestV = sqrt((Tss-params.T0)*GthofTss/GelofTss);        
        steepestI = sqrt((Tss-params.T0)*GthofTss*GelofTss); 
        Ccrit = GelofTss/ (gminus*GthofTss/CthofTss);
        fcrit = (1/(2*pi))*GthofTss/CthofTss*sqrt(gminus*gplus);

        % Print bias state properties to terminal for user convenience
        format longEng;
        fprintf([modelopts.name,'\n']);
        fprintf(['Steepest voltage [V]: ', num2str(steepestV), '\n'])
        fprintf(['Corresponding current [A]: ', num2str(steepestI), '\n'])
        fprintf(['Corresponding temperature [K]: ', num2str(steepestT), '\n'])

        % If material can oscillate (some NDR state), Ccrit>0
        if Ccrit > 0 % Material can oscillate 
        % Print output and save to a file to pass to diffsolve=1 step
            save('biasPointData.mat', 'steepestI', 'steepestV', 'Tss','Ccrit','GelofTss','GthofTss','CthofTss');
            fprintf(['Critical capacitance [F] at chosen bias is ', num2str(Ccrit),'\n']);
            fprintf(['Critical frequency [Hz] at chosen bias is ', num2str(fcrit),'\n']);
        else % Material cannot oscillate
        % Print errror warning and exit
             fprintf(['This material does not have CC-NDR and cannot oscillate! Try a different material''\n']);
        end
      
        % Figures for steady-state IVT properties
        figure() % Temperature-dependent electrical conductance
        plot(SSTemp, Econd); title('Electrical Conductance vs. T');
        set(gcf,'color','w'); xlabel('Temperature (K)'); ylabel('E Conductance (S)');

        figure() % Temperature-dependent thermal conductance
        plot(SSTemp, Tcond); title('Thermal Conductance vs. T');
        set(gcf,'color','w'); xlabel('Temperature (K)'); ylabel('T Conductance (W/K)');

        figure() % Temperature-dependent thermal capacitance
        plot(SSTemp, Tcap); title('Thermal Capacitance vs. T');
        set(gcf,'color','w'); xlabel('Temperature (K)'); ylabel('T Capacitance (J/m^3-K)');

        figure() % Steady state temperature-voltage curve
        plot(SSTemp, SSVolt); title('Steady State V vs T');
        set(gcf,'color','w'); xlabel('Temperature (K)'); ylabel('SS Voltage (V)');

        figure() % Steady-state current-voltage curve
        plot(SSCurrent, SSVolt); title('Steady State V vs I');
        set(gcf,'color','w'); set(gca,'XScale','log');
        xlabel('SS Current (A)'); ylabel('SS Voltage (V)')

        figure() % T-dependent electrical conductance: Arrhenian coordinates
        plot(1000./SSTemp, Econd); title('Arrhenius Plot');
        set(gcf,'color','w'); set(gca,'YScale','log')
        xlabel('1000/T'); ylabel('E Conductance (S)');

        figure() % Slope of steady state V-T curve
        plot(SSTemp,dvdT); title('Slope of Steady State VT')
        set(gcf,'color','w'); set(gca,'XScale','log');
        xlabel('SS Temperature (T)'); ylabel('dV/dT (V/K)');
    
    case 1
        %% diffsolve = 1 performs second step, oscillation simulations 
                
        % Pass bias data computed in diffsolve=0 step
        if exist('biasPointData.mat','file') == 2
            load('biasPointData.mat');
        else
            error(['Steady state data file not found. Be sure to run ' ...
                'diffsolve=0 step before diffsolve=1 step.']);
        end

        % This block is part of renormalizing time and temperature to help
        % with convergence

        % Need to calculate model at T0, but this is a constant, so compute 
        % it just once, outside of for loops
        outT0 = electrothermalmodels(funpars,params.T0);
        params.GthofT0 = outT0.Tcond;
        params.GelofT0 = outT0.Econd;
        params.CthofT0 = outT0.Tcap;
     
        % Normalize by thermal time constant, tth = 2*pi*Cth(T0)*Rth(T0)
        % (Helps with convergence, see SI)
        thermalt = 2*pi*params.CthofT0/params.GthofT0;
        lims = lims/thermalt;

        % Variables needed for renormalizing and denormalizing
        params.V0 = sqrt(params.T0*params.GthofT0/params.GelofT0);
        params.I0 = steepestI;
        params.GelofTss = GelofTss;

        % Convert user-input capacitance factor [1] to actual ckt capacitance [F]
        params.Cel = Ccrit*params.CelFactor;
        % Standard initial conditions at 10% offset from V-T steady state
        initcond = 1.1*[round(Tss/params.T0,5); round(steepestV/params.V0,5)]; % Rounded to 5th decimal place 

        % Function call to RK4Nsolver to implement Runge-Kutta integrator
        % Equations integrated in derivfunc.m are renormalized, see SI
        tic
        [sol,~,err] = RK4Nsolver(lims,1/numsteps,initcond,prev,imgflag,state,params);        
        toc

        % De-normalize the simulated oscillations (see SI)
        sol(1,:) = thermalt*sol(1,:);
        sol(2,:) = params.T0*sol(2,:);
        sol(3,:) = params.V0*sol(3,:);

        figure() % Time-dependent temperature
        plot(sol(1,:),sol(2,:));
        set(gca,'XScale','linear'); set(gca,'YScale','log');
        ylabel('Temperature (K)'); xlabel('Time (s)');
        title('Temperature Trajectory'); set(gcf,'color','w');
        
        figure() % Time-dependent voltage
        plot(sol(1,:),sol(3,:))
        set(gca,'XScale','linear'); set(gca,'YScale','log');
        ylabel('Voltage (V)'); xlabel('Time (s)');
        title('Voltage Trajectory'); set(gcf,'color','w');

        % Compute minT [K], maxT [K], minV [V], maxV [V], and freq [Hz]
        oscResults = oscAnalyzer(sol);
        
        % Calculate power [W] at bias point 
        P = GthofTss*(Tss-params.T0);
        
        % Calculate deltaT and deltaV
        deltaT = oscResults(2)-oscResults(1);
        deltaV = oscResults(4)-oscResults(3);
        
        %  Results: 
        % Cckt, Tmin, Tmax, Vmin, Vmax, Freq, Power, deltaT, deltaV, 
        % Print to output and write to variable "result"
        result = [params.CelFactor*Ccrit, oscResults,P,deltaT,deltaV];
        format longEng;
        fprintf([modelopts.name,'\n']);
        fprintf(['Circuit capacitance [F]: ', num2str(result(1)), '\n']);
        fprintf(['Frequency [Hz]: ', num2str(result(6)), '\n']);
        fprintf(['Power [W]: ', num2str(result(7)), '\n']);
        fprintf(['T amplitude [K]: ', num2str(result(8)), '\n']);
        fprintf(['V amplitude [V]: ', num2str(result(9)), '\n']);        
end

%% Subfunctions Block

function [gminus, gplus] = gammacalculator(econdpars, tcondpars, T0, T)
% Function "gammacalculator.m" v2.0.0, tested 18 September 2024
% Written by T.D. Brown Fall 2023
%
% This function computes what are called "gamma parameters" for the
% electrothermal memristor as a function of temperature. These gammas are a
% unitless rescaling of the stability slopes m_I and m_V, see DOI 10.1063/5.0070558
%
% Material models are differentiated analytically to eliminate numerical
% discretization noise
%
% Input: econdpars = vector of model parameters for T-dependent electrical conductivity
% Input: tcondpars = vector of model parameters for T-dependent thermal conductivity
% Input: T0 = scalar of ambient temperature reservoir (typically 300 K) [K]
% Input: T = vector of temperatures to evaluate gamma parameters at [K]
% Output: gminus = vector of gamma- parameters [1]
% Output: gplus = vector of gamma+ parameters [1]

    % If T is a vector, ensure it is a column
    if size(T,1) < size(T,2)
        T = T';
    end

    % Electrical and thermal log derivatives (computed analytically)
    econddlog = [T, PiecewiseArrhenian(0, econdpars, T, 1)];
    tconddlog = [T, PiecewiseLaurent(0, tcondpars, T, 1)];

    % Gamma parameters
    gplus = (tconddlog(:,1)-T0).*(econddlog(:,2)+tconddlog(:,2)) + 1;
    gminus = (tconddlog(:,1)-T0).*(econddlog(:,2)-tconddlog(:,2)) - 1;

end

function result = oscAnalyzer(in)
% Function "oscAnalyzer.m" v2.0.0, tested 18 September 2024
% Written by T.D. Brown Fall 2023
%
% This function computes relevant oscillation properties from a given
% oscillation time series for temperature and voltage. Returns the
% temperature extrema (min and max) and voltage extrema, and oscillation
% frequency. Works for nonlinear oscillations and automatically waits 90%
% of total time to ensure transients are damped out.
%
% Input: in = array of temperature-voltage time series (rows: time, T, v)
% Output: result = vector of oscillation properties: extrema + frequency

    % Remove all but the last 10% of the data
    cutoff = round(0.9*size(in,2));
    dat = in(:,cutoff:size(in,2));

    % Maxima analysis
    % Thermal oscilations
    [pks,locs] = findpeaks(dat(2,:),dat(1,:));
    perd = mean(diff(locs)); % Oscillation period [s]
    freq = 1/perd;
    tmax = mean(pks);
    % Voltage oscillations
    pks = findpeaks(dat(3,:),dat(1,:));
    vmax = mean(pks);

    % Minima analysis. Have to flip signal upside-down
    dat2 = dat;
    dat2(2,:) = 1.5*max(dat(2,:)) - dat2(2,:);
    dat2(3,:) = 1.5*max(dat(3,:)) - dat2(3,:);
    % Thermal oscillations
    pks = findpeaks(dat2(2,:),dat2(1,:));
    tmin = 1.5*max(dat(2,:))-mean(pks);
    % Voltage oscillations
    pks = findpeaks(dat2(3,:),dat2(1,:));
    vmin = 1.5*max(dat(3,:))-mean(pks);

    % Return oscillation results in "result"
    result = [tmin, tmax, vmin, vmax, freq];
end

function modelopts = selectMaterial(material_select)
% Function "selectMaterial.m" v2.0.0, tested 18 September 2024
% Written by T.D. Brown Spring 2024
%
% This function is essentially a wrapper for a switch-case structure to
% choose between 40 pre-loaded material datasets and models
%
% Input: material_select [Int] = flag corresponding to desired material
% Output: modelopts = struct containing needed variables for fitting model

    % Switch between 40 pre-loaded material datasets and models
    switch material_select         
        case 1 % AlN
        trange = linspace(200, 1400, 5000); name = 'AlN';
        econdfile = 'AlNECond.xlsx'; tcondfile = 'AlNTCond.xlsx'; tcapfile = 'AlNCth.xlsx';
        emodelopts = [1,0.4,0.6]; tmodelopts = [0, 0, 0]; tcapmodelopts = [1, 0.7, 0.99];
    
        case 2 %B-Ga2O3 !!! Electrical model
        trange = linspace(300,1000, 5000); name = 'B-Ga2O3';
        econdfile = 'Ga2O3ECond12.xlsx'; tcondfile = 'B-Ga2O3TCond.xlsx'; tcapfile = 'B-Ga2O3Cth.xlsx';
        emodelopts = [1,0.2,0.6]; tmodelopts = [1, 0.5, 0.6]; tcapmodelopts = [0,0.5,0.7];
    
        case 3 % Bi2Te3
        trange = linspace(100, 1400, 5000); name = 'Bi2Te3';
        econdfile = 'Bi2Te3ECond2.xlsx'; tcondfile = 'Bi2Te3TCond.xlsx'; tcapfile = 'Bi2Te3Cth.xlsx';
        emodelopts = [1,0.55,0.65]; tmodelopts = [1, 0.7, 0.8]; tcapmodelopts = [0, 0.7, 0.99];
    
        case 4 % CdS
        trange = linspace(100, 600, 5000); name = 'CdS';
        econdfile = 'CdSECondInterp.xlsx'; tcondfile = 'CdSTCond.xlsx'; tcapfile = 'CdSCth.xlsx';
        emodelopts = [1,0.2,0.9]; tmodelopts = [1, 0.5, 0.6]; tcapmodelopts = [1,0.5,0.7];
        
        case 5 % CdSe
        trange = linspace(300,1000, 5000); name = 'CdSe';
        tcondfile = 'CdSeTCondBulk.xlsx'; econdfile = 'CdSeECond'; tcapfile = 'CdSeCv.xlsx';
        emodelopts = [1,0.4,0.6]; tmodelopts = [0,0,0]; tcapmodelopts = [1,0.4,0.9];
    
        case 6 % Co3O4 
        trange = linspace(200,1000, 5000); name = 'Co3O4';
        tcondfile = 'Co3O4TCond';  econdfile = 'Co3O4ECond'; tcapfile = 'Co3O4Cth';
        emodelopts = [1,0.3,0.5]; tmodelopts = [0, 0.3, 0.6]; tcapmodelopts = [1,0.5,0.7];
        
        case 7 % CuO
        trange = linspace(300, 1000, 5000); name = 'CuO';
        econdfile = 'CuOECond.xlsx'; tcondfile = 'CuOTCond.xlsx'; tcapfile = 'CuOCth.xlsx';
        emodelopts = [1,0.5,0.55]; tmodelopts = [0,0,0]; tcapmodelopts = [0, 0, 0];
    
        case 8 % Fe3O4
        trange = linspace(200, 1200, 5000); name = 'Fe3O4';
        econdfile = 'Fe3O4ECondCrystal.xlsx'; tcondfile = 'Fe3O4TCondPolycrystal.xlsx'; tcapfile = 'Fe3O4Cth.xlsx';
        emodelopts = [1,0.3,0.6]; tmodelopts = [0,0,0]; tcapmodelopts = [0,0,0];
    
        case 9 % FeO 
        trange = linspace(200,1000, 5000); name = 'FeO';
        tcondfile = 'FeOTCond'; econdfile = 'FeOECond'; tcapfile = 'FeOCth';
        emodelopts = [1,0.3,0.5]; tmodelopts = [0, 0.4, 0.8]; tcapmodelopts = [1,0.5,0.7];
    
        case 10 % GaAs
        trange = linspace(100,1000, 5000); name = 'GaAs';
        econdfile = 'GaAsECondSample3.xlsx'; tcondfile = 'GaAsTcond4.xlsx'; tcapfile = 'GaAsCth1.xlsx';
        emodelopts = [0,0,0.0]; tmodelopts = [0, 0, 0]; tcapmodelopts = [1, 0.85, 0.9];
    
        case 11 % GaN
        trange = linspace(200, 800, 5000); name = 'GaN';
        tcondfile = 'GaNTCondBulk.xlsx'; econdfile = 'GaNECond.xlsx'; tcapfile = 'GaNCth.xlsx';
        emodelopts = [0,0.3,0.4]; tmodelopts = [0,0,0]; tcapmodelopts = [0,0,0];
    
        case 12 % GaP
        trange = linspace(200, 800, 5000); name = 'GaP';
        econdfile = 'GaPECond.xlsx'; tcondfile = 'GaPTCond2.xlsx'; tcapfile = 'GaPCth.xlsx';
        emodelopts = [0,0,0]; tmodelopts = [0,0,0]; tcapmodelopts = [0,0,0];
    
        case 13 % GaSb
        trange = linspace(200, 1200, 5000); name = 'GaSb';
        econdfile = 'GaSbECond.xlsx'; tcondfile = 'GaSbTCond.xlsx'; tcapfile = 'GaSbCth2.xlsx';
        emodelopts = [1,0.4,0.55]; tmodelopts = [0,0,0]; tcapmodelopts = [1,0.65,0.8];
    
        case 14 % GaSe
        trange = linspace(50, 1000, 5000);name = 'GaSe';
        econdfile = 'GaSeECond1.xlsx'; tcondfile = 'GaSeTCondPerp.xlsx'; tcapfile = 'GaSeCth.xlsx';
        emodelopts = [0,0,0]; tmodelopts = [0,0.1,0.6]; tcapmodelopts = [1,0.78,0.88];
    
        case 15 % Ge
        trange = linspace(300,1000, 5000); name = 'Ge';
        econdfile = 'GeEcond.xlsx'; tcondfile = 'GeTCond2.xlsx'; tcapfile = 'GeCth.xlsx';
        emodelopts = [1,0.1,0.6]; tmodelopts = [0, 0.0, 0.0]; tcapmodelopts = [0, 0.55, 0.9];
    
        case 16 % InSb
        trange = linspace(100, 1400, 5000); name = 'InSb';
        econdfile = 'InSbECondC.xlsx'; tcondfile = 'InSbTCond.xlsx'; tcapfile = 'InSbCth.xlsx';
        emodelopts = [1,0.4,0.6]; tmodelopts = [1,0.6,0.7]; tcapmodelopts = [1,0.5,0.7];
    
        case 17 % InSe
        trange = linspace(30,800, 5000); name = 'InSe';
        econdfile = 'InSeECond.xlsx'; tcondfile = 'InSeTCond.xlsx'; tcapfile = 'InSeCp.xlsx';
        emodelopts = [0,0.0,0.0]; tmodelopts = [0, 0.0, 0.0]; tcapmodelopts = [0, 0, 0];
    
        case 18 % LCO
        trange = linspace(200, 1000, 5000); name = 'LaCoO3';
        econdfile = 'LCOECond.xlsx'; tcondfile = 'LCOLitTCond.xlsx'; tcapfile = 'LCO_Cp.xlsx';
        emodelopts = [1,0.1,0.6]; tmodelopts = [1, 0.05, 0.12]; tcapmodelopts = [1,0.55,0.9];
    
        case 19 % MoO3 %
        trange = linspace(300, 800, 5000); name = 'MoO3';
        econdfile = 'MoO3ECondHeatingA.xlsx'; tcondfile = 'MoO3TCond.xlsx'; tcapfile = 'MoO3Cth.xlsx';
        emodelopts = [1,0.45,0.6]; tmodelopts = [0,0,0]; tcapmodelopts = [1,0.6,0.9];
    
        case 20 % MoS2 %
        trange = linspace(50, 1000, 5000); name = 'MoS2';
        econdfile = 'MoS2ECond.xlsx'; tcondfile = 'MoS2Thermal.xlsx'; tcapfile = 'MoS2Cth.xlsx';
        emodelopts = [0,0,0]; tmodelopts = [0,0,0]; tcapmodelopts = [1,0.6,0.9];
    
        case 21 % MoSe2
        trange = linspace(200, 1200, 5000); name = 'MoSe2';
        econdfile = 'MoSe2ECond2.xlsx'; tcondfile = 'MoSe2TCond.xlsx'; tcapfile = 'MoSe2Cth.xlsx';
        emodelopts = [0,0,0]; tmodelopts = [0,0,0]; tcapmodelopts = [0,0,0];
    
        case 22 % NbO2 %
        trange = linspace(300, 1000, 5000); name = 'NbO2';
        econdfile = 'NbO2Electrical.xlsx'; tcondfile = 'NbO2TCond.xlsx'; tcapfile = 'NbO2Cth.xlsx';
        emodelopts = [1, 0.92,0.95]; tmodelopts = [0,0,0]; tcapmodelopts = [1, 0.5, 0.9];
    
        case 23 % NdNiO3 
        trange = linspace(300,1000, 5000); name = 'NdNiO3';
        tcondfile = 'NdNiO3_electrical.xlsx'; econdfile = 'NdNiO3_thermal'; tcapfile = 'NdNiO3Cth';
        emodelopts = [1,0.2,0.9]; tmodelopts = [1, 0.3, 0.8]; tcapmodelopts = [0,0.5,0.7];
    
        case 24 % PbTe
        trange = linspace(200, 1000, 5000); name = 'PbTe';
        econdfile = 'PbTeECond.xlsx'; tcondfile = 'PbTeTCond.xlsx'; tcapfile = 'PbTeCth.xlsx';
        emodelopts = [0,0,0]; tmodelopts = [0, 0, 0]; tcapmodelopts = [0, 0,0];
    
        case 25 % Sb2Te3 
        trange = linspace(300, 700, 5000); name = 'Sb2Te3';
        econdfile = 'Sb2Te3ECondCrystal.xlsx'; tcondfile = 'Sb2Te3TCondCrystal.xlsx'; tcapfile = 'Sb2Te3Cth.xlsx';
        emodelopts = [1,0.2,0.8]; tmodelopts = [1, 0.2, 0.9]; tcapmodelopts = [0, 0, 0];
    
        case 26 % Si
        trange = linspace(300, 1400, 5000); name = 'Si';
        econdfile = 'SiECond_Sample1.xlsx'; tcondfile = 'SiTCondSamples.xlsx'; tcapfile = 'SiCth';
        emodelopts = [1,0.4,0.6]; tmodelopts = [0,0,0]; tcapmodelopts = [0, 0, 0];
    
        case 27 % SiC
        trange = linspace(200, 1400, 5000); name = 'SiC';
        econdfile = 'SiCECond3.xlsx'; tcondfile = 'SiCTCond2.xlsx'; tcapfile = 'SiCCth.xlsx';
        emodelopts = [0,0,0]; tmodelopts = [0, 0, 0]; tcapmodelopts = [0, 0, 0];
    
        case 28 % SnO2
        trange = linspace(100, 1400, 5000); name = 'SnO2';
        econdfile = 'SnO2ECond41.xlsx'; tcondfile = 'SnO2TCond.xlsx'; tcapfile = 'SnO2Cth.xlsx';
        emodelopts = [0,0.2, 0.8]; tmodelopts = [0, 0, 0]; tcapmodelopts = [0,0,0];
    
        case 29 % SnSe
        trange = linspace(300,1000, 5000); name = 'SnSe';
        tcondfile = 'SnSeTCond.xlsx'; econdfile = 'SnSeECond.xlsx'; tcapfile = 'SnSeCth.xlsx';
        emodelopts = [1,0.3,0.7]; tmodelopts = [0, 0, 0]; tcapmodelopts = [1,0.5,0.6];
    
        case 30 % STO
        trange = linspace(200, 1200, 5000); name = 'SrTiO3';
        econdfile = 'SrTiO3ECond.xlsx'; tcondfile = 'SrTiO3TCond.xlsx'; tcapfile = 'SrTiO3Cth.xlsx';
        emodelopts = [1,0.4,0.6]; tmodelopts = [0,0,0]; tcapmodelopts = [1,0.45,0.6];
    
        case 31 % TaS2 !!! Cth file missing
        trange = linspace(300,1000, 5000); name = 'TaS2';
        tcondfile = 'TaS2TCond1'; econdfile = 'TaS2ECond1'; tcapfile = 'TaS2Cth';
        emodelopts = [1,0.2,0.4]; tmodelopts = [0, 0.3, 0.6]; tcapmodelopts = [1,0.5,0.7];
    
        case 32 % TbMnO3
        trange = linspace(100, 1000, 5000); name = 'TbMnO3';
        econdfile = 'TbMnO3ECond.xlsx'; tcondfile = 'TbMnO3Thermal.xlsx'; tcapfile = 'TbMnO3Cth.xlsx';
        emodelopts = [0,0,0]; tmodelopts = [0,0,0]; tcapmodelopts = [1,0.2,0.7];
    
        case 33 % TiO2 !!! Files not uploaded
        trange = linspace(200,1000, 5000); name = 'TiO2';
        tcondfile = 'TiO2TCond'; econdfile = 'TiO2ECond'; tcapfile = 'TiO2Cth';
        emodelopts = [1,0.3,0.5]; tmodelopts = [1, 0.3, 0.6]; tcapmodelopts = [1,0.5,0.7];
    
        case 34 %V3O5
         trange = linspace(200, 500, 5000); name = 'V3O5';
        econdfile = 'V3O5Electrical.xlsx'; tcondfile = 'V3O5Thermal.xlsx'; tcapfile = 'V3O5Cth.xlsx';
        emodelopts = [1,0.1,0.6]; tmodelopts = [1,0.4,0.6]; tcapmodelopts = [1,0.5,0.9];
    
        case 35 % VO2 
        trange = linspace(300, 1000, 5000); name = 'VO2';
        tcondfile = 'VO2TCond.xlsx'; econdfile = 'VO2ECond'; tcapfile = 'VO2Cth';
        emodelopts = [1,0.2,0.9]; tmodelopts = [1, 0.3, 0.8]; tcapmodelopts = [1,0.5,0.7];
    
        case 36 % WO3 !!!
        trange = linspace(300, 1000, 5000); name = 'WO3';
        econdfile = 'WO3ECond.xlsx'; tcondfile = 'WO3TCond.xlsx'; tcapfile = 'WO3Cth.xlsx';
        emodelopts = [0,0,0]; tmodelopts = [0, 0,0]; tcapmodelopts = [1, 0.5, 0.6];
    
        case 37 % WSe2
        trange = linspace(200, 800, 5000); name = 'WSe2';
        econdfile = 'WSe2ECond2.xlsx'; tcondfile = 'WSe2TCond.xlsx'; tcapfile = 'WSe2Cth.xlsx';
        emodelopts = [0,0,0]; tmodelopts = [0,0,0]; tcapmodelopts = [0,0.6,0.9];        
    
        case 38 % ZnO 
        trange = linspace(200,1200, 5000); name = 'ZnO';
        tcondfile = 'ZnOTCond.xlsx'; econdfile = 'ZnOEcond2.xlsx'; tcapfile = 'ZnOCth.xlsx';
        emodelopts = [1,0.3,0.7]; tmodelopts = [0,0,0]; tcapmodelopts = [1, 0.3, 0.8];
    
        case 39 % ZnSe
        trange = linspace(200,1400, 5000); name = 'ZnSe';
        econdfile = 'ZnSeECond1.xlsx'; tcondfile = 'ZnSeTCond.xlsx'; tcapfile = 'ZnSeCth.xlsx';
        emodelopts = [0,0.1,0.3]; tmodelopts = [0, 0, 0]; tcapmodelopts = [1,0.5,0.9];
        
        case 40 % ZnTe !!! Shows NDR predicted
        trange = linspace(100, 800, 5000); name = 'ZnTe';
        econdfile = 'ZnTeECondCrystal.xlsx'; tcondfile = 'ZnTeTCond.xlsx'; tcapfile = 'ZnTeCth.xlsx';
        emodelopts = [1,0.3,0.7]; tmodelopts = [1,0.32,0.42]; tcapmodelopts = [1,0.6,0.8];
    
    
        case 41 % InN
        trange = linspace(100, 1400, 5000);
        econdfile = 'InNECond.xlsx'; tcondfile = 'InNTCond.xlsx'; tcapfile = 'DummyCth.xlsx';
        emodelopts = [0,0,0]; tmodelopts = [0, 0, 0]; tcapmodelopts = [0,0,0];
    end
    
    % Wrap up output in struct modelopts
    modelopts.trange = trange; modelopts.econdfile = econdfile; 
    modelopts.tcondfile = tcondfile; modelopts.tcapfile = tcapfile;
    modelopts.emodelopts = emodelopts; modelopts.tmodelopts = tmodelopts; 
    modelopts.tcapmodelopts = tcapmodelopts; modelopts.name = name;


end