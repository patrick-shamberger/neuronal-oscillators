%  plotAutomatedAshbyDesigns.m

% Script "plotAutomatedAshbyDesigns.m" v2.0.0, tested 25 September
% 2024 
% Written by J. Chong Spring 2024 
% 
% This script plots Ashby-inspired plots based on simulation data generated
% from the compact model, main_NeuronalMaterialsSimulations.m. These plots
% are generated for particular design constraints; all simulations must be
% collected at the exact same ambient temperature T0, and device dimension
% of h (thickness) and R (radius). Later on, the user will have the ability to input new 
% device dimensions to see how performance changes. The user can also control 
% what design capacitance or the design frequency at which meaningful comparisons
% will be made 

% The program requires two spreadsheets: 
% (1) AshbyPlotData.xlsx, containing predetermined simulation data
% for 40 materials systems. It is listed under the variable
% "OriginalData" 
% (2) A spreadsheet containing the user's compiled simulation data from the 
% compact model, which is called "NewData.xlsx". It is listed under
% the variable "filename"
% 
% Both spreadsheets must contain 8 columns, with
% the first row labeled with "Material", "Capacitance", "Frequency", etc. 
% The "result" variable's values from the main script should be copy-pasted
% into the spreadsheet starting at the Capacitance column. The spreadsheet 
% should resemble "AshbyPlotData.xlsx" and follow the format below: 

%           Material    Capacitance     Frequency       Power       DeltaT      DeltaV       MinT       MaxT
%          __________    __________    __________    __________    ________    ________    ________    ________
%
%             'VO2'         10E9              1e+09         100        5           0.5        10           90   
%             'VO2'         20E9              2e+09         150       81           0.6        15           85   
%             'VO2'         15E9            1.5e+09         120       62           0.4        12           88   
%
%             'NbO2'        7E11            4.5e+07         1E3        5           0.5        10           90   
%             'NbO2'        22E11             6e+09         1E5       68           0.6        15           85   
%             'NbO2'        14E12           1.8e+09         3E4       32           0.4        12           88   


% Input [control]: OriginalorNewDimz = Boolean; (0) plot everything in
% original dimensions of simulations, (1) plot in the user's defined dimensions 
% Input [control]: newdims.Radius = user defined radius, R, of cylindrical
% thru-plane device [m]
% Input [control]: newdims.Thickness = user defined thickness, h, of cylindrical
% thru-plane device [m]

% Input [control]: FrequencyorCapacitanceConstraint = Boolean; (0) design
% capacitance is the circuit capacitance. (1) design capacitance is the oscillation
% frequency 
% Input [control]: designConstraint.Capacitance = user defined design
% constraint for circuit capacitance [F]
% Input [control]: designConstraint.Frequency = design constraint for oscillation frequency [Hz]

% When FrequencyorCapacitanceConstraint == 0
% Output: PerformanceMetrics variable contains array of materials and their
% performance metrics at this capacitance constraint 
% Output: plots frequency and deltaV/deltaT as a function of circuit
% capacitance for every material at the specified dimensions. A vertical line
% is drawn where the design capacitance is defined. All predictions of materials 
% intercepting this line are projected onto graphs comparing frequency,
% deltaV/deltaT, and power at the design constraint.

% When FrequencyorCapacitanceConstraint == 1
% Output: PerformanceMetrics variable contains array of materials and their
% performance metrics at this frequency constraint 
% Output: plots frequency as a function of circuit capacitance for every
% material at the specified dimensions. A horizontal line is drawn where
% the design frequency is defined. All predictions of materials intercepting
% this line are projected onto graphs comparing frequency, deltaV/deltaT,
% and power at the design constraint.

clear all;
close all;

%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ USER-MODIFIED BLOCK ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
filename = '\NewData.xlsx'; % spreadsheet must be named "NewData.xlsx" and should be in the same directory as this script 
OriginalorNewDimz = 0; % Boolean: if 0, plot everything in original dims, if 1, plot in new dims 

newdims.Radius = 1e-7; % evaluate properties at this new radius of device using scaling laws
newdims.Thickness = 1e-7; % evaluate properties at this new thickness of device using scaling laws

FrequencyorCapacitanceConstraint = 0; % Boolean: if 0, circuit capacitance is the design constraint.
% if 1, oscillation frequency is the design constraint 

designConstraints.Capacitance = 1e-10; % design constraint for circuit capacitance [F]
designConstraints.Frequency = 5e6; % design constraint for oscillation frequency [Hz]
%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

PerformanceMetrics = plotAshby(filename, newdims,designConstraints,OriginalorNewDimz,FrequencyorCapacitanceConstraint);

function results = plotAshby(filename, newdims, constraints,OriginalorNew,FrequencyorCapacitanceConstraint)

    % Automatically scales original simulated data and creates plots of
    % interest
    
    fileroot = pwd; %! Need root directory for relative path
    rng(55); % Choose any integer as the seed for picking random colors (e.g., 42)
    
    % DO NOT CHANGE OLD DIMENSIONS UNLESS YOU ARE NOT USING
    % AshbyPlotData.XLSX AS ORIGINAL DATA!! 
    olddims.Radius = 1e-6; % Dimensions (device Radius) used for original simulations
    olddims.Thickness = 1e-7; % Dimensions (device Thickness) used for original simulations
    
    % Define the data. Color codes for materials in AshbyPlotData.xlsx. 
    hexcodes = {
             'AlN', '#a5be00','1';       %Thermoelectrics
             'CdSe', '#4e148c','2';      %Chalcogenide
             'Co_3O_4', '#b23a48','3';   %Oxide    
             'CuO', '#d00000','4';       %Oxide
             'FeO', '#38040e','5';       %Oxide
             'GaAs', '#007f5f','6';      %Semiconductor
             'GaP', '#5bba6f','7';       %Semiconductor
             'GaSb', '#fbb02d','8';      %Thermoelectrics
             'GaSe', '#0e6ba8','9';      %Chalcogenide
             'Ge', '#538d22','10';       %Semiconductor
             'InSb', '#80ed99','11';     %Semiconductor
             'InSe', '#0a2472','12';     %Chalcogenide
             'LCO', '#F95738','13';      %perovskite
             'MoO_3','#da627d','14';     %Oxide
             'MoS_2', '#9bb1ff','15';    %Chalcogenide
             'MoSe_2', '#788bff','16';   %Chalcogenide
             'NbO_2', '#e85d04','17';    %Oxide
             'SrTiO_3 ','#684756','18';  %Perovskite             
             'Si', '#1b4332','19';       %Semiconductor
             'SiC', '#D4E755','20';      %Semiconductor
             'SnSe', '#5465ff','21';     %chalcogenide
             'TbMnO_3', '#ec0868','22';  %perovskite
             'V_3O_5', '#EA5F92','23';   %Oxide
             'VO_2', '#ce796b','24';     %Oxide  
             'WO_3', '#ffadc6','25';     %Oxide
             'WSe_2', '#63B0CD','26';    %chalcogenide
             'ZnO', '#690500','27';      %Oxide
             'ZnSe', '#00a6fb','28';     %Chalcogenide
        }; 
    
    % Save the data to a MAT file
    save('hexcodes.mat', 'hexcodes');
    
    % Load, group, and combine simulation data
    originalData = readtable([fileroot,'\AshbyPlotData.xlsx'],ReadVariableNames=true);
    
    newData = readtable([fileroot,filename]);
    OscillatorProperties = GroupData(originalData);
    NewOscillatorProperties = GroupData(newData);
    CombinedOscillatorProperties = [OscillatorProperties;NewOscillatorProperties];
    
    % Initialize cells for 500 K limited and scaled simulation data
    LimOscillatorProperties = cell(length(CombinedOscillatorProperties),1);
    ScaledOscillatorProperties = cell(length(CombinedOscillatorProperties), 1);
    
    % if user does not use AshbyPlotData.xlsx, make new color codes 
    if isempty(originalData)
        colors.hexcodes = {};
    else 
        colors = load('hexcodes.mat');
    end
    
    % % if using AshbyPlotData data, append new color codes to the end of the original list. In alphabetical order 
    for j= length(colors.hexcodes)+1:length(LimOscillatorProperties)
        colors.hexcodes(j,1) = CombinedOscillatorProperties{j,1}.Material(1);
        % Generate a random hex code for the new string
        hex_values = dec2hex(randi([0, 255], 1, 3), 2);
        colors.hexcodes(j,2) ={['#' char(strjoin(cellstr(hex_values), ''))]};
    end 
    
    % Limit deltaT to 1 - 500 K 
    for j = 1:length(LimOscillatorProperties)
        LimOscillatorProperties{j,1} = LimitDeltaTOscillations(CombinedOscillatorProperties{j,1});
        % Extrapolate to higher capacitances 
        LimOscillatorProperties{j,1} = extrapolate(LimOscillatorProperties{j,1},constraints);
    
        if OriginalorNew == 1 
        % Calculate Power, Capacitance, and Frequency for different dimensions 
        ScaledOscillatorProperties{j,1} = calculateScaledProperties(LimOscillatorProperties{j,1},olddims,newdims);
        ScaledOscillatorProperties{j,1} = extrapolate(ScaledOscillatorProperties{j,1},constraints);
        end
    end
    

    if OriginalorNew == 1 
        results = generatePlots(ScaledOscillatorProperties, colors.hexcodes,constraints,FrequencyorCapacitanceConstraint);
    else 
        results = generatePlots(LimOscillatorProperties, colors.hexcodes,constraints,FrequencyorCapacitanceConstraint);
    end
end

%% Functions 
function OscillatorProperties = GroupData(data)
    % If new data file is empty, run plots with only the old data 
    if isempty(data)
        OscillatorProperties = [];
        return
    end
    % Split the data based on the values in the first column
    groups = findgroups(data{:, 1});
    
    % Get unique values in the first column
    unique_values = unique(data{:, 1});
     unique_values = unique_values(~cellfun('isempty', unique_values));

    % Initialize a cell array to store data tables for each unique value
    OscillatorProperties = cell(length(unique_values), 1);
    
    % Loop through each unique value and create a data table for it
    for i = 1:length(unique_values)
        idx = groups == i; % Indices corresponding to the current unique value
        OscillatorProperties{i} = data(idx, :);
    end

end 

function extrap_material = extrapolate(material, designConstraints)
% Extrapolates computed f vs C curves using power law and deltaV/deltaT
% vs C curves 

    % Unpackage input struct
    Capacitance = material.Capacitance;
    Frequency = material.Frequency;
    DeltaV = material.DeltaV;
    DeltaT = material.DeltaT;
    DesignCapacitance = designConstraints.Capacitance;
    DesignFrequency = designConstraints.Frequency;

    % Extrapolate data on log f - log C plot
    % Fit a linear model on the log-transformed data
    logcap = log(Capacitance);
    logfreq = log(Frequency);
    coefficients = polyfit(logcap, logfreq, 1);
    slope = coefficients(1);
    intercept = coefficients(2);

    % Extrapolate deltaV/deltaT 
    logVT = log(DeltaV./DeltaT);
    VTcoefficients = polyfit(logcap, logVT, 1);
    VTslope = VTcoefficients(1);
    VTintercept = VTcoefficients(2);

    
    EndCap = 10^-5; % End of capacitance axis in F
    interpolate_freq = []; % interpolated frequency 
    interpolate_VT = []; % interpolated deltaV/deltaT 

    % Extrapolate using the fitted model
    if Capacitance(end) < DesignCapacitance
        extrapolate_cap = [Capacitance(end), DesignCapacitance, EndCap];
         % if Ccrit is larger than EndCap, then don't extrapolate 
    elseif Capacitance(1) > EndCap
        extrapolate_cap = NaN;
    else
        extrapolate_cap = [Capacitance(end), EndCap];
        % interpolate to design capacitance if within capacitance range
        if Capacitance(1) < DesignCapacitance 
            interpolate_cap = DesignCapacitance;
            interpolate_freq = exp(intercept) * interpolate_cap.^slope;
            interpolate_VT = exp(VTintercept) * interpolate_cap.^VTslope;
        end
    end 
    extrapolate_freq= exp(intercept) * extrapolate_cap.^slope;
    extrapolate_VT = exp(VTintercept) * extrapolate_cap.^VTslope;

     % solve for capacitance at the frequency from design constraint 
    FreqDesignCap = (DesignFrequency/exp(intercept))^(1/slope);
    % solve for vt at the frequency from design constraint 
    FreqDesignVT = exp(VTintercept) * FreqDesignCap.^VTslope;

    % Repackage output into struct
    extrap_material = material;
    extrap_material.extrapolatedcap = extrapolate_cap;
    extrap_material.extrapolatedfreq = extrapolate_freq;
    extrap_material.interpolatedfreq = interpolate_freq;
    extrap_material.extrapolatedVT = extrapolate_VT;
    extrap_material.interpolatedVT = interpolate_VT;
    extrap_material.freqdesigncap = FreqDesignCap;
    extrap_material.freqdesignVT = FreqDesignVT;
end 

function ScaledOscillatorProperties = calculateScaledProperties(OscillatorProperties, olddims, newdims)
% Scales power, freq, and cap for other device dimensions

    % Unpackage input structs
    power = OscillatorProperties.Power;
    freq = OscillatorProperties.Frequency;
    cap = OscillatorProperties.Capacitance;
    deltaV = OscillatorProperties.DeltaV;
    deltaT = OscillatorProperties.DeltaT; 

    R = olddims.Radius; t = olddims.Thickness;
    R_new = newdims.Radius; t_new = newdims.Thickness;
    
    % Frequency ~ 1/R^2
    FScaled = freq / (R_new/R)^2;
    
    % Power ~ t
    PScaled = power * (t_new/t);
    
    % capacitance ~ R^4/t
    CScaled = cap * (R_new/R)^4 / (t_new/t);

    % deltaV ~ t/R 
    deltaVScaled = deltaV * (t_new/t) / (R_new/R);
   
    % deltaT = delta T (ind of dimensions) 
    deltaTScaled = deltaT; 

    % Repackage output struct
    ScaledOscillatorProperties = OscillatorProperties;
    ScaledOscillatorProperties.Power = PScaled;
    ScaledOscillatorProperties.Frequency = FScaled;
    ScaledOscillatorProperties.Capacitance = CScaled;
    ScaledOscillatorProperties.DeltaV = deltaVScaled;
    ScaledOscillatorProperties.DeltaT = deltaTScaled;

end

function LimitedOscillatorProperties = LimitDeltaTOscillations(OscillatorProperties)
% Limits data to < deltaT 500 K 

    % Unpackage Input Struc
    Material = OscillatorProperties.Material;
    Freq = OscillatorProperties.Frequency;
    MaxT = OscillatorProperties.MaxT;
    MinT = OscillatorProperties.MinT;
    MaxV = OscillatorProperties.MaxV;
    MinV = OscillatorProperties.MinV;
    Cckt = OscillatorProperties.Capacitance;
    DeltaV = OscillatorProperties.DeltaV;
    DeltaT = OscillatorProperties.DeltaT;
    Power = OscillatorProperties.Power;
    
    materialName = string(Material(1));
    
    % Find indices with deltaT < 1 K
    indices_under_1 = find(DeltaT < 1);
    
    % Find the first index with deltaT exceeding 900 K
    [~, Index900] = min(abs(DeltaT - 900));    
    
    % if there are no indices under 1 K, start from index 1 
    if isempty(indices_under_1) 
        indices_under_1=0;
    end
    if isempty(Index900) == 1
        deltaTResults = DeltaT(indices_under_1(end)+1:end);
        deltaVResults = DeltaV(indices_under_1(end)+1:end);
        CapResults = Cckt(indices_under_1(end)+1:end);
        MinTResults = MinT(indices_under_1(end)+1:end);
        MaxTResults = MaxT(indices_under_1(end)+1:end);
        MinVResults = MinV(indices_under_1(end)+1:end);
        MaxVResults = MaxV(indices_under_1(end)+1:end);
        freqResults = Freq(indices_under_1(end)+1:end);
    else 
        % Remove indices under 1 and above the first index above 500
        deltaTResults = DeltaT(indices_under_1(end)+1:Index900);
        deltaVResults = DeltaV(indices_under_1(end)+1:Index900);
        CapResults = Cckt(indices_under_1(end)+1:Index900);
        MinTResults = MinT(indices_under_1(end)+1:Index900);
        MaxTResults = MaxT(indices_under_1(end)+1:Index900);
        MinVResults = MinV(indices_under_1(end)+1:Index900);
        MaxVResults = MaxV(indices_under_1(end)+1:Index900);
        freqResults = Freq(indices_under_1(end)+1:Index900);
    
    end


    % Find the index of the temperature closest to 500 K
    [~, Index500] = min(abs(deltaTResults - 500));
    if isempty(Index500) 
        Index500 = Index900;
    end

    % Find the first index closest to 300 K 
    [~, Index300] = min(abs(deltaTResults - 300));
    if isempty(Index300) 
        Index300 = Index500;
    end
    
    % Find the first index closest to 100 K
    [~, Index100] = min(abs(deltaTResults - 100));
    if isempty(Index100) 
        Index100 = Index300;
    end
    
        % Find the first index closest to 100 K
    [~, Index10] = min(abs(deltaTResults - 10));
    if isempty(Index100) 
        Index10 = Index100;
    end

    % Repackage outputs into struct
    LimitedOscillatorProperties.MaterialName = materialName;
    LimitedOscillatorProperties.Frequency = freqResults;
    LimitedOscillatorProperties.DeltaT = deltaTResults;
    LimitedOscillatorProperties.DeltaV = deltaVResults;
    LimitedOscillatorProperties.Capacitance = CapResults;
    LimitedOscillatorProperties.MinT = MinTResults;
    LimitedOscillatorProperties.MaxT = MaxTResults;
    LimitedOscillatorProperties.MinV = MinVResults;
    LimitedOscillatorProperties.MaxV = MaxVResults;
    LimitedOscillatorProperties.TInds = [Index10, Index100, Index300, Index500, Index900];
    LimitedOscillatorProperties.Power = Power(1);

end
%%
function results = generatePlots(OscProperties, hex_codes, designConstraints,FrequencyorCapacitanceConstraint)
    results = cell(length(OscProperties)+1,4);
    
    if FrequencyorCapacitanceConstraint == 0 % capacitance constraint 
        results{1,1} = "Material Name";
        results{1,2}="Frequency (Hz)";
        results{1,3} ="Power (W)";
        results{1,4} = "DeltaV / DeltaT (V/K)";

        fig_width = 6.5;   % Width of the figure in inches
        fig = figure;
        fig.Position(3) = fig_width * 100;  % inches to pixels for figure width
        
        left_subplot1_pos = [0.1 0.56 0.35 0.4];   
        left_subplot2_pos = [0.1 0.12 0.35 0.4];    
        
        right_subplot1_pos = [0.56 0.75 0.42 0.22]; 
        right_subplot2_pos = [0.56 0.43 0.42 0.22]; 
        right_subplot3_pos = [0.56 0.10 0.42 0.22]; 
        
        
        subplot('Position', left_subplot1_pos);
        
        for j = 1:length(OscProperties)
            plot(OscProperties{j,1}.Capacitance, OscProperties{j,1}.Frequency,'-o','markersize',2, 'Color',string(hex_codes(j,2)))
            hold on 
            % Plot extrapolations 
            plot(OscProperties{j,1}.extrapolatedcap, OscProperties{j,1}.extrapolatedfreq, '-.','Color',string(hex_codes(j,2)), 'LineWidth',0.5);
            
            % Unpackage cutoff temperature indices
            TInds = OscProperties{j,1}.TInds;
            
            % Plot first index above 100 K with triangle pointing up
            scatter(OscProperties{j,1}.Capacitance(TInds(2)),OscProperties{j,1}.Frequency(TInds(2)),40,...
                '^','MarkerFaceColor',string(hex_codes(j,2)),'MarkerEdgeColor','none');
            % Plot first index above 300 K with triangle pointing down
            scatter(OscProperties{j,1}.Capacitance(TInds(3)),OscProperties{j,1}.Frequency(TInds(3)),40,...
                'v','MarkerFaceColor',string(hex_codes(j,2)),'MarkerEdgeColor','none');
            % Plot first index above 500 K with square 
            scatter(OscProperties{j,1}.Capacitance(TInds(4)),OscProperties{j,1}.Frequency(TInds(4)),40,...
                's','MarkerFaceColor',string(hex_codes(j,2)),'MarkerEdgeColor','none');
        end
    
        % Plot vertical line intersecting design capacitance
        xline(designConstraints.Capacitance, '--k','LineWidth',1.2) 
        
        ylabel('Frequency (Hz)');
        set(gca, 'YScale', 'log'); 
        set(gca, 'XScale', 'log'); 
    
        xlim([1e-12 1e-8]);
        ylim ([1e0 1e8]);
        yticks([1e0 1e2 1e4 1e6 1e8]);
        xticklabels([]);
        set(gca,'Color', 'none');
        set(gca,'FontName','Arial','FontSize',10,'LineWidth',1);
    
        % DeltaV / DeltaT vs Capacitance 
        subplot('Position', left_subplot2_pos);
        
        for j = 1:length(OscProperties)
        
            dVdT = OscProperties{j,1}.DeltaV./OscProperties{j,1}.DeltaT;
            plot(OscProperties{j,1}.Capacitance,dVdT,'-o','markersize',2,'Color',string(hex_codes(j,2)));
            hold on 
            % Unpackage temperature cutoff indices
            TInds = OscProperties{j,1}.TInds;
        
            % Plot first index above 100 K 
            scatter(OscProperties{j,1}.Capacitance(TInds(2)),dVdT(TInds(2)),40,'^','MarkerFaceColor',string(hex_codes(j,2)),'MarkerEdgeColor','none');
            % Plot first index above 300 K 
            scatter(OscProperties{j,1}.Capacitance(TInds(3)),dVdT(TInds(3)),40,'v','MarkerFaceColor',string(hex_codes(j,2)),'MarkerEdgeColor','none');
            % Plot first index above 500 K 
            scatter(OscProperties{j,1}.Capacitance(TInds(4)),dVdT(TInds(4)),40,'s','MarkerFaceColor',string(hex_codes(j,2)),'MarkerEdgeColor','none');
            % plot extrapolated lines
            plot(OscProperties{j,1}.extrapolatedcap, OscProperties{j,1}.extrapolatedVT, '-.','markersize',5, 'Color',string(hex_codes(j,2)), 'LineWidth',0.5);
        end
    
            xline(designConstraints.Capacitance, '--k','LineWidth',1.2) % capacitance constraint line 
            set(gca, 'YScale', 'log'); 
            set(gca, 'XScale', 'log'); 
            ylabel('\DeltaV/\DeltaT (V/K)');
            xlabel('Circuit Capacitance (F)');
            xlim([1e-12 1e-8]);
            ylim ([1e-4 1e-1]);
            % yticks([1e-4 1e-3 1e-2]);
            set(gca,'Color', 'none');
            set(gca,'FontName','Arial','FontSize',10,'LineWidth',1);
            
          
        %% Scatter Plot frequency vs power for all materials that meet the design capacitance
        subplot('Position', right_subplot1_pos);
        for j = 1:length(OscProperties)
            hold on 
            % only plot if the extrapolated line extends to design
            % capacitance
            if length(OscProperties{j,1}.extrapolatedfreq) == 3
                plot(OscProperties{j,1}.Power,OscProperties{j,1}.extrapolatedfreq(2),'o','markersize',5,'MarkerFaceColor',string(hex_codes(j,2)),'MarkerEdgeColor','none');
                % Add text to plot
                text(OscProperties{j,1}.Power(1), OscProperties{j,1}.extrapolatedfreq(2), hex_codes(j,1), 'verticalAlignment', 'bottom', 'HorizontalAlignment', 'right','FontSize',10,'Color',string(hex_codes(j,2)));
               
                % Add power and frequency to an array
                results{j+1,1} = OscProperties{j,1}.MaterialName;
                results{j+1,2} = OscProperties{j,1}.extrapolatedfreq(2);
                results{j+1,3} = OscProperties{j,1}.Power;

            % check if other lines from simulations already reach design capacitance, and interpolate value 
            elseif  isempty(OscProperties{j,1}.interpolatedfreq) == 0
                plot(OscProperties{j,1}.Power,OscProperties{j,1}.interpolatedfreq,'o','markersize',5,'MarkerFaceColor',string(hex_codes(j,2)),'MarkerEdgeColor','none');
                % Add text to plot
                text(OscProperties{j,1}.Power(1), OscProperties{j,1}.interpolatedfreq, hex_codes(j,1), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right','FontSize',10,'Color',string(hex_codes(j,2)));
                
                results{j+1,1} = OscProperties{j,1}.MaterialName;
                results{j+1,2} = OscProperties{j,1}.interpolatedfreq;
                results{j+1,3} = OscProperties{j,1}.Power;
            end
        end
    
        set(gca, 'YScale', 'log'); 
        set(gca, 'XScale', 'log'); 
        ylabel('Frequency (Hz)');
        xlabel('Power (W)');
        xlim ([1e-5 1.5e-2]); 
        xticks([1e-5 1e-4 1e-3 1e-2]);
        ylim([1e0 1e8]);
        yticks([1e0 1e4 1e8]);
        ax = gca;
        ax.YColor='k'; 
        ax.Box = 'on';
        set(gca,'Color', 'none');
        set(gca,'FontName','Arial','FontSize',10,'LineWidth',1);
        
        
        %% DeltaV/DeltaT vs power at design capacitance
        subplot('Position', right_subplot2_pos);
        for j = 1:length(OscProperties)
            hold on 
            % only plot if the extrapolated line extends to the design capacitance
            if length(OscProperties{j,1}.extrapolatedfreq) == 3
                plot(OscProperties{j,1}.Power,OscProperties{j,1}.extrapolatedVT(2),'o','markersize',5,'MarkerFaceColor',string(hex_codes(j,2)),'MarkerEdgeColor','none');
                % Add text to plot
                text(OscProperties{j,1}.Power(1), OscProperties{j,1}.extrapolatedVT(2), hex_codes(j,1), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right','FontSize',10,'Color',string(hex_codes(j,2)));
                results{j+1,4} = OscProperties{j,1}.extrapolatedVT(2);

            % check if other lines from simulations already reach design capacitance, and interpolate value 
            elseif  isempty(OscProperties{j,1}.interpolatedfreq) == 0
                plot(OscProperties{j,1}.Power,OscProperties{j,1}.interpolatedVT,'o','markersize',5,'MarkerFaceColor',string(hex_codes(j,2)),'MarkerEdgeColor','none');
                % Add text to plot
                text(OscProperties{j,1}.Power, OscProperties{j,1}.interpolatedVT, hex_codes(j,1), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right','FontSize',10,'Color',string(hex_codes(j,2)))
                results{j+1,4} = OscProperties{j,1}.interpolatedVT;
            end
        end
        
        set(gca, 'YScale', 'log'); 
        set(gca, 'XScale', 'log'); 
        xlabel('Power (W)');
        ylabel('\DeltaV/\DeltaT (V/K)');
        xlim ([1e-5 1.5e-2]); 
        xticks([1e-5 1e-4 1e-3 1e-2]);
        ylim([1e-6 1e-2]);
        yticks([1e-6 1e-4 1e-2]);
        ax=gca
        ax.Box='on';
        set(gca,'Color', 'none');
        set(gca,'FontName','Arial','FontSize',10,'LineWidth',1);
     
        
         %% DeltaV/DeltaT vs Frequency at design capacitance 
        subplot('Position', right_subplot3_pos);
        
        for j = 1:length(OscProperties)
            hold on 
            % only plot if the extrapolated line extends to design capacitance
            if length(OscProperties{j,1}.extrapolatedfreq) == 3
                plot(OscProperties{j,1}.extrapolatedfreq(2),OscProperties{j,1}.extrapolatedVT(2),'o','markersize',5,'MarkerFaceColor',string(hex_codes(j,2)),'MarkerEdgeColor','none');
                % Add text to plot
                text(OscProperties{j,1}.extrapolatedfreq(2), OscProperties{j,1}.extrapolatedVT(2), hex_codes(j,1), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right','FontSize',10,'Color',string(hex_codes(j,2)));
            % check if other lines from simulations already reach design capacitance, and interpolate value 
            elseif  isempty(OscProperties{j,1}.interpolatedfreq) == 0
                plot(OscProperties{j,1}.interpolatedfreq,OscProperties{j,1}.interpolatedVT,'o','markersize',5,'MarkerFaceColor',string(hex_codes(j,2)),'MarkerEdgeColor','none');
                % Add text to plot
                text(OscProperties{j,1}.interpolatedfreq, OscProperties{j,1}.interpolatedVT, hex_codes(j,1), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right','FontSize',10,'Color',string(hex_codes(j,2)));
            end
        end
    
        set(gca, 'YScale', 'log'); 
        set(gca, 'XScale', 'log'); 
        xlabel('Frequency (Hz)');
        ylabel('\DeltaV/\DeltaT (V/K)');
        xlim ([1e0 1e8]); 
        xticks([1e0 1e2 1e4 1e6 1e8]);
        set(gca,'Color', 'none');
        ax=gca;
        ax.Box='on';
        set(gca,'FontName','Arial','FontSize',10,'LineWidth',1);
    
    else  % frequency constraint 
        results{1,1} = "Material Name";
        results{1,2}="Capacitance (F)";
        results{1,3} ="Power (W)";
        results{1,4} = "DeltaV / DeltaT (V/K)";
        % Frequency vs Capacitance
        fig_width = 6.5;   % Width of the figure in inches
        fig = figure;
        fig.Position(3) = fig_width * 100;  % inches to pixels for figure width
        
        left_subplot1_pos = [0.1 0.26 0.4 0.5];      
        
        right_subplot1_pos = [0.6 0.7 0.3 0.22]; 
        right_subplot2_pos = [0.6 0.4 0.3 0.22]; 
        right_subplot3_pos = [0.6 0.12 0.3 0.22]; 
        
        subplot('Position', left_subplot1_pos);
        for j = 1:length(OscProperties)
            plot(OscProperties{j,1}.Capacitance, OscProperties{j,1}.Frequency,'-o','markersize',2, 'Color',string(hex_codes(j,2)))
            hold on 
            % Plot extrapolations 
            plot(OscProperties{j,1}.extrapolatedcap, OscProperties{j,1}.extrapolatedfreq, '--','Color',string(hex_codes(j,2)), 'LineWidth',1.2);
            
            % Unpackage cutoff temperature indices
            TInds = OscProperties{j,1}.TInds;
            
            % Plot first index above 10 K with hexagram pointing up
            scatter(OscProperties{j,1}.Capacitance(TInds(1)),OscProperties{j,1}.Frequency(TInds(1)),70,...
                'hexagram','MarkerFaceColor',string(hex_codes(j,2)),'MarkerEdgeColor','none');
            % Plot first index above 100 K with triangle pointing up
            scatter(OscProperties{j,1}.Capacitance(TInds(2)),OscProperties{j,1}.Frequency(TInds(2)),70,...
                '^','MarkerFaceColor',string(hex_codes(j,2)),'MarkerEdgeColor','none');
            % Plot first index above 300 K with triangle pointing down
            scatter(OscProperties{j,1}.Capacitance(TInds(3)),OscProperties{j,1}.Frequency(TInds(3)),70,...
                'v','MarkerFaceColor',string(hex_codes(j,2)),'MarkerEdgeColor','none');
            % Plot first index above 500 K with square 
            scatter(OscProperties{j,1}.Capacitance(TInds(4)),OscProperties{j,1}.Frequency(TInds(4)),70,...
                's','MarkerFaceColor',string(hex_codes(j,2)),'MarkerEdgeColor','none');
            % % Add text to plot
        end
    
        % Plot vertical line intersecting design capacitance
         yline(designConstraints.Frequency, '--k','LineWidth',2.5) % frequency constraint line 
        set(gca,'FontName','Arial','FontSize',10,'LineWidth',1);
        xlabel('Circuit Capacitance (F)');
        ylabel('Frequency (Hz)');
        set(gca, 'YScale', 'log'); 
        set(gca, 'XScale', 'log'); 
        set(gca, 'FontName', 'Arial','LineWidth',1.2);
        set(gca,'Color', 'none');
        xlim([1e-12 1e-8])
        ylim ([1e0 1e10]);
    
        %% Frequency vs Capacitance at design frequency with no extrapolations 
        subplot('Position',right_subplot1_pos);
        yline(designConstraints.Frequency, '--k','LineWidth',1) % frequency constraint line 
        hold on 
        for j = 1:length(OscProperties)
            % if the design frequency is under the critical
            % capacitance's frequency (maximum frequency) and the capacitance is 
            % between 1e-8 1e-12 F (realistic capacitances), plot the point
            if designConstraints.Frequency < OscProperties{j,1}.Frequency(1)...
                    && OscProperties{j,1}.freqdesigncap > 10^-12 &&  OscProperties{j,1}.freqdesigncap < 10^-8
              scatter(OscProperties{j,1}.freqdesigncap, designConstraints.Frequency, 50, 'MarkerFaceColor',string(hex_codes(j,2)),'MarkerEdgeColor','none')
              hold on 
              text(OscProperties{j,1}.freqdesigncap, designConstraints.Frequency, hex_codes(j,1), 'verticalAlignment', 'bottom', 'HorizontalAlignment', 'right','FontSize',10,'Color',string(hex_codes(j,2)));
             
              results{j+1,1} = OscProperties{j,1}.MaterialName;
              results{j+1,2} = OscProperties{j,1}.freqdesigncap;
            end
        end
    
        set(gca,'FontName','Arial','FontSize',10,'LineWidth',1);
        ylabel('Frequency (Hz)');
        set(gca, 'YScale', 'log'); 
        set(gca, 'XScale', 'log'); 
        set(gca, 'FontName', 'Arial','LineWidth',1.2);
        set(gca,'Color', 'none');
        box on 
        xlim([1e-12 1e-8])
        xticklabels([]);
        
        
        %% dvdt vs Capacitance at design frequency 
        subplot('Position',right_subplot2_pos);
        for j = 1:length(OscProperties)
            % if the design frequency is under the critical
            % capacitance's frequency (maximum frequency), plot the point
            if designConstraints.Frequency < OscProperties{j,1}.Frequency(1) && OscProperties{j,1}.freqdesigncap > 10^-12
              scatter(OscProperties{j,1}.freqdesigncap, OscProperties{j,1}.freqdesignVT, 50, 'MarkerFaceColor',string(hex_codes(j,2)),'MarkerEdgeColor','none')
              hold on 
            results{j+1,4} = OscProperties{j,1}.freqdesignVT;

            end
        end
    
        set(gca,'FontName','Arial','FontSize',10,'LineWidth',1);
        ylabel('\DeltaV/\DeltaT (V/K)')
        set(gca, 'YScale', 'log'); 
        set(gca, 'XScale', 'log'); 
        set(gca, 'FontName', 'Arial','LineWidth',1.2);
        set(gca,'Color', 'none');
        box on 
        xlim([1e-12 1e-8])
        xticklabels([]);
        
    
        %% Power vs Capacitance at design frequency with no extrapolations 
        subplot('Position',right_subplot3_pos);
        for j = 1:length(OscProperties)
            % if the design frequency is under the critical
            % capacitance's frequency (maximum frequency), plot the point
            if designConstraints.Frequency < OscProperties{j,1}.Frequency(1) && OscProperties{j,1}.freqdesigncap > 10^-12
              scatter(OscProperties{j,1}.freqdesigncap, OscProperties{j,1}.Power, 50, 'MarkerFaceColor',string(hex_codes(j,2)),'MarkerEdgeColor','none')
              hold on 
              results{j+1,3} = OscProperties{j,1}.Power;

            end
        end
        set(gca,'FontName','Arial','FontSize',10,'LineWidth',1);
        xlabel('Circuit Capacitance (F)');
        ylabel('Power (W)')
        set(gca, 'YScale', 'log'); 
        set(gca, 'XScale', 'log'); 
        set(gca, 'FontName', 'Arial','LineWidth',1.2);
        set(gca,'Color', 'none');
        box on 
        xlim([1e-12 1e-8])
    end

    filteredData = {};
    % Remove any empty cells from results 
    for j = 1:length(results)
         if all(~cellfun(@isempty, results(j, :)))
              filteredData = [filteredData; results(j, :)]; 
         end
    end
    results = filteredData;
end