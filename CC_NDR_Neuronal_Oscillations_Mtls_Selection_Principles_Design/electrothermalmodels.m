function out = electrothermalmodels(funpars,in)
% Function "electrothermalmodels.m" v3.0.0, tested 18 September 2024
% Written by T.D. Brown Fall 2023
%
% This function computes all electrothermal conductance and capacity
% functions at the desired temperatures "in" using the model parameters
% stored in "funpars", with calls to the modeling functions
% "PiecewiseArrhenian.m", "PiecewiseLaurent.m", and "DebyewithFOPT.m"
%
% Model parameters in "funpars" can be created with the right syntax using
% "modelfitter.m"
%
% Intrinsic functions (conductivity / specific heat) is automatically
% re-scaled to extrinsic functions (conductance, capacitance)
%
% Input: funpars = struct of all needed model parameters. 
% See "modelfitter.m", "PiecewiseArrhenian.m", "PiecewiseLaurent.m",
% and "DebyewithFOPT.m"
    % Epars = struct of electrical conductivity model parameters
    % Tpars = struct of thermal conductivity model parameters
    % Tcappars = struct of specific heat model parameters
    % dims = vector of dimensional factors 
% Input: in = vector of temperatures [K] for evaluating models
% Output: out = struct of modeled T-dependent electrothermal properties
    % Econd = T-dependent electrical conductance [S]
    % Tcond = T-dependent thermal conductance [W/K]
    % Tcap = T-dependent thermal capacitance [J/K]

    % Unpackage model parameters
    Epars = funpars.Epars;
    Tpars = funpars.Tpars;
    Tcappars = funpars.Tcappars;

    dims = funpars.dims;

    % Evaluate T-dependent electrical model
    % Input: T; Output: Gel(T) [S]
    % Implements piecewise Arrhenian function
    econd = PiecewiseArrhenian(0, Epars, in, 0);

    % Input: T; Output: Gth(T) [W/K]    
    % Implements a piecewise Laurent series model
    tcond = max(0,PiecewiseLaurent(0, Tpars, in, 0));

    % Input: T; Output: Cth(T) [J/K]    
    % Implements a piecewise Debye modelwith FOPT peak
    tcap = DebyewithFOPT(0, Tcappars, in);

    % Rescale parameters and output
    % Conductivity -> Conductance; Specific heat -> thermal capacitance
    out.Econd = econd*dims(1);
    out.Tcond = tcond*dims(2);
    out.Tcap = tcap*dims(3);

end