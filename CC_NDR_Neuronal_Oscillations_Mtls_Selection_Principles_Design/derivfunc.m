function [deriv,state] = derivfunc(t,vec_x,state,params)

% Function "derivfunc.m" v3.0.0, tested 18 September 2024
% Written by T.D. Brown Fall 2023
%
% This function defines the particular differential equation to be
% integrated within RK4Nsolver.m, and is totally user-defined. However, a
% good general work flow appears to be to import a params struct from
% main.m for parametric analyses (and a state struct when necessary). The
% user simply needs to define how the output of each state variable
% derivative depends on the current value of the state variables.
%
% For example, the differential equation system x'=y, y'=-by-kx can be
% implemented as
% deriv = [0 1; -k -b]*vec_x, for some parameters k and b
%
% There is also the option to solve forced or non-homogeneous differential
% equation systems using functions of time t, for example
% deriv = [0 1; -k -b]*vec_x + [force1; force2]
% force1 = cos(t);
% force2 = sin(4*t);
% 
%
% Input: t = time, for inhomogeneous forcing functions of time
% Input: vec_x = vector of system variables at times t, [ARRAY, vert: vars, horiz: times]
% Input: state = struct of local state variables, for example, deriv at last timestep
% Input: params = struct of params defining differential equation system
% Output: deriv = vector of system derivatives at time t, [ARRAY, vert: derivatives, horiz: times]
% Output: state = struct of updated local state variables

%% Import internal and external parameters from main.m

    state.blank = []; % Not currently using state-dependence of solver
    Cel = params.Cel; % Circuit capacitance [F]
    T0 = params.T0; % Ambient temperature [K]
    I0 = params.I0; % Bias current [A]
    gt0 = params.GthofT0; % T Conductance at T0 [W/K] 
    g0 = params.GelofT0; % E Conductance at T0 [S]
    c0 = params.CthofT0; % T Capacitance at T0 [J/K]
    gg = params.GelofTss; % E Conductance at Tss
    rhoSS = params.rhoSS; % Limiter resistance as fraction of R(Tss)

    funpars = params.funpars; % Import material model parameters

%% Define the Differential Equation System

    % Defines x'=f1(x,y) y'=f2(x,y), x = temperature, y = voltage
 
    % Diff eqs: 
    % Thermal dynamics (Joule heating - Newton cooling)
    % Voltage dynamics (Capacitor charging)

    % The below implementation is in SI units, but can easily have
    % convergence issues
    %{
    deriv1 = (1/Cel)*(Econd.*vec_x(2,:).^2 - (vec_x(1,:)-T0).*Tcond);
    deriv2 = -(1/Cth)*(Econd.*vec_x(2,:) - I0);
    deriv = [deriv1; deriv2];
    %}

    % The implementation below comes from the first, with the formal
    % change of variables t -> t/t_th with t_th = (2pi*Cth(T0)*Rth(T0)) and
    % T -> T/T0 and v-> V/V0 (V0 = sqrt(T0*Gth(T0)/Gel(T0) )
    % It's a little more complicated, but is much nicer to make sure the
    % simulation converges. See SI.
    
    % Evaluate T-dependent electrical and thermal conductance / capacitance
    ET = electrothermalmodels(funpars, vec_x(1,:)' *T0);
    Econd = ET.Econd;
    Tcond = ET.Tcond;
    Cth = ET.Tcap;

    % Renormalized variables
    normC = c0/Cel;
    scaleC = Cth/c0;
    normI = I0/sqrt(T0*gt0*g0);
    normg = Econd/g0;
    normgt = Tcond/gt0;

    % Include simulation with limiter resistor, if desired
    rho = 0;
    if rhoSS > 0
        rho = Econd./(Econd + gg.*(1/rhoSS-1));
    end

    % Return the derivatives dT/dt and dv/dt at the specified temperature
    % (Renormalized variables; see SI)
    deriv1 = 2*pi*(1/scaleC)*(normg.*vec_x(2,:).^2.*(1-rho).^2 - (vec_x(1,:)-1).*normgt);
    deriv2 = -2*pi*normC*g0/gt0*(normg.*vec_x(2,:).*(1-rho) - normI);
    deriv = [deriv1; deriv2];     

end