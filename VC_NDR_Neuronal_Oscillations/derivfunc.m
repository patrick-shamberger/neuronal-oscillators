function [deriv,state] = derivfunc(t,vec_x,state,params)

% Function "derivfunc.m" v3.0.0, tested 29 June 2021
% Written by T.D. Brown Fall 2019
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
    L = params.L;
    T0 = params.T0;
    V0 = params.V0;
    gt0 = params.GthofT0;
    g0 = params.GelofT0;
    c0 = params.CthofT0;
    %gg = params.GelofTss;
    %rhoSS = params.rhoSS;

    funpars = params.funpars;

%% Define the Differential Equation System

    % Defines x'=f1(x,y,z) y'=f2(x,y,z), z'=f3(x,y,z)
    % x = temperature, y = cap 1 voltage, z = cap 2 voltage
    % To-Do: Vectorize for speed (Is this even possible with non-linearity??) 

    % The below implementation is in SI units, but can easily have
    % convergence issues
    % vec_x = [Temperature, Current] (in rows) 
    %{
    deriv1 = (1/Cth)*(vec_x(2,:).^2/Econd - (vec_x(1,:)-T0).*Tcond);
    deriv2 = -(1/L)*(vec_x(2,:)/Econd - V0);
    deriv = [deriv1; deriv2];
    %}

    % This implementation can be obtained from the first with the formal
    % change of variables t -> t/t_th with t_th = 1/(2pi*Cth*Rth(T0)) and
    % T -> T/T0 and v-> I/I0 (I0 = sqrt(T0*Gth(T0)*Gel(T0) )
    % It's a little more complicated, but is much nicer to make sure the
    % simulation converges
    %
    
    % Tough part first: Evaluate both noninear electrical and thermal
    % conductance functions
    ET = electrothermalmodels(funpars, vec_x(1,:)' *T0);
    Econd = ET.Econd;
    Tcond = ET.Tcond;
    Cth = ET.Tcap;

    normLC = c0/L;
    scaleC = Cth/c0;
    normV = V0/sqrt(T0*gt0/g0);
    normg = Econd/g0;
    normgt = Tcond/gt0;

    % Series inductor + electrothermal memristor
    deriv1 = 2*pi*(1/scaleC)*(vec_x(2,:).^2./normg - (vec_x(1,:)-1).*normgt);
    deriv2 = -2*pi*normLC/(g0*gt0)*(vec_x(2,:)./normg - normV);
    deriv = [deriv1; deriv2];     

end
   

% Another formula for normI needed for extended memristors (there may be 
% multiple Vss for each Tss):
% normI = (Tss/T0-1)*Gth(Tss)/Gth(T0)*V0/Vss