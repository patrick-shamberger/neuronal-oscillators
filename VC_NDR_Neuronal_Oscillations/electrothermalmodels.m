function out = electrothermalmodels(funpars,in)

    % Unpackage model parameters
    Epars = funpars.Epars;
    Tpars = funpars.Tpars;
    Tcappars = funpars.Tcappars;

    dims = funpars.dims;

    % Evaluate T-dependent electrical conductance model
    % Input: T; Output: Gel(T) [S]
    % This is modeled based on Alexandrov et al.,PRL 96, 117003 (2006) 
    % Total Conductivity of system (S/cm):
    econductivity = 1./TotalResistivity(Epars, in); % Total electrical Conductivity of system (S/cm):
    econductivity = econductivity*100; % Total electrical Conductivity of system (S/m):
    % Rescale parameters Conductivity -> Conductance;
    econd=econductivity*dims(1); %Electrical Conductance
   
    
   
    % Input: T; Output: Gth(T) [W/K]    
    % Implements a piecewise Laurent series model
    tcond = max(0,PiecewiseLaurent(0, Tpars, in, 0));
    tcond=tcond*dims(2);

    % Input: T; Output: Cth(T) [J/K]    
    % Implements a piecewise Debye modelwith FOPT peak
    tcap = DebyewithFOPT(0, Tcappars, in);
    tcap=tcap*dims(3);


    % Rescaled parameters and output
    out.Econd = econd;
    out.Tcond = tcond;
    out.Tcap = tcap;

end

function out = TotalResistivity(Epars,in)
    lowTpars = Epars.lowT;
    highTpars = Epars.highT;
    sigpars = Epars.sigmoid;
    lowTmodel = rho2(lowTpars,in);
    highTmodel = rho1(highTpars,in);
    s = volumefrac(sigpars,in); %(1+exp(-( (T-Tc)/sigc ))^-1 with x=T-Tc/sigC
% 
% 
    out = (highTmodel.^(1-s)).*(lowTmodel.^s);   
end

function out = rho2(lowTpars,in)
    %rho2= rho02+a*T.^2+b.*T.^4.5;     % Low temp Ferromagnetic phase
    rho02=lowTpars(1);
    a=lowTpars(2);
    b=lowTpars(3);
    
    out=rho02+a*in.^2+b.*in.^4.5;
end

function out = rho1(highTpars,in)
    %rho1= rho0.*T.*exp(Ep./(kB*T));       % High Temp Paramagnetic phase
    rho0=highTpars(1);
    Ep = highTpars(2);                          
    kB = highTpars(3);                           
    
    out= rho0.*in.*exp(Ep./(kB*in));
end

function out = volumefrac(sigpars,in)
    %s =(1/2).*erfc((T-TC)/Gamma);
    TC=sigpars(1);
    Gamma = sigpars(2);                          
    out= (1/2).*erfc((in-TC)/Gamma);
end

function out = shiftedreciprocal(params,in)
    offset = params(1);
    factor = params(2);
    out = offset+factor./in;
end

% function out = toyCthmodel(Tcappars,in)
%     out = Tcappars.val*ones(size(in));    
% end

