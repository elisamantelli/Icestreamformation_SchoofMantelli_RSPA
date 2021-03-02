function [fout, Dfout] = regularizedfriction_hydrology(ub,Pi, parameters)
%For a friction law of the form
%tau_b = (gamma.*(1+gamma_pert).*f_slide_Tbed.*f_Pi.*u_bed)
%regularizedfriction_hydrology.m provides f_Pi(ub,Pi) and its
%derivative for a variety of functional forms.


c= parameters.reg.c;
gamma_0 = parameters.gamma;
mu_0 = parameters.mu;
Pi_0 = parameters.reg.Pi_0;
Pi_s = parameters.reg.Pi_s;
if isfield( parameters.reg, 'u_0') == 1
    u_0 = parameters.reg.u_0;
else
    u_0 = 0.01;
end

index_p = find(Pi>=0);
index_n = Pi<0;

Pi_slide = zeros(length(Pi),1);
Pi_slide(index_p) = Pi(index_p);
Pi_slide(index_n) = 0;

dPislide_dPi = ones(length(Pi),1);
dPislide_dPi(index_n) = 0;

if strcmp(parameters.flag.hydrology, 'weertman')==1
    fout = ones(length(ub),1);
    Dfout.dPi = zeros(length(ub),1);
    Dfout.dub = zeros(length(ub),1);
    
elseif strcmp(parameters.flag.hydrology, 'coulomb')==1
    fout = 1./(gamma_0/mu_0.*sqrt(ub.^2 + u_0.^2).*(Pi_slide.^2 + Pi_0.^2).^((c-1)/2).*Pi_slide + 1);
    Dfout.dPi = -((Pi_slide.^2 + Pi_0.^2).^((c-1)/2).*(c.*Pi_slide.^2 + Pi_0.^2).*sqrt(ub.^2 + u_0.^2).*gamma_0.*mu_0./...
        (Pi_slide.*(Pi_slide.^2 + Pi_0.^2).^(c/2).*sqrt(ub.^2 + u_0.^2).*gamma_0 + (Pi_slide.^2 + Pi_0.^2).^(1/2).*mu_0 ).^2).*dPislide_dPi;
    Dfout.dub = - Pi_slide.*(Pi_slide.^2 + Pi_0.^2).^((c-1)/2).*(gamma_0/mu_0).*ub./...
        (sqrt(ub.^2 + u_0.^2).*(1+Pi_slide.*(Pi_slide.^2 + Pi_0.^2).^((c-1)/2).*sqrt(ub.^2 + u_0.^2).*gamma_0./mu_0 ).^2);
    
elseif strcmp(parameters.flag.hydrology, 'budd')==1
    %tau_b = gamma_0 u / (1 +(Pi/Pi_s)^c)
    if isfield(parameters.flag, 'hydropermnew')==0
        fout = 1./(1+(Pi_slide./Pi_s).^c);%1./((Pi.^2 + Pi_s.^2).^(c/2));
        Dfout.dub = zeros(length(ub),1);
        Dfout.dPi =-(c.*(Pi_s).^(-c).*Pi_slide.^(c-1))./(1+(Pi_slide./Pi_s).^c ).^2.*dPislide_dPi;   %-c.*Pi.*((Pi.^2 + Pi_s.^2).^(-1-c/2));
    elseif isfield(parameters.flag, 'hydropermnew')==1
        %tau_b = gamma_0 u / (1 +(Pi/Pi_s)^c)
        fout = 1./(1+ (Pi_slide.^2 + Pi_0.^2).^((c-1)/2).*Pi_slide./Pi_0.^c);%1./((Pi.^2 + Pi_s.^2).^(c/2));
        Dfout.dub = zeros(length(ub),1);
        Dfout.dPi = -(2*((c-1)/2).*(Pi_slide.^2 + Pi_0.^2).^((c-1)/2-1).*Pi_slide.^2./Pi_0.^c + (Pi_slide.^2 + Pi_0.^2).^((c-1)/2)./Pi_0.^c)  ./(1+ (Pi_slide.^2 + Pi_0.^2).^((c-1)/2).*Pi_slide./Pi_0.^c).^2 .*dPislide_dPi;   %-c.*Pi.*((Pi.^2 + Pi_s.^2).^(-1-c/2));
    end
end
end