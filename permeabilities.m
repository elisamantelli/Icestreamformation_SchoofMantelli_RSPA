function [fout, Dfout] = permeabilities(x, parameters)

c= parameters.reg.c;
Pi_0 = parameters.reg.Pi_0;

if isfield(parameters.flag, 'hydropermnew')==0
    kappa =  (x.^2 + Pi_0^2).^((c-1)/2).*x;
    dkappa = (x.^2 + Pi_0^2).^((c-1)/2) + 1/2*(c-1).*(Pi_0^2 + x.^2).^(-1+1/2*(c-1)).*x.^2.*2;
elseif isfield(parameters.flag, 'hydropermnew')==1
    kappa =  (x.^2 + Pi_0^2).^((c)/2).*x;
    dkappa = (x.^2 + Pi_0^2).^((c)/2) + 1/2*(c).*(Pi_0^2 + x.^2).^(-1+1/2*(c)).*x.^2.*2;    
end
%kappa(index_neg) = 0;
fout.kappa = kappa;
fout.kappa_2 = c*ones(length(x),1);

%dkappa(index_neg) = 0;
Dfout.kappa = dkappa;
Dfout.kappa2 = zeros(length(x),1);
end