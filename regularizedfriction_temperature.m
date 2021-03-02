function [fout, Dfout] = regularizedfriction_temperature(T_bed, parameters)
%For a friction law of the form 
%tau_b = (gamma.*(1+gamma_pert).*f_slide_Tbed.*.*f_Pi.*u_bed)
%regularizedfriction_temperature provides f_slide_Tbed(T_bed) and its
%derivative

epsilon = parameters.reg.epsilon; 
epsilon_f = parameters.reg.epsilon_f;

T_bed = full(T_bed);

H0 =  1;
H01 = -0.1;
x0 = 0;

H = @(x)  1/2*(1 + erf(pi/2.*(x- H0-x0)/epsilon ));
dH = @(x) sqrt(pi)/(2*epsilon).*(exp(-(pi/2.*(x- H0-x0)./epsilon).^2));  %(epsilon+epsilon*cosh(2*x/epsilon)).^(-1);

H1 = @(x)  1/2*(1 + erf(pi/2.*(x-H01-x0)/epsilon ));
dH1 = @(x) sqrt(pi)/(2*epsilon).*(exp(-(pi/2.*(x-H01-x0)./epsilon ).^2));  %(epsilon+epsilon*cosh(2*x/epsilon)).^(-1);

f = @(x) exp(-x/epsilon_f);
df = @(x) -1/epsilon_f*exp(-x/epsilon_f);

f_high =  H1(T_bed) + (1-H(T_bed)).*f(T_bed); %1+ H_ev.*(f(T_bed)-1);
df_high = dH1(T_bed)+ df(T_bed).*(1-H(T_bed)) - dH(T_bed).*f(T_bed);
 
x0 = -10*parameters.reg.epsilon_f;

if epsilon_f == 0.03
    x01 = 0;  
    x02 = 0.0003;  
elseif epsilon_f == 0.01
    x01 = 0.00001;  
    x02 = 0.00001;
elseif epsilon_f == 0.005
    x01 = 0.00005 ;
    x02 = 0.00001;
elseif epsilon_f == 0.001
    x01 = 0.0001;
    x02 = 0.0001;
elseif epsilon_f == 0.0005
    x01 = -01e-05;
    x02 = 0.0001;
elseif epsilon_f == 0.00005
    x01 = -01e-04;
    x02 = 0.0001;
elseif epsilon_f == 0.000005
    x01 = -07e-05;
    x02 = 0.00005;
end
threshold = H1(x0) + (1-H(x0)).*f(x0);

H = @(x)  1/2*(1 + erf(pi/2.*(x + x0 + x02)/(epsilon) ));
dH = @(x) sqrt(pi)/(2*(epsilon)).*(exp(-(pi/2.*(x + x0 + x02)./(epsilon)).^2));  
H1 = @(x)  1/2*(1 + erf(pi/2.*(x+x0 -x01)/(epsilon) ));
dH1 = @(x) sqrt(pi)/(2*(epsilon)).*(exp(-(pi/2.*(x+x0-x01)./(epsilon) ).^2)); 

fout = threshold.*H1(-T_bed) + (1-H(-T_bed)).*f_high;
Dfout =  threshold.*(-dH1(-T_bed)) + (1-H(-T_bed)).*df_high +dH(-T_bed).*f_high;

fout(isnan(fout)) = threshold;
Dfout(isnan(Dfout)) = 0;
