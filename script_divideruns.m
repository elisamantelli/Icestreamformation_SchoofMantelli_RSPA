N = 8;
maxNumCompThreads(N)

%% set physical parameters
parameters.T_s = -1; 
parameters.a = 1;
parameters.r = 917/1000;
parameters.beta = 1;
parameters.mu = 1;
%% bed geometry
parameters. bed. b0 = 0; %elevation at the divide, x= 0
parameters. bed. b1 = 0;   %0.15; %bed slope, must be strictly positive

%% physics & perturbation
parameters.flag.heat_full = 1;
parameters.flag.plug = 0;
parameters.flag_Tdep = 1; %turn on (1) and off(0) T_dep sliding
parameters.flag.hydrology = 'budd'; %'weertman', 'budd'

parameters.flag.pert = 'rand'; %'mono', 'gamma_rand'
parameters.amplitude = 0;    %amplitude of bed temperature perturbation; set to zero to compute steady states
parameters.n_wavelengths = 1; %number of wavelengths within assigneddomain size
parameters.gamma_pert = 0;
%% regularization parameters 
%friction coefficient
parameters.reg.epsilon = 0.005; %cutoff for regularised heaviside function. Should be smaller than epsilon_f
parameters.reg.epsilon_f = 0.03; %delta. keep O(1) to start with

%permeability
parameters.reg.c = 1.2; %keep as strictly larger than 1. Exponent in the permeability power law
parameters.reg.Pi_0 = 1; %tthis sets the derivative of the regularization for Pi=0. Eventually should be small. 
parameters.reg.Pi_s = 1;

%% timestepping
parameters.timestep.x_init = 0.001;
parameters.timestep.dx = 0.001; %(i+1)
parameters.timestep.dx_prev =  0.001; %i
parameters.timestep.dx_pprev =  0.001; %i-1
parameters.timestep.x = parameters.timestep.x_init; % update at each timestep
parameters.tplot = parameters.timestep.dx;
parameters.stepmax = 10000;
parameters.iter.x_current = parameters.timestep.x_init;

parameters.da = 0;    %parameters.timestep.dx;

%% grid generation
parameters.n_nodes_transverse = 1;             %#nodes in the transverse direction in the uppermost layer. Remaining grid parameters are set in the grid routine. Must be even!
parameters.ratio_hor = 5;  %dy = dz*ratio
parameters.ratio_vert = 1;
parameters.flag1d = 1;
if parameters.flag1d == 0
    parameters.grid = fv_grid_transverse_v6(parameters);
elseif parameters.flag1d == 1
    parameters.grid = fv_grid_transverse_1d_v4(parameters);
end

%% parameter sweep
h_list = [0.8 1 1.5 2];
gamma_list = [0.1 0.5 1 2];
nu_list = [0.1 0.5 1 2];
Pe_list = [1 10];
alpha_list = [1 5];

for k = 1: length(h_list)
    parameters.h_init = h_list(k); 
    parameters.h = parameters.h_init;
    
    for g = 1:length(gamma_list)
        parameters.gamma = gamma_list(g);
        
        for n = 1:length(nu_list)
            parameters.nu = nu_list(n);
            
            for a = 1:length(alpha_list)
                parameters.alpha = alpha_list(a);
                
                for p = 1:length(Pe_list)
                    parameters.Pe = Pe_list(p);
                    timestepping_divide(parameters)
                end
            end
            
        end
        
    end
    
   
end






