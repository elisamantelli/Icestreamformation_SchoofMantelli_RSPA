N = 8;
maxNumCompThreads(N)

%set physical parameters
parameters.T_s = -1;                                                       % surface temperature
parameters.nu = 0.5;                                                       % geothermal heat flux
parameters.alpha = 1;                                                      % Brinkman number
parameters.gamma = 0.1;                                                    % friction coefficient
parameters.mu = 0.1;                                                       % yield strength in Coulomb friction
parameters.Pe = 1;                                                         % Peclet number
parameters.a = 1;                                                          % accumulation rate
parameters.r = 917/1000;                                                   
parameters.beta = 0.1;                                                     % bed permeability
parameters.T_bed = -0.2;                                                   % basal temperature for slab-type simulations. Irrelevant for standard set up                                                

parameters.flag.heat_full = 1;
parameters.flag.plug = 0;
parameters.flag.pert = 'rand_gamma';                                       % set style of perturbation. 'rand_gamma' perturbs friction coefficient; 'rand' perturbs basal temperature; 'mono' perturbs basal temperature with a monochromatic wave
parameters.flag.hydrology = 'weertman';                                    % set type of friction. Other choices are: 'budd', 'coulomb'
parameters.flag.fluxy_sigma = 'Pionly';                                    % set controls on transverse drainage flux. Keep as is.
parameters.flag.hydropermnew = 1;                                          

%regularization parameters 
%friction coefficient
parameters.reg.epsilon = 0.005;                                            % cutoff for regularised heaviside function. Must be smaller than epsilon_f
parameters.reg.epsilon_f = 0.03;                                           % This sets the temperature dependence of friction. Called '\delta' in paper

 %permeability
parameters.reg.c = 3;                                                      % Exponent in the permeability power law. Must be strictly larger than 1. 
parameters.reg.Pi_0 = 0.1;                                                 % Regularizer for Coulomb friction law. Sets the derivative of the regularization for Pi=0. Must be small
parameters.reg.Pi_s = 0.1;                                                 % Regularizer for Buss friction law. Sets the derivative of the regularization for Pi=0. Must be small

%bed elevation
parameters. bed. b0 = 0;                                                   % Bed elevation at the divide
parameters. bed. b1 = 0.05;                                                % Bed slope, must be strictly positive

%% timestepping
parameters.timestep.x_init = 0.00036;                                      % x coordinate where the simulation is started
parameters.timestep.dx = 0.00036;                                          % dx at next time step, i+1. Must remain comparable with/smaller than epsilon_f. Set at least dx=parameters.reg.epsilon_f/10
parameters.timestep.dx_prev =  0.00036;                                    % dx at previous timestep, i
parameters.timestep.dx_pprev =  0.00036;                                   % dx two time steps before current, i-1
parameters.timestep.x = 0.00036;                                           % current x coordinate, x(i+1)
parameters.tplot = parameters.timestep.dx;                                 % output writing intervals
parameters.stepmax = 3125000;                                              % max # of timesteps

parameters.da = 0;                                                         %parameters.timestep.dx;

%% grid generation
parameters.n_nodes_transverse = 1;                                         %#nodes in the transverse direction in the uppermost layer. Remaining grid parameters are set in the grid routine. Must be even!
parameters.ratio_hor = 5;                                                  % This sets the cell aspect ratio, dy = dz*ratio_hor
parameters.ratio_vert = 1;                                                 % This sets the number of nodes along the vertical direction as n_nodes = ratio_ver*[4 4 4 8]
parameters.flag1d = 1;                                                     % set to 1 to have 1 cell laterally, set to 0 to have many cells laterally (y)
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






