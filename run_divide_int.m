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

parameters.timestep.x_init = 0.00036;                                      % x coordinate where the simulation is started
parameters.timestep.dx = 0.00036;                                          % dx at next time step, i+1. Must remain comparable with/smaller than epsilon_f. Set at least dx=parameters.reg.epsilon_f/10
parameters.timestep.dx_prev =  0.00036;                                    % dx at previous timestep, i
parameters.timestep.dx_pprev =  0.00036;                                   % dx two time steps before current, i-1
parameters.timestep.x = 0.00036;                                           % current x coordinate, x(i+1)
parameters.tplot = parameters.timestep.dx;                                 % output writing intervals
parameters.stepmax = 3125000;                                              % max # of timesteps

parameters.da = 0;                                                         %parameters.timestep.dx;


%grid generation
parameters.flag1d = 0;                                                     % set to 1 to have 1 cell laterally, set to 0 to have many cells laterally (y)
if parameters.flag1d == 1
    parameters.n_nodes_transverse = 1;
else
    parameters.n_nodes_transverse = 8;                                     % #nodes in the transverse direction in the uppermost layer. Remaining grid parameters are set in the grid routine
    if rem(parameters.n_nodes_transverse,2)>0
        warning('transverse cell number must be even')
    end
end
parameters.ratio_hor = 7;                                                  % This sets the cell aspect ratio, dy = dz*ratio_hor
parameters.ratio_vert = 4;                                                 % This sets the number of nodes along the vertical direction as n_nodes = ratio_ver*[4 4 4 8]

if parameters.flag1d == 0
    parameters.grid = fv_grid_transverse_v6(parameters);
elseif parameters.flag1d == 1
    parameters.grid = fv_grid_transverse_1d_v6(parameters);
end

%set up initial condition and perturbation
parameters.h_init = 1.5;                                                   % thickness at the ice divide                     
parameters.h = parameters.h_init;
parameters.iter.x_current = parameters.timestep.x_init;
parameters.flag_Tdep = 1;                                                  % turn on (1) and off(0) temperature dependent sliding
parameters.amplitude = 0;                                                  % amplitude of bed temperature perturbation; set to zero to compute unperturbed steady states; relevant to parameters.flag.pert = 'rand'
parameters.amplitude_gamma = 0;                                            % amplitude of bed friction perturbation; set to zero to compute unperturbed steady states; relevant to parameters.flag.pert = 'rand_gamma'
parameters.n_wavelengths = 2;                                              % number of wavelengths within assigneddomain size; relevant to parameters.flag.pert = 'mono'

%% DIVIDE INTEGRATION -> compute first timestep downstream of ice divide, assuming symmetry at the divide itself
T_nodes = parameters.grid.N.n_nodes.tot;                                
psi_nodes = parameters.grid.psi.n_nodes.tot; 

%construct initial guess
h_in = parameters.h;
Q_in = parameters.a*parameters.grid.N.extra.bd_y*parameters.timestep.dx/(2);
u_in = sparse(T_nodes,1);
T_in = parameters.T_bed - parameters.nu * parameters.grid.N.coor_nodes.z; 
phi_in = 1-parameters.a/2*(parameters.grid.N.coor_nodes.z).^2;

v_in = [phi_in; u_in; T_in; h_in; Q_in];

%set parameters for Newton solver
srch.tolF= 2e-04;
srch.verbose= 1;
srch.itmax = 50;
srch.toldelta = 1e-06;

[vout,error_flag,faux] = Newton_v2(@network_timestep_v5_divide,@network_timestep_v5_divide_jacobian,v_in,parameters,srch);

phi_index = 1:T_nodes;
u_index = T_nodes+1:2*T_nodes;
T_index = 2*T_nodes+1: 3*T_nodes;
h_index = 3*T_nodes+1;
Q_index = 3*T_nodes+2;

phi = vout(phi_index);
u = vout(u_index);
T = vout(T_index);
h = vout(h_index);
Q = vout(Q_index);

%% set up timestepping with vout as initial condition

%set variables at previous time step
parameters.v_in_prev.h_prev = h;
parameters.v_in_prev.u_prev = u;
parameters.v_in_prev.Q_prev = Q;
parameters.v_in_prev.qx_prev = sparse(length(parameters.grid.N.bdy_nodes.bottom),1);
parameters.v_in_prev.T_prev = T; 
parameters.v_in_prev.h_pprev = parameters.h_init;
parameters.v_in_prev.h_av_prev = (parameters.v_in_prev.h_pprev + parameters.v_in_prev.h_prev)/2;
parameters.v_in_prev.du_dz_centre_full_prev = faux.du_dz_centre_full;
parameters.v_in_prev.du_dy_centre_prev = sparse(T_nodes,1);
parameters.v_in_prev.u_vert_prev = faux.u_vert;
parameters.v_in_prev.I_prev = sparse(length(parameters.grid.N.bdy_nodes.bottom),1);

%set up initial condition/guess
h_in = h;
Q_in = Q;
u_in = u;
p_in = rand(T_nodes,1);
T_in = T; 
Tb_in = T(end) - parameters.nu* parameters.grid.Tb.coor_nodes.z;
phi_in = phi;
omega_in = rand(psi_nodes,1);
psi_in = rand(psi_nodes,1);
psi_ghost = 0;
Pi_in = zeros(length(parameters.grid.N.bdy_nodes.bottom),1);

v_in = [psi_in; omega_in; phi_in; u_in; p_in; T_in; Tb_in; psi_ghost; h_in; Q_in; Pi_in];

dir_out = '/Users/mantelli/Documents/ACADEMIC/research/postdoc_vancouver/code_local/code_thermal_finger/development/scratch';    %set directory where output is saved
timestepping_cluster(parameters, v_in, dir_out)





