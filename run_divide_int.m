N = 8;
maxNumCompThreads(N)

%set physical parameters
parameters.T_s = -1; 
parameters.nu = 0.5;
parameters.alpha = 1;
parameters.gamma = 0.1;
parameters.mu = 0.1;
parameters.Pe = 1;
parameters.a = 1;
parameters.r = 917/1000;
parameters.beta = 0.1;
parameters.T_bed = -0.2;

parameters.flag.heat_full = 1;
parameters.flag.plug = 0;
parameters.flag.pert = 'rand_gamma'; % set style of perturbation. 'rand_gamma' perturbs friction coefficient
parameters.flag.hydrology = 'weertman';% 'budd', 'coulomb'
parameters.flag.fluxy_sigma = 'Pionly';
parameters.flag.hydropermnew = 1;

%regularization parameters 
%friction coefficient
parameters.reg.epsilon = 0.005; %cutoff for regularised heaviside function. Should be smaller than epsilon_f
parameters.reg.epsilon_f = 0.03; %delta

 %permeability
parameters.reg.c = 3; %keep as strictly larger than 1. Exponent in the permeability power law
parameters.reg.Pi_0 = 0.1; %this sets the derivative of the regularization for Pi=0. Eventually should be small. 
parameters.reg.Pi_s = 0.1;

%bed elevation
parameters. bed. b0 = 0; %elevation at the divide, x= 0
parameters. bed. b1 = 0.05; %bed slope, must be strictly positive

parameters.timestep.x_init = 0.00036;%parameters.reg.epsilon_f/10;
parameters.timestep.dx = 0.00036;%parameters.reg.epsilon_f/10; %(i+1)
parameters.timestep.dx_prev =  0.00036;%parameters.reg.epsilon_f/10; %i
parameters.timestep.dx_pprev =  0.00036;%parameters.reg.epsilon_f/10; %i-1
parameters.timestep.x = 0.00036;%parameters.timestep.x_init; % update at each timestep
parameters.tplot = parameters.timestep.dx;
parameters.stepmax = 3125000;

parameters.da = 0;    %parameters.timestep.dx;


%grid generation
parameters.flag1d = 1;
if parameters.flag1d == 1
    parameters.n_nodes_transverse = 1;
else
    parameters.n_nodes_transverse = 8;   %#nodes in the transverse direction in the uppermost layer. Remaining grid parameters are set in the grid routine. Must be even!
    if rem(parameters.n_nodes_transverse,2)>0
        warning('transverse cell number must be even')
    end
end
parameters.ratio_hor = 7;  %dy = dz*ratio
parameters.ratio_vert = 4; 

if parameters.flag1d == 0
    parameters.grid = fv_grid_transverse_v6(parameters);
elseif parameters.flag1d == 1
    parameters.grid = fv_grid_transverse_1d_v6(parameters);
end

%initial condition with perturbation
parameters.h_init = 1.5; 
parameters.h = parameters.h_init;
parameters.iter.x_current = parameters.timestep.x_init;
parameters.flag_Tdep = 1; %turn on (1) and off(0) T_dep sliding
parameters.amplitude = 0;    %amplitude of bed temperature perturbation; set to zero to compute steady states
parameters.amplitude_gamma = 0;    %amplitude of bed temperature perturbation; set to zero to compute steady states
parameters.n_wavelengths = 2; %number of wavelengths within assigneddomain size

%% DIVIDE INTEGRATION
T_nodes = parameters.grid.N.n_nodes.tot;                                
psi_nodes = parameters.grid.psi.n_nodes.tot; 

%construct initial guess
h_in = parameters.h;
Q_in = parameters.a*parameters.grid.N.extra.bd_y*parameters.timestep.dx/(2);
u_in = sparse(T_nodes,1);
T_in = parameters.T_bed - parameters.nu * parameters.grid.N.coor_nodes.z; 
phi_in = 1-parameters.a/2*(parameters.grid.N.coor_nodes.z).^2;

v_in = [phi_in; u_in; T_in; h_in; Q_in];

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

%% set up timestepping with this as initial condition

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

dir_out = '/Users/mantelli/Documents/ACADEMIC/research/postdoc_vancouver/code_local/code_thermal_finger/development/scratch';
timestepping_cluster(parameters, v_in, dir_out)





