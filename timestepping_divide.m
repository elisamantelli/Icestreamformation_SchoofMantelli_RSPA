function timestepping_divide(parameters)

%% DIVIDE INTEGRATION
T_nodes = parameters.grid.N.n_nodes.tot;                                
psi_nodes = parameters.grid.psi.n_nodes.tot; 

%construct initial guess
h_in = parameters.h;
Q_in = parameters.a*parameters.grid.N.extra.bd_y*parameters.timestep.dx/(2);
u_in = sparse(T_nodes,1);
T_in = parameters.T_s - parameters.nu * parameters.grid.N.coor_nodes.z; 
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
p_in = sparse(T_nodes,1);
T_in = T; 
Tb_in = T(end) - parameters.nu* parameters.grid.Tb.coor_nodes.z;
phi_in = phi;
omega_in = sparse(psi_nodes,1);
psi_in = sparse(psi_nodes,1);
psi_ghost = 0;
Pi_in = zeros(length(parameters.grid.N.bdy_nodes.bottom),1);

v_in = [psi_in; omega_in; phi_in; u_in; p_in; T_in; Tb_in; psi_ghost; h_in; Q_in; Pi_in];

timestepping_cluster(parameters, v_in)
end