function [fout, faux] = network_timestep_v5_mech_Tfixed_v2(v_in,parameters)
%finite-volume, steady solver for 3D thermal finger problem with periodic boundary conditions
%laterally. Same as network_timestep_v5, but for the fact that temperature
%is considered as known

%Jacobian is in network_timestep_v5_mech_Tfixed_v2
%Elisa Mantelli, Apr 2020

%parameters 
gamma = parameters.gamma;
a = parameters.a;
alpha = parameters.alpha;

%psi grid
psi_nodes = parameters.grid.psi.n_nodes.tot;                                   %number of nodes
psi_nodes_hor = parameters.grid.psi.n_nodes.hor;                               %number of nodes

psi_bdy_nodes_top = parameters.grid.psi.bdy_nodes.top;                         %list of nodes adjacent to the top of the box;
psi_bdy_nodes_bed = parameters.grid.psi.bdy_nodes.bed;

psi_Delta_z_cell = parameters.grid.psi.Delta_z_cell;                           %length of cells, ver (list)
psi_Delta_y_cell = parameters.grid.psi.Delta_y_cell;                           %length of cells, hor (list)
psi_Delta_y_partial_cell = parameters.grid.psi.delta_y_vflux;                  %length of hor cell esges crossed by vertical network edges, hor (list)

index_up_node_top_psi = 1:length(psi_bdy_nodes_top);
index_down_node_top_psi = circshift(1:length(psi_bdy_nodes_top), -1);

index_up_node_bed_psi = 1:length(psi_bdy_nodes_bed);
index_down_node_bed_psi = circshift(1:length(psi_bdy_nodes_bed), -1);

psi_connect_ver =  parameters.grid.psi.connect_ver;
psi_connect_hor = parameters.grid.psi.connect_hor;

%T grid
T_nodes = parameters.grid.N.n_nodes.tot;                                   %number of nodes
T_up_node_ver = parameters.grid.N.up_node.vert;                            %list (n_edges-by-1 vector) of 'upstream' node for each vertical edge

T_bdy_nodes_top = parameters.grid.N.bdy_nodes.top;                         %list of nodes adjacent to the top of the box;
T_bdy_nodes_bed = parameters.grid.N.bdy_nodes.bottom;

T_Delta_z_cell_volume = parameters.grid.N.Delta_z_cell_volume;             %length of cells, ver (list)
T_Delta_y_cell = parameters.grid.N.Delta_y_cell;                           %length of cells, hor (list)
T_Delta_y_partial_cell = parameters.grid.N.delta_y_vflux;                  %length of hor cell esges crossed by vertical network edges, hor (list)
T_coor_hor_cell_edges_z = parameters.grid.N.coor_hor_cell_edges.z;         %coordinate of horizontal cell edges (or, identically, vertical network edges)

index_up_node_top_T = circshift(1:length(T_bdy_nodes_top),1);
index_down_node_top_T = 1:length(T_bdy_nodes_top);

T_connect_ver =  parameters.grid.N.connect_ver;
T_connect_hor = parameters.grid.N.connect_hor;

%Tb grid
Tb_nodes = parameters.grid.Tb.n_nodes.tot;                                   %number of nodes

%unpack input variable v_in 
%v_in_pert_mech = [psi_in; omega_in; phi_in; u_in; p_in; psi_ghost; h_in; Q_in];

psi_index = 1:psi_nodes;
omega_index = psi_nodes+1: 2*psi_nodes;
phi_index = 2*psi_nodes+1:2*psi_nodes+T_nodes;
u_index = 2*psi_nodes+T_nodes+1:2*psi_nodes+2*T_nodes;
p_index = 2*psi_nodes+2*T_nodes+1: 2*psi_nodes+3*T_nodes;
psighost_index = 2*psi_nodes+3*T_nodes+1;
Q_index = 2*psi_nodes+3*T_nodes+2;

psi = v_in(psi_index);
omega = v_in(omega_index);
phi = v_in(phi_index);
u = v_in(u_index);
p = v_in(p_index);
h = parameters.h;
psi_ghost = v_in(psighost_index);

%unpack variable at previous timestep
h_prev = parameters.v_in_prev.h_prev;
u_prev = parameters.v_in_prev.u_prev;
h_pprev = parameters.v_in_prev.h_pprev;
h_av_prev = parameters.v_in_prev.h_av_prev;
u_vert_prev = parameters.v_in_prev.u_vert_prev;
du_dz_centre_full = parameters.v_in_prev.du_dz_centre_full;

%b = parameters.b; 

%initialize output
fout = zeros(length(v_in),1);
%%  PRELIMINARIES 

%DISCRETE DEP. VARIABLES
v_in_disc = [v_in(1:p_index(end)); sparse(T_nodes+Tb_nodes,1); psi_ghost; 0; 0; sparse(length(T_bdy_nodes_bed),1)];
[discvar, Ddiscvar] = discretisation_hanging_full_v4(v_in_disc, parameters);

dpsi_dy = discvar.dpsi_dy;
dpsi_dz = discvar.dpsi_dz;

domega_dy = discvar.domega_dy;
domega_dz = discvar.domega_dz;

u_vert = discvar.u_vert;
du_dy = discvar.du_dy;
du_dz = discvar.du_dz;
du_dy_centre = discvar.du_dy_centre;

dphi_dy = discvar.dphi_dy;
dphi_dz = discvar.dphi_dz;

dp_dy = discvar.dp_dy;
dp_dz = discvar.dp_dz;

%construct regularized bedwater content
T_bed = parameters.T_bed;
gamma_pert = parameters.gamma_pert;

Tbed_upnode = circshift(1:length(T_bdy_nodes_bed),1,2);
Tbed_downnode = 1:length(T_bdy_nodes_bed);
gamma_pert_psigrid = (gamma_pert(Tbed_upnode)+gamma_pert(Tbed_downnode))./2;
T_bed_psigrid = (T_bed(Tbed_upnode)+T_bed(Tbed_downnode))./2;

if parameters.flag_Tdep == 1
    [f_slide_Tbed, ~] = regularizedfriction_temperature(T_bed, parameters);
    [f_slide_Tbedpsigrid, ~] = regularizedfriction_temperature(T_bed_psigrid, parameters);
else
    f_slide_Tbed = ones(length(T_bdy_nodes_bed),1);
    f_slide_Tbedpsigrid = ones(length(psi_bdy_nodes_bed),1);
end


%TIMESTEPPING
% dx = parameters.timestep.dx;
% dx_prev =  parameters.timestep.dx_prev;
% dx_int =  parameters.timestep.dx_int;

x_current = parameters.timestep.x;
x_prev = x_current - (1/2*parameters.timestep.dx + 1/2*parameters.timestep.dx_prev);
x_pprev = x_prev - (1/2*parameters.timestep.dx_prev + 1/2*parameters.timestep.dx_pprev);

b = bed_finger(x_current, parameters); 
b_prev = bed_finger(x_prev, parameters); 
dx_prev = x_prev-x_pprev;  % x_i - x_{i-1}
dx_int = parameters.timestep.dx_prev; % x_{i+1/2} - x_{i-1/2}
dx = x_current -x_prev; % x_{i+1} - x_{i}

h_av = (h*parameters.timestep.dx_prev+h_prev*parameters.timestep.dx)./(parameters.timestep.dx_prev + parameters.timestep.dx); %eq. 16
%% STREAM FUNCTION, eq. 23 with bcs 28, 30
net_psi_horflux = psi_connect_hor*dpsi_dy;
net_psi_verflux = psi_connect_ver*(h_prev^(-1).*dpsi_dz.*psi_Delta_y_partial_cell);
% dirichlet conditions are applied at top and bottom. Implement second
% order accurate approximation

%top
psi_top = 0;
flux_psi_top = (8*psi_top -9*psi(psi_bdy_nodes_top) + psi(psi_bdy_nodes_top+psi_nodes_hor(1))).*psi_Delta_y_partial_cell(psi_bdy_nodes_top)./(3*psi_Delta_z_cell(psi_bdy_nodes_top).*h_prev);
net_psi_verflux(psi_bdy_nodes_top) = net_psi_verflux(psi_bdy_nodes_top) + flux_psi_top;

%bottom 
psi_bottom = psi_ghost;
flux_psi_bottom = -(8*psi_bottom - 9*psi(psi_bdy_nodes_bed) + psi(psi_bdy_nodes_bed-psi_nodes_hor(end)))./(3*psi_Delta_z_cell(psi_bdy_nodes_bed).*h_prev);
%flux_psi_bottom = h_prev.^(-1).*(psi(psi_bdy_nodes_bed)./psi_Delta_z_cell(psi_bdy_nodes_bed) - psi_bottom./psi_Delta_z_cell(psi_bdy_nodes_bed));
net_psi_verflux(psi_bdy_nodes_bed) = net_psi_verflux(psi_bdy_nodes_bed) - flux_psi_bottom.*psi_Delta_y_cell(psi_bdy_nodes_bed);

%conservation law
div_psifluxes =   net_psi_verflux./(h_prev.*psi_Delta_z_cell.*psi_Delta_y_cell)+ net_psi_horflux./psi_Delta_y_cell ;
fout(psi_index) = div_psifluxes - omega;

%% SLIDING LAW
%at T cell centres
u_bed = u(T_bdy_nodes_bed);

%at psi cell centres
index_up_node_bed = circshift(1:length(T_bdy_nodes_bed),1);
index_down_node_bed = 1:length(T_bdy_nodes_bed);
dphi_dy_bed = -(phi(T_bdy_nodes_bed(index_up_node_bed))-phi(T_bdy_nodes_bed(index_down_node_bed)))./psi_Delta_y_cell(psi_bdy_nodes_bed);
omega_bed = (gamma.*(1+gamma_pert_psigrid).*f_slide_Tbedpsigrid).*(dphi_dy_bed + flux_psi_bottom); 

%% additional condition for psi_ghost, int_width omega_bed = 0
fout(psighost_index) = sum(omega_bed.*psi_Delta_y_cell(psi_bdy_nodes_bed));  

%% VORTICITY, eq. 24 with bcs 29, 31
net_omega_horflux = psi_connect_hor * domega_dy;
net_omega_verflux = psi_connect_ver*(h_prev.^(-1).*domega_dz.*psi_Delta_y_partial_cell); 

%enforce Dirichlet conditions on top and bottom boundaries (psi_x=0 on
%inflow and outflow bdies)

%ice surface: omega = - (h-h_prev)/dx*1/2*(du/dy_{prev} +du/dy_{current})
du_dy_top = (u(T_bdy_nodes_top(index_up_node_top_T))-u(T_bdy_nodes_top(index_down_node_top_T)))./psi_Delta_y_cell(psi_bdy_nodes_top);
du_dy_top_prev = (u_prev(T_bdy_nodes_top(index_up_node_top_T))-u_prev(T_bdy_nodes_top(index_down_node_top_T)))./psi_Delta_y_cell(psi_bdy_nodes_top);

omega_top = -(h-h_prev)/(2*dx).*(du_dy_top_prev + du_dy_top);   
domegadztop = (-9*omega(psi_bdy_nodes_top) + omega(psi_bdy_nodes_top + length(psi_bdy_nodes_top)) +8*omega_top)./(3*psi_Delta_z_cell(psi_bdy_nodes_top));
net_omega_verflux(psi_bdy_nodes_top) = net_omega_verflux(psi_bdy_nodes_top) + h_prev.^(-1).*domegadztop.*psi_Delta_y_cell(psi_bdy_nodes_top);

%bed: omega = omega_bed -((8 a - 9 f1 + f2)/(3 z))
domegadz_bed =  (9*omega(psi_bdy_nodes_bed) - omega(psi_bdy_nodes_bed - length(psi_bdy_nodes_bed))-8*omega_bed)./(3*psi_Delta_z_cell(psi_bdy_nodes_bed));
net_omega_verflux(psi_bdy_nodes_bed) = net_omega_verflux(psi_bdy_nodes_bed) - h_prev.^(-1).*domegadz_bed.*psi_Delta_y_cell(psi_bdy_nodes_bed);

%conservation law
div_omegafluxes =   net_omega_verflux./(h_prev.*psi_Delta_z_cell.*psi_Delta_y_cell)+ net_omega_horflux./psi_Delta_y_cell ;
fout(omega_index) = div_omegafluxes;

%% PHI, eq. 22 with bcs 26, 27

%diffusive fluxes
net_phi_horflux = T_connect_hor * dphi_dy;
net_phi_diff_verflux = T_connect_ver * (h_prev.^(-1).*dphi_dz.*T_Delta_y_partial_cell);

%advective flux
flux_phi_adv = -1/2*T_coor_hor_cell_edges_z.*(u_vert.*(h - h_prev)/dx + (h_prev - h_pprev)/dx_prev*u_vert_prev);
net_phi_adv_verflux = T_connect_ver * (flux_phi_adv.*T_Delta_y_partial_cell);

net_phi_verflux =  net_phi_diff_verflux + net_phi_adv_verflux;

%surface condition, eq. 26
flux_phi_top = - a;
net_phi_verflux(T_bdy_nodes_top) = net_phi_verflux(T_bdy_nodes_top) + flux_phi_top.*T_Delta_y_cell(T_bdy_nodes_top);

div_phi_fluxes = net_phi_verflux./(h_prev*T_Delta_z_cell_volume.*T_Delta_y_cell)+ net_phi_horflux./T_Delta_y_cell;

%source term
S_phi = (h_av*u - h_av_prev*u_prev)./dx_int;

%conservation law 
fout(phi_index) = div_phi_fluxes + h_prev.^(-1)*S_phi;
fout(phi_index(1)) = phi(1)-1;

%% ALONG FLOW VELOCITY, U, eq. 19 with 20-21 (11/10)
%diffusive fluxes
net_u_horflux = T_connect_hor * du_dy;
net_u_verflux = T_connect_ver * (h_av.^(-1).*du_dz.*T_Delta_y_partial_cell);

%eq. 21: prescribed flux at cell centre at the bed
flux_u_bottom  = (gamma.*(1+gamma_pert).*f_slide_Tbed.*u_bed).*T_Delta_y_cell(T_bdy_nodes_bed); 
net_u_verflux(T_bdy_nodes_bed) = net_u_verflux(T_bdy_nodes_bed) - flux_u_bottom;

div_u_fluxes =   net_u_verflux./(h_av.*T_Delta_z_cell_volume.*T_Delta_y_cell)+ net_u_horflux./T_Delta_y_cell; %nb: h_int depends on h, so must be differentiated!

%source term
S_u = (h + b - h_prev - b_prev)/(dx).* ones(T_nodes, 1);

%conservation law
fout(u_index) = div_u_fluxes - S_u;

%% P', eq. 25 with 32 (11/10)
% NOTE: prescribe dirichlet condition on *one* surface node so as to fix p

%diffusive fluxes
net_p_horflux = T_connect_hor * dp_dy;
net_p_verflux = T_connect_ver * (h_prev^(-1).*dp_dz.*T_Delta_y_partial_cell);

%eq.32: prescribed flux at top nodes
domegady_top = (omega_top(index_down_node_top_psi) - omega_top(index_up_node_top_psi))./T_Delta_y_cell(T_bdy_nodes_top);
flux_p_top = -domegady_top.*T_Delta_y_cell(T_bdy_nodes_top); 
net_p_verflux(T_bdy_nodes_top) = net_p_verflux(T_bdy_nodes_top) + flux_p_top;

%eq. 32: prescribed flux at cell centre, bed
domegady_bottom = (omega_bed(index_down_node_bed_psi) - omega_bed(index_up_node_bed_psi))./T_Delta_y_cell(T_bdy_nodes_bed);
flux_p_bottom  = -domegady_bottom.*T_Delta_y_cell(T_bdy_nodes_bed); 
net_p_verflux(T_bdy_nodes_bed) = net_p_verflux(T_bdy_nodes_bed) - flux_p_bottom;

div_p_fluxes =   net_p_verflux./(h_prev.*T_Delta_z_cell_volume.*T_Delta_y_cell)+ net_p_horflux./T_Delta_y_cell;

%conservation law 
fout(p_index) = div_p_fluxes;
fout(p_index(1)) = p(1);

%% AUXILIARY VARIABLES
faux.du_dz_centre_full = du_dz_centre_full;
faux.du_dy_centre = du_dy_centre;
faux.u_vert = u_vert;
faux.h_prev = h_prev;
faux.u_b = u(T_bdy_nodes_bed);
faux.friction = alpha*(gamma.*(1+gamma_pert).*f_slide_Tbed).*((u(T_bdy_nodes_bed)+u_prev(T_bdy_nodes_bed))/2).^(2);
faux.surfaceslope = (h + b - h_prev - b_prev)/(dx);
faux.S = du_dz_centre_full;
faux.h_av = h_av;
faux.Q = h*sum(u.*T_Delta_z_cell_volume.*T_Delta_y_cell);

%% compute transverse velocities
%transverse velocity, eq. 34
dpsi_dz_full = [flux_psi_top; dpsi_dz; flux_psi_bottom];
faux.v = dpsi_dz; 
faux.dpsi_dy = dpsi_dy; 

[dpsi_dy_reindexed,~] = verticalvelocity(parameters,dpsi_dy,Ddiscvar.ddpsidy_dpsi);
faux.w = -dpsi_dy_reindexed ;
faux.f = (gamma.*(1+gamma_pert).*f_slide_Tbed.*u_bed);








