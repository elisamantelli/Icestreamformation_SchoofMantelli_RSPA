function [fout, faux] = network_timestep_v5_divide(v_in,parameters)
%divide calculation with the same code as in network_timestep_v5
%nb: as set up, the code assumes the bed to be at most subtemperate at the divide and in its proximity, that is, Pi = 0 and no derivatives wrt Pi.

%parameters 
gamma = parameters.gamma;
a = parameters.a;
Pe = parameters.Pe;
alpha = parameters.alpha;
T_s = parameters.T_s;

%psi grid
psi_nodes = parameters.grid.psi.n_nodes.tot;                                   %number of nodes

%T grid
T_nodes = parameters.grid.N.n_nodes.tot;                                   %number of nodes
T_nodes_hor = parameters.grid.N.n_nodes.hor;                               %number of nodes, hor
T_bdy_nodes_top = parameters.grid.N.bdy_nodes.top;                         %list of nodes adjacent to the top of the box;
T_bdy_nodes_bed = parameters.grid.N.bdy_nodes.bottom;

T_Delta_z_cell = parameters.grid.N.Delta_z_cell;                           %length of cells, ver (list)
T_Delta_z_cell_volume = parameters.grid.N.Delta_z_cell_volume;             %length of cells, ver (list)
T_Delta_y_cell = parameters.grid.N.Delta_y_cell;                           %length of cells, hor (list)
T_Delta_y_partial_cell = parameters.grid.N.delta_y_vflux;                  %length of hor cell esges crossed by vertical network edges, hor (list)
T_coor_hor_cell_edges_z = parameters.grid.N.coor_hor_cell_edges.z;         %coordinate of horizontal cell edges (or, identically, vertical network edges)

T_connect_ver =  parameters.grid.N.connect_ver;

%Tb grid
Tb_nodes = parameters.grid.Tb.n_nodes.tot;                                   %number of nodes

%unpack input variable v_in 

phi_index = 1:T_nodes;
u_index = T_nodes+1:2*T_nodes;
T_index = 2*T_nodes+1: 3*T_nodes;
h_index = 3*T_nodes+1;
Q_index = 3*T_nodes+2;

phi = v_in(phi_index);
u = v_in(u_index);
T = v_in(T_index);
h = v_in(h_index);
Q = v_in(Q_index);

h_prev = parameters.h_init;  %this is h at the divide, fixed; we're solving for h_2

%initialize output
fout = zeros(length(v_in),1);
%%  PRELIMINARIES 

%DISCRETE DEP. VARIABLES
v_in_disc = [sparse(2*psi_nodes,1); phi; u; sparse(T_nodes,1); T; sparse(Tb_nodes+1,1); h; Q; sparse(length(T_bdy_nodes_bed),1)];
[discvar, ~] = discretisation_hanging_full_v4(v_in_disc, parameters);

u_vert = discvar.u_vert;
du_dz = discvar.du_dz;
du_dz_centre = discvar.du_dz_centre;

dphi_dz = discvar.dphi_dz;

dT_dz = discvar.dT_dz;
T_vert = discvar.T_vert;


%construct regularized bedwater content
T_bed = T(T_bdy_nodes_bed);
u_bed = u(T_bdy_nodes_bed);

if isfield(parameters,'gamma_pert')==1
    gamma_pert = parameters.gamma_pert;
elseif isfield(parameters,'gamma_pert')==0
    gamma_pert = zeros(length(T_bdy_nodes_bed),1);
end

if parameters.flag_Tdep == 1
    [f_slide_Tbed, ~] = regularizedfriction_temperature(T_bed, parameters);
elseif parameters.flag_Tdep == 0
    f_slide_Tbed = ones(length(T_bdy_nodes_bed),1);
end

[f_slide_hydro, ~] = regularizedfriction_hydrology(u_bed,zeros(length(u_bed),1), parameters);

%TIMESTEPPING
x_current = parameters.timestep.x; %i+1
x_prev = 0; %x_i
x_pprev = -x_current; %x_i-1

b = bed_finger(x_current, parameters); 
b_prev = bed_finger(x_prev, parameters); 
dx_prev = x_prev-x_pprev;  % x_i - x_{i-1}
dx_int = parameters.timestep.dx_prev; % x_{i+1/2} - x_{i-1/2}
dx = x_current -x_prev; % x_{i+1} - x_{i}

h_av = (h*parameters.timestep.dx_prev+h_prev*parameters.timestep.dx)./(parameters.timestep.dx_prev + parameters.timestep.dx); %eq. 16

%% SIMMETRY ACROSS THE DIVIDE

h_av_prev = h_av;
u_prev = -u;  
Q_prev = -Q;

T_prev = T; 
h_pprev = h; 
u_vert_prev = -u_vert; 
%% MASS CONSERVATION AND EQ. FOR h, eqs. 17-18
fout(Q_index) = (Q- Q_prev)./(dx_int)- (a*parameters.grid.N.extra.bd_y);

u_int = sum(u.*T_Delta_z_cell_volume.*T_Delta_y_cell);

fout(h_index) = Q - h_av.*u_int;

%% PHI, eq. 22 with bcs 26, 27
%diffusive flux
net_phi_diff_verflux = T_connect_ver * (h_prev.^(-1).*dphi_dz.*T_Delta_y_partial_cell);

%advective flux
flux_phi_adv = -1/2*T_coor_hor_cell_edges_z.*(u_vert.*(h - h_prev)/dx + (h_prev - h_pprev)/dx_prev*u_vert_prev);
net_phi_adv_verflux = T_connect_ver * (flux_phi_adv.*T_Delta_y_partial_cell);

net_phi_verflux =  net_phi_diff_verflux + net_phi_adv_verflux;

%surface condition, eq. 26
flux_phi_top = - a;
net_phi_verflux(T_bdy_nodes_top) = net_phi_verflux(T_bdy_nodes_top) + flux_phi_top.*T_Delta_y_cell(T_bdy_nodes_top);

div_phi_fluxes = net_phi_verflux./(h_prev*T_Delta_z_cell_volume.*T_Delta_y_cell);

%source term
S_phi = (h_av*u - h_av_prev*u_prev)./dx_int;

%conservation law 
fout(phi_index) = div_phi_fluxes + h_prev.^(-1)*S_phi;
fout(phi_index(T_bdy_nodes_top)) = phi(T_bdy_nodes_top)-1;

%% ALONG FLOW VELOCITY, U, eq. 19 with 20-21 (11/10)
%diffusive fluxes
net_u_verflux = T_connect_ver *( h_av.^(-1).*du_dz.*T_Delta_y_partial_cell);

%eq. 21: prescribed flux at cell centre at the bed
flux_u_bottom  = (gamma.*(1+gamma_pert).*f_slide_Tbed.*f_slide_hydro.*u_bed).*T_Delta_y_cell(T_bdy_nodes_bed); 
net_u_verflux(T_bdy_nodes_bed) = net_u_verflux(T_bdy_nodes_bed) - flux_u_bottom;

div_u_fluxes =   net_u_verflux./(h_av.*T_Delta_z_cell_volume.*T_Delta_y_cell); %nb: h_int depends on h, so must be differentiated!

%source term
S_u = (h + b - h_prev - b_prev)/(dx).* ones(T_nodes, 1);

%conservation law
fout(u_index) = div_u_fluxes - S_u;

%% CONSTITUTIVE RELATION for w_eff

w_eff =  h_prev.^(-1).*dphi_dz - T_coor_hor_cell_edges_z./2.*(u_vert *(h-h_prev)/dx + u_vert_prev .*(h_prev-h_pprev)/dx_prev); 

%% HEAT EQUATION, ICE, eq. 36-37 with bcs. 39 and 41 
%diffusive fluxes
net_T_diffverflux = T_connect_ver *( dT_dz.*T_Delta_y_partial_cell);

%advective fluxes 
T_adv_flux_ver = w_eff.*T_vert;

net_T_advfluxver = T_connect_ver *(T_adv_flux_ver.*T_Delta_y_partial_cell);

net_T_advfluxalong = (u.*h_prev.*T -u_prev.*h_pprev.*T_prev)./dx_int;
%source term
du_dz_centre_top = sparse(length(T_bdy_nodes_top),1);
du_dz_centre_bottom =  h*gamma.*(1+gamma_pert).*f_slide_Tbed.*f_slide_hydro.*u_bed;
du_dz_centre_full = [du_dz_centre_top; du_dz_centre; du_dz_centre_bottom];

S_T = (h_prev)^(-2).*((du_dz_centre_full).^2); 

%conservation law
fout(T_index) = ...
    Pe/h_prev.*( net_T_advfluxalong + net_T_advfluxver./(T_Delta_z_cell.*T_Delta_y_cell))- ...
    net_T_diffverflux./(h_prev.^2*T_Delta_z_cell.*T_Delta_y_cell) - alpha*S_T;

%surface and bed boundary condition
fout(T_index(T_bdy_nodes_top)) = T(T_bdy_nodes_top) - T_s;
%fout(2*psi_nodes+3*T_nodes+T_bdy_nodes_bed) = T(T_bdy_nodes_bed) - T_bed;

%% BASAL ENERGY BUDGET, eq. 41 
%this is effectively a constraint on basal temperature, and should be
%treated as an additional equation in the unknown T_bed

bedflux_above = - (T(T_bdy_nodes_bed-T_nodes_hor(end))- T(T_bdy_nodes_bed))./(h_prev*T_Delta_z_cell(T_bdy_nodes_bed));   
bedflux_below = parameters.nu;   
m = -(bedflux_above - bedflux_below) + alpha*(gamma.*(1+gamma_pert).*f_slide_Tbed.*f_slide_hydro).*u(T_bdy_nodes_bed).^(2);

enth_eq = - m;

%constrain Pi where bed was temperate at previous timestep and Tbed where
%the bed was frozen
fout(T_index(T_bdy_nodes_bed)) = enth_eq;

%% AUXILIARY VARIABLES
faux.du_dz_centre_full = du_dz_centre_full;
faux.u_vert = u_vert;
faux.u_b = u(T_bdy_nodes_bed);
faux.bedflux = bedflux_above;
faux.bedflux_bed = bedflux_below;
faux.friction = alpha*(gamma.*f_slide_Tbed).*((u(T_bdy_nodes_bed))).^(2);
faux.surfaceslope = (h + b - h_prev - b_prev)/(dx);
faux.m = m;
faux.S = du_dz_centre_full;
faux.w = w_eff;
faux.h_av = h_av;













