function [fout, faux] = network_timestep_v5(v_in,parameters)
%finite-volume, steady solver for 3D thermal finger problem with periodic boundary conditions
%laterally 

%Considers half cells on top
%and bottom layer for T in the ice and mechanical variables living on T
%grid. Tthe computation of the flux quadrature is formulated in such
%a way to satisphy the solvability condition for phi in discrete form. Second-order accurate implementation of boundary
%condition for psi at the bed

%Jacobian is in network_timestep_v5_jacobian
%Elisa Mantelli, Jan 2021

%parameters 
gamma = parameters.gamma;
a = parameters.a;
Pe = parameters.Pe;
alpha = parameters.alpha;
nu = parameters.nu;
T_s = parameters.T_s;
r = parameters.r;
beta = parameters.beta;

%psi grid
psi_nodes = parameters.grid.psi.n_nodes.tot;                                   %number of nodes
psi_nodes_hor = parameters.grid.psi.n_nodes.hor;                               %number of nodes

psi_bdy_nodes_top = parameters.grid.psi.bdy_nodes.top;                         %list of nodes adjacent to the top of the box;
psi_bdy_nodes_bed = parameters.grid.psi.bdy_nodes.bed;

psi_Delta_z_cell = parameters.grid.psi.Delta_z_cell;                           %length of cells, ver (list)
psi_Delta_y_cell = parameters.grid.psi.Delta_y_cell;                           %length of cells, hor (list)
psi_Delta_y_partial_cell = parameters.grid.psi.delta_y_vflux;                  %length of hor cell esges crossed by vertical network edges, hor (list)

index_up_node_top_psi = 1:length(psi_bdy_nodes_top);
index_down_node_top_psi = circshift(1:length(psi_bdy_nodes_top), -1,2);

index_up_node_bed_psi = 1:length(psi_bdy_nodes_bed);
index_down_node_bed_psi = circshift(1:length(psi_bdy_nodes_bed),-1,2);

psi_connect_ver =  parameters.grid.psi.connect_ver;
psi_connect_hor = parameters.grid.psi.connect_hor;

%T grid
T_nodes = parameters.grid.N.n_nodes.tot;                                   %number of nodes
T_nodes_hor = parameters.grid.N.n_nodes.hor;                                   %number of nodes

T_bdy_nodes_top = parameters.grid.N.bdy_nodes.top;                         %list of nodes adjacent to the top of the box;
T_bdy_nodes_bed = parameters.grid.N.bdy_nodes.bottom;

T_Delta_z_cell = parameters.grid.N.Delta_z_cell;                           %length of cells, ver (list)
T_Delta_z_cell_volume = parameters.grid.N.Delta_z_cell_volume;             %length of cells, ver (list)
T_Delta_y_cell = parameters.grid.N.Delta_y_cell;                           %length of cells, hor (list)
T_Delta_y_partial_cell = parameters.grid.N.delta_y_vflux;                  %length of hor cell esges crossed by vertical network edges, hor (list)
T_coor_hor_cell_edges_z = parameters.grid.N.coor_hor_cell_edges.z;         %coordinate of horizontal cell edges (or, identically, vertical network edges)

index_up_node_top_T = circshift(1:length(T_bdy_nodes_top),1,2);
index_down_node_top_T = 1:length(T_bdy_nodes_top);

T_connect_ver =  parameters.grid.N.connect_ver;
T_connect_hor = parameters.grid.N.connect_hor;

%Tb grid
Tb_nodes = parameters.grid.Tb.n_nodes.tot;                                   %number of nodes

Tb_bdy_nodes_top = parameters.grid.Tb.bdy_nodes.top;                         %list of nodes adjacent to the top of the box;
Tb_bdy_nodes_bed = parameters.grid.Tb.bdy_nodes.bottom;

Tb_Delta_z_cell = parameters.grid.Tb.Delta_z_cell;                           %length of cells, ver (list)
Tb_Delta_y_cell = parameters.grid.Tb.Delta_y_cell;                           %length of cells, hor (list)
Tb_Delta_y_partial_cell = parameters.grid.Tb.delta_y_vflux;                  %length of hor cell esges crossed by vertical network edges, hor (list)

Tb_connect_ver =  parameters.grid.Tb.connect_ver;
Tb_connect_hor = parameters.grid.Tb.connect_hor;

%unpack input variable v_in 

psi_index = 1:psi_nodes;
omega_index = psi_nodes+1: 2*psi_nodes;
phi_index = 2*psi_nodes+1:2*psi_nodes+T_nodes;
u_index = 2*psi_nodes+T_nodes+1:2*psi_nodes+2*T_nodes;
p_index = 2*psi_nodes+2*T_nodes+1: 2*psi_nodes+3*T_nodes;
T_index = 2*psi_nodes+3*T_nodes+1: 2*psi_nodes+4*T_nodes;
Tb_index = 2*psi_nodes+4*T_nodes+1: 2*psi_nodes+4*T_nodes+Tb_nodes;
h_index = 2*psi_nodes+4*T_nodes+Tb_nodes+2;
Q_index = 2*psi_nodes+4*T_nodes+Tb_nodes+3;
psighost_index = 2*psi_nodes+4*T_nodes+Tb_nodes+1;
Pi_index =2*psi_nodes+4*T_nodes+Tb_nodes+4: 2*psi_nodes+4*T_nodes+Tb_nodes+3+length(T_bdy_nodes_bed);

psi = v_in(psi_index);
omega = v_in(omega_index);
phi = v_in(phi_index);
u = v_in(u_index);
p = v_in(p_index);
T = v_in(T_index);
Tb = v_in(Tb_index);
h = v_in(h_index);
Q = v_in(Q_index);
psi_ghost = v_in(psighost_index);
Pi = v_in(Pi_index);

%unpack variable at previous timestep
h_prev = parameters.v_in_prev.h_prev;
u_prev = parameters.v_in_prev.u_prev;
Q_prev = parameters.v_in_prev.Q_prev;
%b_prev = parameters.v_in_prev.b_prev;
qx_prev = parameters.v_in_prev.qx_prev;

if parameters.flag.plug == 1
    u_prev_T = parameters.v_in_prev.u_prev_T;
    h_prev_T = parameters.v_in_prev.h_T;
end

T_prev = parameters.v_in_prev.T_prev; 
h_pprev = parameters.v_in_prev.h_pprev;
h_av_prev = parameters.v_in_prev.h_av_prev;
du_dz_centre_full_prev = parameters.v_in_prev.du_dz_centre_full_prev;
du_dy_centre_prev = parameters.v_in_prev.du_dy_centre_prev;
u_vert_prev = parameters.v_in_prev.u_vert_prev;
I_prev = parameters.v_in_prev.I_prev;
if isfield(parameters.flag, 'sigma') == 1
    f_prev = parameters.v_in_prev.f_prev;
end

%b = parameters.b; 

%initialize output
fout = zeros(length(v_in),1);
%%  PRELIMINARIES 

%DISCRETE DEP. VARIABLES
[discvar, Ddiscvar] = discretisation_hanging_full_v4(v_in, parameters);

dpsi_dy = discvar.dpsi_dy;
dpsi_dz = discvar.dpsi_dz;

domega_dy = discvar.domega_dy;
domega_dz = discvar.domega_dz;

u_vert = discvar.u_vert;
du_dy = discvar.du_dy;
du_dz = discvar.du_dz;
du_dz_centre = discvar.du_dz_centre;
du_dy_centre = discvar.du_dy_centre;

dphi_dy = discvar.dphi_dy;
dphi_dz = discvar.dphi_dz;

dp_dy = discvar.dp_dy;
dp_dz = discvar.dp_dz;

dT_dy = discvar.dT_dy;
dT_dz = discvar.dT_dz;
T_vert = discvar.T_vert;
T_hor = discvar.T_hor;

dTb_dy = discvar.dTb_dy;
dTb_dz = discvar.dTb_dz;

%construct regularized friction coefficient
T_bed = T(T_bdy_nodes_bed);
u_bed = u(T_bdy_nodes_bed);

%T_bed, u_bed, Pi at psi cell centres
Tbed_upnode = circshift(1:length(T_bdy_nodes_bed),1,2);
Tbed_downnode = 1:length(T_bdy_nodes_bed);
T_bed_psigrid = (T_bed(Tbed_upnode)+T_bed(Tbed_downnode))./2;
u_bed_psigrid = (u_bed(Tbed_upnode)+u_bed(Tbed_downnode))./2;
Pi_psigrid = (Pi(Tbed_upnode)+Pi(Tbed_downnode))./2;

if isfield(parameters,'gamma_pert')==1
    gamma_pert = parameters.gamma_pert;
    gamma_pert_psigrid = (gamma_pert(Tbed_upnode)+gamma_pert(Tbed_downnode))./2;
elseif isfield(parameters,'gamma_pert')==0
    gamma_pert = zeros(length(T_bdy_nodes_bed),1);
    gamma_pert_psigrid = zeros(length(psi_bdy_nodes_bed),1);
end

if parameters.flag_Tdep == 1
    [f_slide_Tbed, ~] = regularizedfriction_temperature(T_bed, parameters);
    [f_slide_Tbedpsigrid, ~] = regularizedfriction_temperature(T_bed_psigrid, parameters);
elseif parameters.flag_Tdep == 0
    f_slide_Tbed = ones(length(T_bdy_nodes_bed),1);
    f_slide_Tbedpsigrid = ones(length(psi_bdy_nodes_bed),1);
end

[f_slide_hydro, ~] = regularizedfriction_hydrology(u_bed,Pi, parameters);
[f_slide_hydro_psigrid, ~] = regularizedfriction_hydrology(u_bed_psigrid,Pi_psigrid, parameters);

%TIMESTEPPING
x_current = parameters.timestep.x;
x_prev = x_current - (1/2*parameters.timestep.dx + 1/2*parameters.timestep.dx_prev);
x_pprev = x_prev - (1/2*parameters.timestep.dx_prev + 1/2*parameters.timestep.dx_pprev);

b = bed_finger(x_current, parameters); 
b_prev = bed_finger(x_prev, parameters); 
dx_prev = x_prev-x_pprev;  % x_i - x_{i-1}
dx_int = parameters.timestep.dx_prev; % x_{i+1/2} - x_{i-1/2}
dx = x_current -x_prev; % x_{i+1} - x_{i}

h_av = (h*parameters.timestep.dx_prev+h_prev*parameters.timestep.dx)./(parameters.timestep.dx_prev + parameters.timestep.dx); %eq. 16

%% MASS CONSERVATION AND EQ. FOR h
fout(Q_index) = (Q -Q_prev)./(dx_int)- (a*parameters.grid.N.extra.bd_y);

u_int = sum(u.*T_Delta_z_cell_volume.*T_Delta_y_cell);

fout(h_index) = Q - h_av.*u_int;

%% STREAM FUNCTION

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
net_psi_verflux(psi_bdy_nodes_bed) = net_psi_verflux(psi_bdy_nodes_bed) - flux_psi_bottom.*psi_Delta_y_cell(psi_bdy_nodes_bed);

%conservation law
div_psifluxes =   net_psi_verflux./(h_prev.*psi_Delta_z_cell.*psi_Delta_y_cell)+ net_psi_horflux./psi_Delta_y_cell ;
fout(psi_index) = div_psifluxes - omega;


%% SLIDING LAW

%at psi cell centres
index_up_node_bed = circshift(1:length(T_bdy_nodes_bed),1,2);
index_down_node_bed = 1:length(T_bdy_nodes_bed);
dphi_dy_bed = (-phi(T_bdy_nodes_bed(index_up_node_bed))+phi(T_bdy_nodes_bed(index_down_node_bed)))./psi_Delta_y_cell(psi_bdy_nodes_bed);
omega_bed = (gamma.*(1+gamma_pert_psigrid).*f_slide_Tbedpsigrid.*f_slide_hydro_psigrid).*(dphi_dy_bed + flux_psi_bottom); 

%% additional condition for psi_ghost, int_width omega_bed = 0
fout(psighost_index) = sum(omega_bed.*psi_Delta_y_cell(psi_bdy_nodes_bed));  

%% VORTICITY
net_omega_horflux = psi_connect_hor * domega_dy;
net_omega_verflux = psi_connect_ver*(h_prev.^(-1).*domega_dz.*psi_Delta_y_partial_cell);

%enforce Dirichlet conditions on top and bottom boundaries (psi_x=0 on
%inflow and outflow bdies)

%ice surface: omega = - (h-h_prev)/dx*1/2*(du/dy_{prev} +du/dy_{current})
du_dy_top = (u(T_bdy_nodes_top(index_up_node_top_T))-u(T_bdy_nodes_top(index_down_node_top_T)))./psi_Delta_y_cell(psi_bdy_nodes_top);
du_dy_top_prev = (u_prev(T_bdy_nodes_top(index_up_node_top_T))-u_prev(T_bdy_nodes_top(index_down_node_top_T)))./psi_Delta_y_cell(psi_bdy_nodes_top);

omega_top = (h-h_prev)/(2*dx).*(du_dy_top_prev + du_dy_top);   
domegadztop = (-9*omega(psi_bdy_nodes_top) + omega(psi_bdy_nodes_top + length(psi_bdy_nodes_top)) +8*omega_top)./(3*psi_Delta_z_cell(psi_bdy_nodes_top));
net_omega_verflux(psi_bdy_nodes_top) = net_omega_verflux(psi_bdy_nodes_top) + h_prev.^(-1).*domegadztop.*psi_Delta_y_cell(psi_bdy_nodes_top);

%bed: omega = omega_bed -((8 a - 9 f1 + f2)/(3 z))
domegadz_bed =  (9*omega(psi_bdy_nodes_bed) - omega(psi_bdy_nodes_bed - length(psi_bdy_nodes_bed))-8*omega_bed)./(3*psi_Delta_z_cell(psi_bdy_nodes_bed));
net_omega_verflux(psi_bdy_nodes_bed) = net_omega_verflux(psi_bdy_nodes_bed) - h_prev.^(-1).*psi_Delta_y_cell(psi_bdy_nodes_bed).*domegadz_bed;

%conservation law
div_omegafluxes =   net_omega_verflux./(h_prev.*psi_Delta_z_cell.*psi_Delta_y_cell)+ net_omega_horflux./psi_Delta_y_cell ;
fout(omega_index) = div_omegafluxes;

%% PHI

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
%dirichlet condition on top 
fout(phi_index(1)) = phi(1)-1;


%% ALONG FLOW VELOCITY, U
%diffusive fluxes
net_u_horflux = T_connect_hor * du_dy;
net_u_verflux = T_connect_ver * (h_av.^(-1).*du_dz.*T_Delta_y_partial_cell);

%eq. 21: prescribed flux at cell centre at the bed
flux_u_bottom  = (gamma.*(1+gamma_pert).*f_slide_Tbed.*f_slide_hydro.*u_bed).*T_Delta_y_cell(T_bdy_nodes_bed); 
net_u_verflux(T_bdy_nodes_bed) = net_u_verflux(T_bdy_nodes_bed) - flux_u_bottom;

div_u_fluxes =   net_u_verflux./(h_av.*T_Delta_z_cell_volume.*T_Delta_y_cell)+ net_u_horflux./T_Delta_y_cell; 
%source term
S_u = (h + b - h_prev - b_prev)/(dx).* ones(T_nodes, 1);

%conservation law
fout(u_index) = div_u_fluxes - S_u;

%% P'

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

%dirichlet condition
fout(p_index(1)) = p(1);

%% CONSTITUTIVE RELATIONS for p(p'), v, w_eff

%notation is p -> p', pr->p, eq. 33
pr = p - (u-u_prev)./dx_int;

%transverse velocity, eq. 34
dpsi_dz_full = [sparse(length(psi_bdy_nodes_top),1); dpsi_dz; sparse(length(psi_bdy_nodes_bed),1)];
v = dphi_dy + h_prev^(-1).*dpsi_dz_full; 

%vertical velocity (defined at T hor cell edges), eq. 35
[dpsi_dy_reindexed,~] = verticalvelocity(parameters,dpsi_dy,Ddiscvar.ddpsidy_dpsi);
w_eff = -dpsi_dy_reindexed + h_prev.^(-1).*dphi_dz - T_coor_hor_cell_edges_z./2.*(u_vert *(h-h_prev)/dx + u_vert_prev .*(h_prev-h_pprev)/dx_prev);

%% HEAT EQUATION, ICE
%diffusive fluxes
net_T_diffhorflux = T_connect_hor * dT_dy;
net_T_diffverflux = T_connect_ver*( dT_dz.*T_Delta_y_partial_cell);

%advective fluxes 
T_adv_flux_ver = w_eff.*T_vert;
T_adv_flux_transverse = v.*T_hor.*h_prev;

net_T_advfluxver = T_connect_ver*(T_adv_flux_ver.*T_Delta_y_partial_cell);
net_T_advfluxtransverse = T_connect_hor *T_adv_flux_transverse;

if parameters.flag.plug == 0
    net_T_advfluxalong = (u.*h_prev.*T -u_prev.*h_pprev.*T_prev)./dx_int;
elseif parameters.flag.plug == 1
    net_T_advfluxalong = (u_prev_T.* h_prev_T.*T -u_prev_T.* h_prev_T.*T_prev)./dx_int;
end 

%source term
du_dz_centre_top = sparse(length(T_bdy_nodes_top),1);
du_dz_centre_bottom = h*gamma.*(1+gamma_pert).*f_slide_Tbed.*f_slide_hydro.*u_bed;
du_dz_centre_full = [du_dz_centre_top; du_dz_centre; du_dz_centre_bottom];

S_T = (2*h_prev)^(-2).*((du_dz_centre_full).^2 + (du_dz_centre_full_prev).^2 + 2.*du_dz_centre_full.*du_dz_centre_full_prev)+...
    1/4*((du_dy_centre).^2 + (du_dy_centre_prev).^2 + 2.*du_dy_centre.*du_dy_centre_prev);

%conservation law

if parameters.flag.heat_full == 1 || isfield(parameters,'flag') == 0
fout(T_index) = ...
    Pe/h_prev.*(net_T_advfluxtransverse./T_Delta_y_cell + net_T_advfluxalong + net_T_advfluxver./(T_Delta_z_cell.*T_Delta_y_cell))- ...
    net_T_diffverflux./(h_prev.^2*T_Delta_z_cell.*T_Delta_y_cell)- net_T_diffhorflux./T_Delta_y_cell - alpha*S_T;
elseif parameters.flag.heat_full == 0 || parameters.flag.plug == 0
    fout(T_index) = ...
    Pe/h_prev.*( net_T_advfluxalong)- ...
    net_T_diffverflux./(h_prev.^2*T_Delta_z_cell.*T_Delta_y_cell)- net_T_diffhorflux./T_Delta_y_cell;
end

%surface and bed boundary condition
fout(T_index(T_bdy_nodes_top)) = T(T_bdy_nodes_top) - T_s;
%fout(2*psi_nodes+3*T_nodes+T_bdy_nodes_bed) = T(T_bdy_nodes_bed) - T_bed;

%% HEAT EQUATION, BED
net_Tb_diffhorflux =  Tb_connect_hor *( -dTb_dy);
net_Tb_diffverflux = Tb_connect_ver *( -dTb_dz.*Tb_Delta_y_partial_cell);

%neumann conditions at top and bottom
bedflux_below = - (T(T_bdy_nodes_bed) - Tb(Tb_bdy_nodes_top)).*Tb_Delta_y_cell(Tb_bdy_nodes_top)./(Tb_Delta_z_cell(Tb_bdy_nodes_top));   
net_Tb_diffverflux(Tb_bdy_nodes_top) = net_Tb_diffverflux(Tb_bdy_nodes_top) + bedflux_below;

net_Tb_diffverflux(Tb_bdy_nodes_bed) = net_Tb_diffverflux(Tb_bdy_nodes_bed) - (nu).*Tb_Delta_y_cell(Tb_bdy_nodes_bed);

%conservation law
fout(Tb_index) = ...
    net_Tb_diffverflux./(Tb_Delta_z_cell.*Tb_Delta_y_cell)+ net_Tb_diffhorflux./Tb_Delta_y_cell;

%% DRAINAGE 

%defined modified effective pressure (this is an auxiliary variable, living
%at T_nodes)

%enforce frozen/temperate bed with indicator function from previous
%timestep, eq. 42-43
index_temperate_prev = find(I_prev == 1); 
index_frozen_prev = find(I_prev == 0); 

%nb: the two equations below all together take n_rows =
%length(T_bdy_nodes_bed);

fout(Pi_index(index_frozen_prev)) = Pi(index_frozen_prev);
fout(T_index(T_bdy_nodes_bed(index_temperate_prev))) = T(T_bdy_nodes_bed(index_temperate_prev));

%bed state indicator at current time step, eq. 47
I = sparse(length(T_bdy_nodes_bed),1);
index_temperate_T = find(T_bed>0);
index_incipient = find(T_bed==0);
index_temperate_P = index_incipient(Pi(index_incipient)>0);

I(index_temperate_T) = ones(length(index_temperate_T),1);
I(index_temperate_P) = ones(length(index_temperate_P),1);

%normal stress at the base
%second order accurate, one-sided interpolation formula for second
%derivative of phi at the base
ddpsy_dz_dy = (flux_psi_bottom(circshift(1:length(psi_bdy_nodes_bed),-1,2)) - flux_psi_bottom(1:end))./psi_Delta_y_cell(psi_bdy_nodes_bed);

%spectral computation for phi_zz

if isfield(parameters.flag, 'sigma') == 0 ||  strcmp(parameters.flag.sigma, 'FV') == 1
    ddphi_dz_bed = 2./(2*(T_Delta_z_cell(T_bdy_nodes_bed)).^2).*(-2*phi(T_bdy_nodes_bed-T_nodes_hor(end))+ phi(T_bdy_nodes_bed-2*T_nodes_hor(end)) +phi(T_bdy_nodes_bed));
elseif strcmp(parameters.flag.sigma, 'spectral') == 1
    fout_phizz = conv_phizz(length(T_bdy_nodes_bed), parameters.grid.N.extra.bd_y, parameters.grid.N.coor_nodes.y(T_bdy_nodes_bed), h_prev);
    Gf_conv = fout_phizz.Gf_conv;
    Gfx_conv = fout_phizz.Gfx_conv;
    
    %define f, f_x and u_x; h should be h_i
    f = (gamma.*(1+gamma_pert).*f_slide_Tbed.*f_slide_hydro.*u_bed);
    f_x = (f-f_prev)./ dx_prev;
    u_x = (u_bed-u_prev(T_bdy_nodes_bed))./dx;
    h_x = (h-h_prev)./dx;
    
    ddphi_dz_bed = Gf_conv*(f.*h_x) + Gfx_conv*f_x + 1/parameters.grid.N.extra.bd_y*sum(u_x.*T_Delta_y_cell(T_bdy_nodes_bed));
end
    
sigma = pr(T_bdy_nodes_bed)  - 2*h_prev^(-2).*ddphi_dz_bed + 2*h_prev^(-1).*ddpsy_dz_dy;

%meltwater fluxes
%permeabilities
[kappa_fout, ~] = permeabilities(Pi, parameters);
kappa = kappa_fout.kappa;
kappa_2 = kappa_fout.kappa_2;

qx = -kappa.*(h-h_prev+r^(-1)*(b-b_prev))/dx;

index_sigma_edge_down = circshift(1:length(sigma),-1, 2);
index_sigma_edge_up = 1:length(sigma);

if isfield(parameters.flag, 'flux_y') == 0
    dsigma_dy = (sigma(index_sigma_edge_down)-sigma(index_sigma_edge_up))./T_Delta_y_cell(T_bdy_nodes_bed);
elseif strcmp(parameters.flag.flux_y, 'upwind') == 1
    dsigma_dy = parameters.v_in_prev.dsigma_dy_prev;  
else
    dsigma_dy = (sigma(index_sigma_edge_down)-sigma(index_sigma_edge_up))./T_Delta_y_cell(T_bdy_nodes_bed);
end
dPi_dy = (Pi(index_sigma_edge_down)-Pi(index_sigma_edge_up))./T_Delta_y_cell(T_bdy_nodes_bed);

%define I_prev at lateral cell edges
I_prev_edge = 1/2*(I_prev(index_sigma_edge_up) + I_prev(index_sigma_edge_down));
I_prev_edge(I_prev_edge<1) = 0;

if isfield(parameters.flag, 'flux_y') == 0 %pick scheme for transverse permeability in kappa d sigma/dy
    
    kappa_edge = kappa;
    kappa_2_edge = kappa_2;
    
elseif strcmp(parameters.flag.flux_y, 'centered') == 1
    Pi_edge = 1/2*(Pi(index_sigma_edge_up) + Pi(index_sigma_edge_down));
    [kappa_edge_fout, ~] = permeabilities(Pi_edge, parameters);
    kappa_edge = kappa_edge_fout.kappa;
    kappa_2_edge = kappa_edge_fout.kappa_2;
    
elseif strcmp(parameters.flag.flux_y, 'centered_cutoff') == 1
    Pi_edge = 1/2*(Pi(index_sigma_edge_up) + Pi(index_sigma_edge_down));
    [kappa_edge_fout, ~] = permeabilities(Pi_edge, parameters);
    kappa_edge = max(kappa_edge_fout.kappa,0);
    kappa_2_edge = kappa_edge_fout.kappa_2;
    
elseif strcmp(parameters.flag.flux_y, 'centered_harmonic') == 1
    
    Pi_up = Pi(index_sigma_edge_up);
    Pi_down = Pi(index_sigma_edge_down);
    [kappa_up_fout, ~] = permeabilities(Pi_up, parameters);
    [kappa_down_fout, ~] = permeabilities(Pi_down, parameters);
    
    kappa_up = max(kappa_up_fout.kappa,0);
    kappa_down =  max(kappa_down_fout.kappa,0);
    kappa_edge = 2.* kappa_up.*kappa_down./(kappa_up + kappa_down);
    kappa_edge(kappa_up.*kappa_down==0) = 0;
    kappa_2_edge = kappa_up_fout.kappa_2;
    
elseif strcmp(parameters.flag.flux_y, 'upwind') == 1
    
    adv_speed = -parameters.v_in_prev.dsigma_dy_prev; 
    index_p = find(adv_speed>=0);
    index_n = find(adv_speed<0);
    
    index_Pi_edge = zeros(length(index_p),1);
    index_Pi_edge(index_p) = index_sigma_edge_up(index_p);
    index_Pi_edge(index_n) = index_sigma_edge_down(index_n);
    
    Pi_edge = Pi(index_Pi_edge);
    [kappa_edge_fout, ~] = permeabilities(Pi_edge, parameters);
    
    kappa_edge = kappa_edge_fout.kappa;
    kappa_2_edge = kappa_edge_fout.kappa_2;
       
end


if isfield(parameters.flag, 'fluxy_sigma') == 0
    qy = -I_prev_edge.*(kappa_edge.*dsigma_dy + beta.*kappa_2_edge.*dPi_dy);
elseif strcmp(parameters.flag.fluxy_sigma, 'full')== 1
    qy = -I_prev_edge.*(kappa_edge.*dsigma_dy + beta.*kappa_2_edge.*dPi_dy);
elseif strcmp(parameters.flag.fluxy_sigma, 'Pionly')== 1
    qy = -I_prev_edge.*(beta.*kappa_2_edge.*dPi_dy);
end


%% BASAL ENERGY BUDGET
%this is effectively a constraint on basal temperature, and should be
%treated as an additional equation in the unknown T_bed

bedflux_above = - (T(T_bdy_nodes_bed-T_nodes_hor(end))- T(T_bdy_nodes_bed))./(h_prev*T_Delta_z_cell(T_bdy_nodes_bed));   
bedflux_below = - (T(T_bdy_nodes_bed) - Tb(Tb_bdy_nodes_top))./(Tb_Delta_z_cell(Tb_bdy_nodes_top));   

m = -(bedflux_above - bedflux_below) + alpha*(gamma.*(1+gamma_pert).*f_slide_Tbed.*f_slide_hydro).*u(T_bdy_nodes_bed).^(2);

index_sigma_down = 1:length(sigma);
index_sigma_up = circshift(1:length(sigma),1,2);
div_qy = (qy(index_sigma_down)-qy(index_sigma_up))./(T_Delta_y_cell(T_bdy_nodes_bed));
enth_eq = I_prev.*(qx-qx_prev)/dx_int  - m + div_qy;

%constrain Pi where bed was temperate at previous timestep and Tbed where
%the bed was frozen
fout(Pi_index(index_temperate_prev)) = enth_eq(index_temperate_prev);
fout(T_index(T_bdy_nodes_bed(index_frozen_prev))) = enth_eq(index_frozen_prev);

%% AUXILIARY VARIABLES
faux.du_dz_centre_full = du_dz_centre_full;
faux.du_dy_centre = du_dy_centre;
faux.u_vert = u_vert;
faux.I = I;
faux.qx= qx;
faux.qy = qy;
faux.sigma = sigma;
faux.dsigma_dy = (sigma(index_sigma_edge_down)-sigma(index_sigma_edge_up))./T_Delta_y_cell(T_bdy_nodes_bed);

faux.sigma_phi =  - 2*h_prev^(-2).*ddphi_dz_bed;
faux.sigma_psi = + 2*h_prev^(-1).*ddpsy_dz_dy;
faux.sigma_p = pr(T_bdy_nodes_bed);
faux.dsigmaphi_dy = (faux.sigma_phi(index_sigma_edge_down)-faux.sigma_phi(index_sigma_edge_up))./T_Delta_y_cell(T_bdy_nodes_bed);
faux.dsigmapsi_dy = (faux.sigma_psi(index_sigma_edge_down)-faux.sigma_psi(index_sigma_edge_up))./T_Delta_y_cell(T_bdy_nodes_bed);

faux.kappa = kappa_edge;
faux.kappa_edge = kappa_edge;
faux.qy_sigma = -I_prev_edge.*(kappa_edge.*dsigma_dy);
faux.qy_Pi = -I_prev_edge.*(beta.*kappa_2_edge.*dPi_dy);
faux.div_qy_sigma = (faux.qy_sigma(index_sigma_down)-faux.qy_sigma(index_sigma_up))./(T_Delta_y_cell(T_bdy_nodes_bed));


faux.h_prev = h_prev;
faux.u_b = u(T_bdy_nodes_bed);
faux.tau_b = flux_u_bottom;
faux.bedflux = bedflux_above;
faux.bedflux_bed = bedflux_below;
faux.friction = alpha*(gamma.*(1+gamma_pert).*f_slide_Tbed.*f_slide_hydro).*u(T_bdy_nodes_bed).^(2);
faux.surfaceslope = (h + b - h_prev - b_prev)/(dx);
faux.m = m;
faux.S = S_T;
faux.w =  h_prev.^(-1).*dphi_dz;
faux.wtop = -a + 1/2.*(u(1).*(h - h_prev)/dx + (h_prev - h_pprev)/dx_prev*u_prev(1));
faux.h_av = h_av;

%first order correction to surface elevation
ddphi_dz_top = (phi(T_bdy_nodes_top) -2* phi(T_bdy_nodes_top + length(T_bdy_nodes_top))+ phi(T_bdy_nodes_top + 2*length(T_bdy_nodes_top)))./(T_Delta_z_cell(T_bdy_nodes_top).^2) ; 
ddpsi_dz_dy_top = (flux_psi_top(circshift(1:length(psi_bdy_nodes_top),-1,2)) - flux_psi_top(1:end))./psi_Delta_y_cell(psi_bdy_nodes_top);
faux.s1 = pr(T_bdy_nodes_top) + ddpsi_dz_dy_top - h_prev^(-2).*ddphi_dz_top;

faux.v = v;
faux.v(T_bdy_nodes_bed) = faux.v(T_bdy_nodes_bed) + flux_psi_bottom;
faux.v(T_bdy_nodes_top) = faux.v(T_bdy_nodes_top) + flux_psi_top;
faux.w = w_eff;
faux.v_bed = faux.v(T_bdy_nodes_bed);
faux.v_top = faux.v(T_bdy_nodes_top);
faux.f = (gamma.*(1+gamma_pert).*f_slide_Tbed.*f_slide_hydro.*u_bed);
faux.ddphi_dzeta_bed = ddphi_dz_bed;














