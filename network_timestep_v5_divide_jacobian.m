function [fout] = network_timestep_v5_divide_jacobian(v_in,parameters)
%Jacobian of network_timestep_v5_divide.m
%Tested against analytical jacobian. Elisa Mantelli, 22 Oct 2019

%parameters 
gamma = parameters.gamma;
Pe = parameters.Pe;
alpha = parameters.alpha;

%psi grid
psi_nodes = parameters.grid.psi.n_nodes.tot;                                   %number of nodes
psi_bdy_nodes_bed = parameters.grid.psi.bdy_nodes.bed;


%T grid
T_nodes = parameters.grid.N.n_nodes.tot;                                   %number of nodes
T_nodes_hor = parameters.grid.N.n_nodes.hor;                               %number of nodes
T_nedges_ver = parameters.grid.N.n_edges.vert;

T_up_node_ver = parameters.grid.N.up_node.vert;                            %list (n_edges-by-1 vector) of 'upstream' node for each vertical edge
T_down_node_ver = parameters.grid.N.down_node.vert;                        %list (n_edges-by-1 vector) of 'downstream' nodes for each vertical edge

T_bdy_nodes_top = parameters.grid.N.bdy_nodes.top;                         %list of nodes adjacent to the top of the box;
T_bdy_nodes_bed = parameters.grid.N.bdy_nodes.bottom;

T_Delta_z_cell = parameters.grid.N.Delta_z_cell;                           %length of cells, ver (list)
T_Delta_z_cell_volume = parameters.grid.N.Delta_z_cell_volume;             %length of cells, ver (list)
T_Delta_y_cell = parameters.grid.N.Delta_y_cell;                           %length of cells, hor (list)
T_Delta_y_partial_cell = parameters.grid.N.delta_y_vflux;                  %length of hor cell esges crossed by vertical network edges, hor (list)
T_coor_hor_cell_edges_z = parameters.grid.N.coor_hor_cell_edges.z;

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
fout = sparse(length(v_in),length(v_in));

%%  PRELIMINARIES 

%DISCRETE DEP. VARIABLES
v_in_disc = [sparse(2*psi_nodes,1); phi; u; sparse(T_nodes,1); T; sparse(Tb_nodes+1,1); h; Q; sparse(length(T_bdy_nodes_bed),1)];
[discvar, Ddiscvar] = discretisation_hanging_full_v4(v_in_disc, parameters);

u_vert = discvar.u_vert;
duvert_du = Ddiscvar.duvert_du;

du_dz_centre = discvar.du_dz_centre;   % derivative at u cell centre; note that the indexing is the same as for u_vert_centre
ddudzcentre_du = Ddiscvar.d_dudzcentre_du;

du_dz = discvar.du_dz; 
ddudz_du =  Ddiscvar.ddudz_du;

dphi_dz = discvar.dphi_dz;
ddphidz_dphi = Ddiscvar.ddphidz_dphi;

ddTdz_dT = Ddiscvar.ddTdz_dT;

T_vert = discvar.T_vert;
dTvert_dT = Ddiscvar.dTver_dT;

%construct regularized bedwater content
T_bed = T(T_bdy_nodes_bed);
u_bed = u(T_bdy_nodes_bed);

%TIMESTEPPING
x_current = parameters.timestep.x; %i+1
x_prev = 0; %x_i
x_pprev = -x_current; %x_i-1

dx_prev = x_prev-x_pprev;  % x_i - x_{i-1}
dx_int = parameters.timestep.dx_prev; % x_{i+1/2} - x_{i-1/2}
dx = x_current -x_prev; % x_{i+1} - x_{i}

h_av = (h*parameters.timestep.dx_prev+h_prev*parameters.timestep.dx)./(parameters.timestep.dx_prev + parameters.timestep.dx); %eq. 16

dhav_dh = parameters.timestep.dx_prev./(parameters.timestep.dx_prev + parameters.timestep.dx)*ones(length(h),1);

%construct regularized bedwater content
if parameters.flag_Tdep == 1
    [f_slide_Tbed, DfslideTbed] = regularizedfriction_temperature(T_bed, parameters);
else
    f_slide_Tbed = ones(length(T_bdy_nodes_bed),1);
    DfslideTbed = zeros(length(T_bdy_nodes_bed),1);
end

if isfield(parameters,'gamma_pert')==1
    gamma_pert = parameters.gamma_pert;
elseif isfield(parameters,'gamma_pert')==0
    gamma_pert = zeros(length(T_bdy_nodes_bed),1);
end

[f_slide_hydro, Df_slide_hydro] = regularizedfriction_hydrology(u_bed,zeros(length(u_bed),1), parameters);
Df_slide_hydro_dub = Df_slide_hydro.dub;

%% SYMMETRY ACROSS THE DIVIDE

h_av_prev = h_av;
dhavprev_dhav = 1;

u_prev = -u;  
duprev_du = -1;

dQprev_dQ = -1;

T_prev = T; 
dTprev_dT = 1;

h_pprev = h; 
dhpprev_dh = 1;

u_vert_prev = -u_vert; 
duvertprev_duvert = -1;

%% MASS CONSERVATION AND EQ. FOR h, eqs. 16-18

%mass conservation
%fout(Q_index) = (Q- Q_prev)./(dx_int)- (a*parameters.grid.N.extra.bd_y);
fout(Q_index,Q_index) = (1 - dQprev_dQ)./dx_int;

%conservation law
%fout(2*psi_nodes+4*T_nodes+Tb_nodes+length(psi_bdy_nodes_bed)+1) = Q - h_av.*u_int;
u_int = sum(u.*T_Delta_z_cell_volume.*T_Delta_y_cell);
duint_du = (T_Delta_z_cell_volume.*T_Delta_y_cell).';

fout(h_index, Q_index) = 1;
fout(h_index, u_index) = -h_av.*duint_du;
fout(h_index, h_index) = -dhav_dh.*u_int;

 
%% SLIDING LAW
%at T cell centres
dubed_du = ones(length(T_bdy_nodes_bed),1);

%% PHI, eq. 22 with bcs 26, 27

%diffusive fluxes
dphiverfluxdiff_dphi = T_connect_ver * (spdiags(h_prev^(-1).*T_Delta_y_partial_cell,0,T_nedges_ver, T_nedges_ver)* ddphidz_dphi);

%advective flux
dfluxphiadv_du = -1/2*(h - h_prev)/dx*sparse(1:length(u_vert), 1:length(u_vert),T_coor_hor_cell_edges_z.*T_Delta_y_partial_cell,length(u_vert),length(u_vert)) *duvert_du; 
dfluxphiadv_dh = -1/2*T_coor_hor_cell_edges_z.*u_vert.*T_Delta_y_partial_cell./dx;

dnetfluxphiadv_du = T_connect_ver * dfluxphiadv_du;
dnetfluxphiadv_dh = T_connect_ver *dfluxphiadv_dh;

%net_phi_verflux =  net_phi_diff_verflux + net_phi_adv_verflux;
dnetphiverflux_dphi = dphiverfluxdiff_dphi;
dnetphiverflux_du = dnetfluxphiadv_du;
dnetphiverflux_dh = dnetfluxphiadv_dh;

%S_phi = (h_av*u - h_av_prev*u_prev)./dx_int; 
dSphi_du = sparse(1:T_nodes,1:T_nodes, (h_av - h_av_prev*duprev_du)/dx_int*ones(length(T_nodes),1), T_nodes,T_nodes);
dSphi_dh = sparse(1:T_nodes,ones(T_nodes,1), dhav_dh/dx_int*u - dhavprev_dhav* dhav_dh*u_prev./dx_int, T_nodes,1);

%conservation law 
%div_phi_fluxes = net_phi_verflux./(h_prev*T_Delta_z_cell_volume.*T_Delta_y_cell)+ net_phi_horflux./T_Delta_y_cell;
%fout(2*psi_nodes+1:2*psi_nodes+T_nodes) = div_phi_fluxes + h_prev.^(-1)*S_phi;
fout(phi_index,phi_index) = spdiags(1./(h_prev.*T_Delta_z_cell_volume.*T_Delta_y_cell),0,T_nodes,T_nodes)*dnetphiverflux_dphi;
fout(phi_index,u_index) = spdiags(1./(h_prev.*T_Delta_z_cell_volume.*T_Delta_y_cell),0,T_nodes,T_nodes)*dnetphiverflux_du + h_prev.^(-1)*dSphi_du;
fout(phi_index,h_index) = spdiags(1./(h_prev.*T_Delta_z_cell_volume.*T_Delta_y_cell),0,T_nodes,T_nodes)*dnetphiverflux_dh + h_prev.^(-1)*dSphi_dh;

fout(phi_index(T_bdy_nodes_top), :) = sparse(1:length(T_bdy_nodes_top), phi_index(T_bdy_nodes_top), ones(length(T_bdy_nodes_top),1), length(T_bdy_nodes_top), length(v_in));
%% ALONG FLOW VELOCITY, U, eq. 19 with 20-21
%diffusive fluxes
net_u_verflux = T_connect_ver *( h_av.^(-1).*du_dz.*T_Delta_y_partial_cell);
duverflux_du = T_connect_ver * (spdiags(h_av^(-1).*T_Delta_y_partial_cell,0,T_nedges_ver, T_nedges_ver)* ddudz_du);
duverflux_dh = T_connect_ver * (-h_av^(-2).*T_Delta_y_partial_cell.*du_dz .*dhav_dh);

%eq. 21: prescribed flux at cell centre at the bed
% flux_u_bottom  = (gamma.*f_slide_Tbed.*u_bed).*T_Delta_y_cell(T_bdy_nodes_bed); 
% net_u_verflux(T_bdy_nodes_bed) = net_u_verflux(T_bdy_nodes_bed) - flux_u_bottom;

flux_u_bottom  = (gamma.*(1+gamma_pert).*f_slide_Tbed.*f_slide_hydro.*u_bed).*T_Delta_y_cell(T_bdy_nodes_bed); 
dfluxubottom_du = sparse(T_bdy_nodes_bed, T_bdy_nodes_bed,gamma.*(1+gamma_pert).*f_slide_Tbed.*f_slide_hydro.*T_Delta_y_cell(T_bdy_nodes_bed), T_nodes, T_nodes) + ...
    sparse(T_bdy_nodes_bed, T_bdy_nodes_bed,gamma.*(1+gamma_pert).*f_slide_Tbed.*Df_slide_hydro_dub.*u_bed.*T_Delta_y_cell(T_bdy_nodes_bed), T_nodes, T_nodes);
dfluxubottom_dT = sparse(T_bdy_nodes_bed, T_bdy_nodes_bed, gamma.*(1+gamma_pert).*f_slide_hydro.*u_bed.*DfslideTbed.*T_Delta_y_cell(T_bdy_nodes_bed),T_nodes,T_nodes);

net_u_verflux(T_bdy_nodes_bed) = net_u_verflux(T_bdy_nodes_bed) - flux_u_bottom;
duverflux_du = duverflux_du - dfluxubottom_du;
duverflux_dT = -dfluxubottom_dT;

%source term
%S_u = (h + b - h_prev - b_prev)/(dx).* ones(T_nodes, 1);
dSu_dh = (1)./(dx).* ones(T_nodes, 1);

%conservation law
%fout(2*psi_nodes+T_nodes+1:2*psi_nodes+2*T_nodes) = div_u_fluxes - S_u;
%div_u_fluxes =   net_u_verflux./(h_av.*T_Delta_z_cell.*T_Delta_y_cell)+ net_u_horflux./T_Delta_y_cell; %nb: h_int depends on h, so must be differentiated!
fout(u_index,u_index) = spdiags(1./(h_av.*T_Delta_z_cell_volume.*T_Delta_y_cell),0,T_nodes,T_nodes)*duverflux_du;
fout(u_index,h_index) = spdiags(1./(h_av.*T_Delta_z_cell_volume.*T_Delta_y_cell),0,T_nodes,T_nodes)*duverflux_dh - net_u_verflux./(h_av.^2.*T_Delta_z_cell_volume.*T_Delta_y_cell)*dhav_dh - dSu_dh;
fout(u_index,T_index) = spdiags(1./(h_av.*T_Delta_z_cell_volume.*T_Delta_y_cell),0,T_nodes,T_nodes)*duverflux_dT;

%% CONSTITUTIVE RELATION for  w_eff, 
w_eff =  h_prev.^(-1).*dphi_dz - T_coor_hor_cell_edges_z./2.*(u_vert *(h-h_prev)/dx + u_vert_prev .*(h_prev-h_pprev)/dx_prev); 
dweff_dh = - T_coor_hor_cell_edges_z./2.*(u_vert/dx) - T_coor_hor_cell_edges_z./2.*(u_vert_prev .*(-dhpprev_dh)/dx_prev);
dweff_ddphidz = h_prev.^(-1);
dweff_duvert = - T_coor_hor_cell_edges_z./2.*(h-h_prev)/dx - T_coor_hor_cell_edges_z./2.*(duvertprev_duvert.*(h_prev-h_pprev)/dx_prev);

%% HEAT EQUATION, ICE, eq. 36-37 with bcs. 39 and 41 (13/10)
%diffusive fluxes
%net_T_diffverflux = T_connect_ver *( dT_dz.*T_Delta_y_partial_cell);
dTdiffverflux_dT =  T_connect_ver * (spdiags(T_Delta_y_partial_cell,0,T_nedges_ver, T_nedges_ver)* ddTdz_dT);
%advective fluxes 

%vertical
%T_adv_flux_ver = w_eff.*T_vert;
dTadverflux_ddphidz = dweff_ddphidz.*T_vert;
dTadverflux_du = spdiags(T_Delta_y_partial_cell.*dweff_duvert.*T_vert,0, length(T_vert), length(T_vert))*duvert_du;
dTadverflux_dTver =  w_eff;
dTadverflux_dh =  dweff_dh.*T_vert;

%net_T_advfluxver = T_connect_ver*(T_adv_flux_ver.*T_Delta_y_partial_cell);
dnetTadvverflux_dT = T_connect_ver * ( spdiags(T_Delta_y_partial_cell.*dTadverflux_dTver ,0,T_nedges_ver, T_nedges_ver)*dTvert_dT);
dnetTadvverflux_dh =  T_connect_ver * (T_Delta_y_partial_cell.*dTadverflux_dh);
dnetTadvverflux_dphi = T_connect_ver * ( spdiags(T_Delta_y_partial_cell.*dTadverflux_ddphidz ,0,T_nedges_ver, T_nedges_ver)*ddphidz_dphi);
dnetTadvverflux_du = T_connect_ver * dTadverflux_du;

%downstream
%net_T_advfluxalong = (u.*h_prev.*T -u_prev.*h_pprev.*T_prev)./dx_int;
dnetTadvalongflux_dT = sparse(1:T_nodes, 1:T_nodes,u.*h_prev./dx_int,T_nodes,T_nodes) - sparse(1:T_nodes, 1:T_nodes,u_prev.*h_pprev.*dTprev_dT./dx_int,T_nodes,T_nodes);
dnetTadvalongflux_du = sparse(1:T_nodes, 1:T_nodes,T.*h_prev./dx_int,T_nodes,T_nodes) - sparse(1:T_nodes, 1:T_nodes,T_prev.*h_pprev.*duprev_du./dx_int,T_nodes,T_nodes);
dnetTadvalongflux_dh = - T_prev.*dhpprev_dh.*u_prev./dx_int;

%source term
%S_T =  (h_prev)^(-2).*((du_dz_centre_full).^2);
    %source term

du_dz_centre_top = sparse(length(T_bdy_nodes_top),1);
ddudzcentretop_du = sparse(length(T_bdy_nodes_top),T_nodes);
du_dz_centre_bottom =  h*gamma.*(1+gamma_pert).*f_slide_Tbed.*f_slide_hydro.*u_bed;
ddudzcentrebottom_du = sparse(1:length(T_bdy_nodes_bed), T_bdy_nodes_bed, h*gamma.*(1+gamma_pert).*f_slide_Tbed.*f_slide_hydro.*dubed_du, length(T_bdy_nodes_bed), T_nodes) + ...
    sparse(1:length(T_bdy_nodes_bed), T_bdy_nodes_bed, h*gamma.*(1+gamma_pert).*f_slide_Tbed.*u_bed.*Df_slide_hydro_dub.*dubed_du, length(T_bdy_nodes_bed), T_nodes); 
ddudzcentrebottom_dh = sparse(T_bdy_nodes_bed, ones(length(T_bdy_nodes_bed),1), gamma.*(1+gamma_pert).*f_slide_Tbed.*f_slide_hydro.*u_bed, T_nodes, 1);
ddudzcentrebottom_dT = sparse(T_bdy_nodes_bed, T_bdy_nodes_bed, h*gamma.*(1+gamma_pert).*f_slide_hydro.*u_bed.*DfslideTbed, T_nodes,T_nodes);

du_dz_centre_full = [du_dz_centre_top; du_dz_centre; du_dz_centre_bottom];
ddudzcentrefull_du = [ddudzcentretop_du; ddudzcentre_du; ddudzcentrebottom_du];
ddudzcentrefull_dh = ddudzcentrebottom_dh;

dST_du = (sparse(1:T_nodes, 1:T_nodes,(h_prev)^(-2).*2*du_dz_centre_full,T_nodes,T_nodes))*ddudzcentrefull_du;
dST_dh = (sparse(1:T_nodes,1:T_nodes,(h_prev)^(-2).*2*du_dz_centre_full,T_nodes,T_nodes))*ddudzcentrefull_dh;
%dST_dT = (spdiags((2*h_prev)^(-2).*2*du_dz_centre_full,0,T_nodes) + spdiags((2*h_prev)^(-2).*2*du_dz_centre_full_prev,0,T_nodes))*ddudzcentrefull_dTbed;

%conservation law
%Pe/h_prev.*( net_T_advfluxalong + net_T_advfluxver./(T_Delta_z_cell.*T_Delta_y_cell))- ...
%    net_T_diffverflux./(h_prev.^2*T_Delta_z_cell.*T_Delta_y_cell) - alpha*S_T;

fout(T_index,T_index) = Pe/h_prev*(sparse(1:T_nodes,1:T_nodes,1./(T_Delta_z_cell.*T_Delta_y_cell), T_nodes,T_nodes)*dnetTadvverflux_dT + ...
    + dnetTadvalongflux_dT ) - ...
    sparse(1:T_nodes,1:T_nodes,1./(h_prev.^2*T_Delta_z_cell.*T_Delta_y_cell), T_nodes, T_nodes)*dTdiffverflux_dT;% - alpha*dST_dT;
fout(T_index,u_index) = Pe/h_prev*(sparse(1:T_nodes, 1:T_nodes, 1./(T_Delta_z_cell.*T_Delta_y_cell), T_nodes,T_nodes)*dnetTadvverflux_du + dnetTadvalongflux_du )- alpha*dST_du;
fout(T_index,phi_index) = Pe/h_prev*(sparse(1:T_nodes,1:T_nodes,1./(T_Delta_z_cell.*T_Delta_y_cell), T_nodes,T_nodes)*dnetTadvverflux_dphi);
fout(T_index,h_index) = Pe/h_prev*(sparse(1:T_nodes,1:T_nodes,1./(T_Delta_z_cell.*T_Delta_y_cell), T_nodes,T_nodes)*dnetTadvverflux_dh + dnetTadvalongflux_dh)- alpha*dST_dh;

%surface and bed boundary condition
%fout(2*psi_nodes+3*T_nodes+T_bdy_nodes_top) = T(T_bdy_nodes_top) - T_s;
fout(T_index(T_bdy_nodes_top),:) = sparse(length(T_bdy_nodes_top),size(fout,2));
fout(T_index(T_bdy_nodes_top),T_index) = sparse(1:length(T_bdy_nodes_top),T_bdy_nodes_top,ones(length(T_bdy_nodes_top),1), length(T_bdy_nodes_top),T_nodes);

fout(T_index(T_bdy_nodes_bed),:) = sparse(length(T_bdy_nodes_bed),length(v_in));

%% BASAL ENERGY BUDGET, eq. 41 (22/10)
%this is effectively a constraint on basal temperature, and should be
%treated as an additional equation in the unknown T_bed

%bedflux_above = - (T(T_bdy_nodes_bed-T_nodes_hor(end))- T(T_bdy_nodes_bed))./(h_prev*T_Delta_z_cell(T_bdy_nodes_bed));   
dbedfluxabove_dT = sparse(1:length(T_bdy_nodes_bed),T_bdy_nodes_bed-T_nodes_hor(end), -1./(h_prev*T_Delta_z_cell(T_bdy_nodes_bed)), length(T_bdy_nodes_bed), T_nodes)+ ...
    sparse(1:length(T_bdy_nodes_bed),T_bdy_nodes_bed, 1./(h_prev*T_Delta_z_cell(T_bdy_nodes_bed)), length(T_bdy_nodes_bed), T_nodes);

%m = -(bedflux_above - bedflux_below) + alpha*(gamma.*(1+gamma_pert).*f_slide_Tbed.*f_slide_hydro).*u(T_bdy_nodes_bed).^(2);
%friction = alpha*(gamma.*(1+gamma_pert).*f_slide_Tbed.*f_slide_hydro).*u(T_bdy_nodes_bed).^(2);

dfriction_dT = sparse(1:length(T_bdy_nodes_bed),T_bdy_nodes_bed,  alpha*(gamma.*(1+gamma_pert)).*f_slide_hydro.*u(T_bdy_nodes_bed).^(2).*DfslideTbed,length(T_bdy_nodes_bed),T_nodes);
dfriction_du = sparse(1:length(T_bdy_nodes_bed),T_bdy_nodes_bed,alpha*(gamma.*(1+gamma_pert).*f_slide_hydro.*f_slide_Tbed).*2.*u(T_bdy_nodes_bed) + ...
    alpha*(gamma.*(1+gamma_pert).*f_slide_Tbed.*Df_slide_hydro_dub).*u(T_bdy_nodes_bed).^(2),length(T_bdy_nodes_bed),T_nodes);

dm_du = dfriction_du;
dm_dT = -(dbedfluxabove_dT) + dfriction_dT;

%conservation law
%enth_eq = (qx-qx_prev)/dx_int + div_qy - m;
dentheq_du =   - dm_du;%
dentheq_dT = -dm_dT;

%constrain T_bed where bed was frozen at current timestep
%fout(T_index(T_bdy_nodes_bed(index_frozen_prev))) = enth_eq(index_frozen_prev);
fout(T_index(T_bdy_nodes_bed), u_index) = dentheq_du;
fout(T_index(T_bdy_nodes_bed), T_index) = dentheq_dT;

































