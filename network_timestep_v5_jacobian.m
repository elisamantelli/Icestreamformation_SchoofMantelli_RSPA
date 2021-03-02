function [fout] = network_timestep_v5_jacobian(v_in,parameters)
%Analytical jacobian od network_timestep_v5
%Tested against finite difference jacobian. Elisa Mantelli, Jan 2021

%parameters 
gamma = parameters.gamma;
Pe = parameters.Pe;
alpha = parameters.alpha;
r = parameters.r;
beta = parameters.beta;

%psi grid
psi_nodes = parameters.grid.psi.n_nodes.tot;                                %number of nodes
psi_nedges_ver = parameters.grid.psi.n_edges.vert;

psi_nodes_hor = parameters.grid.psi.n_nodes.hor;                           %number of nodes

psi_bdy_nodes_top = parameters.grid.psi.bdy_nodes.top;                     %list of nodes adjacent to the top of the box;
psi_bdy_nodes_bed = parameters.grid.psi.bdy_nodes.bed;

psi_Delta_z_cell = parameters.grid.psi.Delta_z_cell;                       %length of cells, ver (list)
psi_Delta_y_cell = parameters.grid.psi.Delta_y_cell;                       %length of cells, hor (list)
psi_Delta_y_partial_cell = parameters.grid.psi.delta_y_vflux;              %length of hor cell esges crossed by vertical network edges, hor (list)

index_up_node_top_psi = 1:length(psi_bdy_nodes_top);
index_down_node_top_psi = circshift(1:length(psi_bdy_nodes_top), -1,2);

index_up_node_bed_psi = 1:length(psi_bdy_nodes_bed);
index_down_node_bed_psi = circshift(1:length(psi_bdy_nodes_bed), -1,2);

psi_connect_ver =  parameters.grid.psi.connect_ver;
psi_connect_hor = parameters.grid.psi.connect_hor;

%T grid
T_nodes = parameters.grid.N.n_nodes.tot;                                   %number of nodes
T_nodes_hor = parameters.grid.N.n_nodes.hor;                               %number of nodes
T_nedges_ver = parameters.grid.N.n_edges.vert;

T_bdy_nodes_top = parameters.grid.N.bdy_nodes.top;                         %list of nodes adjacent to the top of the box;
T_bdy_nodes_bed = parameters.grid.N.bdy_nodes.bottom;

T_Delta_z_cell = parameters.grid.N.Delta_z_cell;                           %length of cells, ver (list)
T_Delta_z_cell_volume = parameters.grid.N.Delta_z_cell_volume;             %length of cells, ver (list)
T_Delta_y_cell = parameters.grid.N.Delta_y_cell;                           %length of cells, hor (list)
T_Delta_y_partial_cell = parameters.grid.N.delta_y_vflux;                  %length of hor cell esges crossed by vertical network edges, hor (list)
T_coor_hor_cell_edges_z = parameters.grid.N.coor_hor_cell_edges.z;

index_up_node_top_T = circshift(1:length(T_bdy_nodes_top),1,2);
index_down_node_top_T = 1:length(T_bdy_nodes_top);

T_connect_ver =  parameters.grid.N.connect_ver;
T_connect_hor = parameters.grid.N.connect_hor;

%Tb grid
Tb_nodes = parameters.grid.Tb.n_nodes.tot;                                 %number of nodes
Tb_nedges_ver = parameters.grid.Tb.n_edges.vert;

Tb_bdy_nodes_top = parameters.grid.Tb.bdy_nodes.top;                       %list of nodes adjacent to the top of the box;

Tb_Delta_z_cell = parameters.grid.Tb.Delta_z_cell;                         %length of cells, ver (list)
Tb_Delta_y_cell = parameters.grid.Tb.Delta_y_cell;                         %length of cells, hor (list)
Tb_Delta_y_partial_cell = parameters.grid.Tb.delta_y_vflux;                %length of hor cell esges crossed by vertical network edges, hor (list)

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
phi = v_in(phi_index);
u = v_in(u_index);
p = v_in(p_index);
T = v_in(T_index);
h = v_in(h_index);
psi_ghost = v_in(psighost_index);
Pi = v_in(Pi_index);

%unpack variable at previous timestep
h_prev = parameters.v_in_prev.h_prev;
u_prev = parameters.v_in_prev.u_prev;

if parameters.flag.plug == 1
    u_prev_T = parameters.v_in_prev.u_prev_T;
    h_prev_T = parameters.v_in_prev.h_T;
end

h_pprev = parameters.v_in_prev.h_pprev;
du_dz_centre_full_prev = parameters.v_in_prev.du_dz_centre_full_prev;
du_dy_centre_prev = parameters.v_in_prev.du_dy_centre_prev;
u_vert_prev = parameters.v_in_prev.u_vert_prev;
I_prev = parameters.v_in_prev.I_prev;
if isfield(parameters.flag, 'sigma') == 1
    f_prev = parameters.v_in_prev.f_prev;
end

%initialize output
fout = sparse(length(v_in),length(v_in));

%%  PRELIMINARIES 

%DISCRETE DEP. VARIABLES
[discvar, Ddiscvar] = discretisation_hanging_full_v4(v_in, parameters);

dpsi_dy = discvar.dpsi_dy;
ddpsidy_dpsi = Ddiscvar.ddpsidy_dpsi;

dpsi_dz = discvar.dpsi_dz;
ddpsidz_dpsi = Ddiscvar.ddpsidz_dpsi;

ddomegady_domega = Ddiscvar.ddomegady_domega;
ddomegadz_domega = Ddiscvar.ddomegadz_domega;

u_vert = discvar.u_vert;
duvert_du = Ddiscvar.duvert_du;

ddudy_du =  Ddiscvar.ddudy_du;
du_dz = discvar.du_dz;
ddudz_du =  Ddiscvar.ddudz_du;
du_dz_centre = discvar.du_dz_centre;
d_dudzcentre_du = Ddiscvar.d_dudzcentre_du;
du_dy_centre = discvar.du_dy_centre;
ddudycentre_du = Ddiscvar.ddudycentre_du;

dphi_dy = discvar.dphi_dy;
ddphidy_dphi = Ddiscvar.ddphidy_dphi;
dphi_dz = discvar.dphi_dz;
ddphidz_dphi = Ddiscvar.ddphidz_dphi;

ddpdy_dp = Ddiscvar.ddpdy_dp;
ddpdz_dp = Ddiscvar.ddpdz_dp;

ddTdy_dT = Ddiscvar.ddTdy_dT;
ddTdz_dT = Ddiscvar.ddTdz_dT;

T_vert = discvar.T_vert;
dTvert_dT = Ddiscvar.dTver_dT;
T_hor = discvar.T_hor;
dThor_dT = Ddiscvar.dThor_dT;

ddTbdy_dTb = Ddiscvar.ddTbdy_dTb;
ddTbdz_dTb = Ddiscvar.ddTbdz_dTb;

%construct regularized friction coefficient
T_bed = T(T_bdy_nodes_bed);
u_bed = u(T_bdy_nodes_bed);

%T_bed, u_bed, Pi at psi cell centres
Tbed_upnode = circshift(1:length(T_bdy_nodes_bed),1,2);
Tbed_downnode = 1:length(T_bdy_nodes_bed);

T_bed_psigrid = (T_bed(Tbed_upnode)+T_bed(Tbed_downnode))./2;
dTbedpsigrid_dTbed_up = 1/2*ones(length(Tbed_upnode),1);
dTbedpsigrid_dTbed_down = 1/2*ones(length(Tbed_downnode),1);
index_dTbedpsigriddTbedup = T_bdy_nodes_bed(Tbed_upnode);
index_dTbedpsigriddTbeddown = T_bdy_nodes_bed(Tbed_downnode);

u_bed_psigrid = (u_bed(Tbed_upnode)+u_bed(Tbed_downnode))./2;
dubedpsigrid_dubed_up = dTbedpsigrid_dTbed_up;
dubedpsigrid_dubed_down = dTbedpsigrid_dTbed_down;
index_dubedpsigriddubedup = T_bdy_nodes_bed(Tbed_upnode);
index_dubedpsigriddubeddown = T_bdy_nodes_bed(Tbed_downnode);

Pi_psigrid = (Pi(Tbed_upnode)+Pi(Tbed_downnode))./2;
dPipsigrid_dPi_up = dTbedpsigrid_dTbed_up;
dPipsigrid_dPi_down = dTbedpsigrid_dTbed_down;
index_dPipsigriddPiup = Tbed_upnode;
index_dPipsigriddPidown = Tbed_downnode;

if isfield(parameters,'gamma_pert')==1
    gamma_pert = parameters.gamma_pert;
    gamma_pert_psigrid = (gamma_pert(Tbed_upnode)+gamma_pert(Tbed_downnode))./2;
elseif isfield(parameters,'gamma_pert')==0
    gamma_pert = zeros(length(T_bdy_nodes_bed),1);
    gamma_pert_psigrid = zeros(length(psi_bdy_nodes_bed),1);
end

[f_slide_hydro, Df_slide_hydro] = regularizedfriction_hydrology(u_bed,Pi, parameters);
Df_slide_hydro_dub = Df_slide_hydro.dub;
Df_slide_hydro_dPi = Df_slide_hydro.dPi;

[f_slide_hydro_psigrid, Df_slide_hydro_psigrid] = regularizedfriction_hydrology(u_bed_psigrid,Pi_psigrid, parameters);
Df_slide_hydro_psigrid_dub = Df_slide_hydro_psigrid.dub;
Df_slide_hydro_psigrid_dPi = Df_slide_hydro_psigrid.dPi;

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
dhav_dh = parameters.timestep.dx_prev./(parameters.timestep.dx_prev + parameters.timestep.dx)*ones(length(h),1);

%construct regularized bedwater content
if parameters.flag_Tdep == 1
    [f_slide_Tbed, DfslideTbed] = regularizedfriction_temperature(T_bed, parameters);
    [f_slide_Tbedpsigrid, DfslideTbedpsigrid] = regularizedfriction_temperature(T_bed_psigrid, parameters);
elseif parameters.flag_Tdep == 0
    f_slide_Tbed = ones(length(T_bdy_nodes_bed),1);
    DfslideTbed = zeros(length(T_bdy_nodes_bed),1);
    
    f_slide_Tbedpsigrid = ones(length(psi_bdy_nodes_bed),1);
    DfslideTbedpsigrid = zeros(length(psi_bdy_nodes_bed),1);
end

%% MASS CONSERVATION AND EQ. FOR h

%mass conservation
%fout(2*psi_nodes+4*T_nodes+Tb_nodes+length(psi_bdy_nodes_bed)+2) = Q - (a*dx_int + Q_prev);
fout(Q_index,Q_index) = 1./dx_int;

%conservation law
%fout(2*psi_nodes+4*T_nodes+Tb_nodes+length(psi_bdy_nodes_bed)+1) = Q - h_av.*u_int;
u_int = sum(u.*T_Delta_z_cell_volume.*T_Delta_y_cell);
duint_du = (T_Delta_z_cell_volume.*T_Delta_y_cell).';

fout(h_index, Q_index) = 1;
fout(h_index, u_index) = -h_av.*duint_du;
fout(h_index, h_index) = -dhav_dh.*u_int;

 %% STREAM FUNCTION
dpsihorflux_dpsi = psi_connect_hor*ddpsidy_dpsi;
dpsiverflux_dpsi = psi_connect_ver* (spdiags(h_prev^(-1).*psi_Delta_y_partial_cell,0,psi_nedges_ver, psi_nedges_ver)*ddpsidz_dpsi);

% Conditions on the tangential velocity, dpsi/dy = 0 is applied at top and bottom
% boundaries. We prescribe a Dirichlet condition at the top and allow for
% dpsi/dy = 0, otherwise psi is defined up to an unspecified constant

%top
%flux_psi_top = (8*psi_top -9*psi(psi_bdy_nodes_top) +...
    %psi(psi_bdy_nodes_top+psi_nodes_hor(1))).*psi_Delta_y_partial_cell(psi_bdy_nodes_top)./(3*psi_Delta_z_cell(psi_bdy_nodes_top).*h_prev);

dfluxpsitop_dpsi = sparse(psi_bdy_nodes_top, psi_bdy_nodes_top, -9*psi_Delta_y_partial_cell(psi_bdy_nodes_top)./(3*psi_Delta_z_cell(psi_bdy_nodes_top).*h_prev),psi_nodes,psi_nodes) + ...
    sparse(psi_bdy_nodes_top, psi_bdy_nodes_top+psi_nodes_hor(1), psi_Delta_y_partial_cell(psi_bdy_nodes_top)./(3*psi_Delta_z_cell(psi_bdy_nodes_top).*h_prev),psi_nodes,psi_nodes);
dpsiverflux_dpsi = dpsiverflux_dpsi + dfluxpsitop_dpsi;

%bottom 
% psi_bottom = psi_ghost;
% flux_psi_bottom = -(8*psi_bottom - 9*psi(psi_bdy_nodes_bed) + psi(psi_bdy_nodes_bed-psi_nodes_hor(end)))./(3*psi_Delta_z_cell(psi_bdy_nodes_bed).*h_prev);
% net_psi_verflux(psi_bdy_nodes_bed) = net_psi_verflux(psi_bdy_nodes_bed) - flux_psi_bottom.*psi_Delta_y_partial_cell(psi_bdy_nodes_bed);

psi_bottom = psi_ghost;
dpsibottom_dpsighost = ones(length(psi_ghost),1);

flux_psi_bottom = -(8*psi_bottom - 9*psi(psi_bdy_nodes_bed) + psi(psi_bdy_nodes_bed-psi_nodes_hor(end)))./(3*psi_Delta_z_cell(psi_bdy_nodes_bed).*h_prev);
dfluxpsibottom_dpsi = sparse(1:length(psi_bdy_nodes_bed),psi_bdy_nodes_bed, 9./(3*psi_Delta_z_cell(psi_bdy_nodes_bed).*h_prev), length(psi_bdy_nodes_bed), psi_nodes)+...
    sparse(1:length(psi_bdy_nodes_bed),psi_bdy_nodes_bed-psi_nodes_hor(end), -1./(3*psi_Delta_z_cell(psi_bdy_nodes_bed).*h_prev), length(psi_bdy_nodes_bed), psi_nodes);

dfluxpsibottom_dpsibed_vec = 9./(3*psi_Delta_z_cell(psi_bdy_nodes_bed).*h_prev);
dfluxpsibottom_dpsibed_vec_index = psi_bdy_nodes_bed;
dfluxpsibottom_dpsiabed_vec =-1./(3*psi_Delta_z_cell(psi_bdy_nodes_bed).*h_prev);
dfluxpsibottom_dpsiabed_vec_index = psi_bdy_nodes_bed-psi_nodes_hor(end);

dfluxpsibottom_dpsighost = -8*dpsibottom_dpsighost./(3*psi_Delta_z_cell(psi_bdy_nodes_bed).*h_prev);

%sum the pieces
%net_psi_verflux(psi_bdy_nodes_bed) = net_psi_verflux(psi_bdy_nodes_bed) - flux_psi_bottom.*psi_Delta_y_cell(psi_bdy_nodes_bed);
dpsiverflux_dpsi(psi_bdy_nodes_bed,:) = dpsiverflux_dpsi(psi_bdy_nodes_bed,:) - sparse(1:length(psi_bdy_nodes_bed),1:length(psi_bdy_nodes_bed),psi_Delta_y_cell(psi_bdy_nodes_bed),length(psi_bdy_nodes_bed),length(psi_bdy_nodes_bed))*dfluxpsibottom_dpsi;
dpsiverflux_dpsighost = - sparse(1:psi_nodes,1:psi_nodes,psi_Delta_y_cell,psi_nodes,psi_nodes)*sparse(psi_bdy_nodes_bed,1:length(psi_ghost),dfluxpsibottom_dpsighost, psi_nodes,length(psi_ghost));

%conservation law
% div_psifluxes =   net_psi_verflux./(h_prev.*psi_Delta_z_cell.*psi_Delta_y_cell)+ net_psi_horflux./psi_Delta_y_cell ;
% fout(1:psi_nodes) = div_psifluxes - omega;
fout(psi_index,psi_index) = spdiags(1./(h_prev.*psi_Delta_z_cell.*psi_Delta_y_cell),0,psi_nodes,psi_nodes)*dpsiverflux_dpsi + ...
    spdiags(1./(psi_Delta_y_cell),0,psi_nodes,psi_nodes)*dpsihorflux_dpsi;
fout(psi_index,psighost_index) = spdiags(1./(h_prev.*psi_Delta_z_cell.*psi_Delta_y_cell),0,psi_nodes,psi_nodes)*dpsiverflux_dpsighost;
fout(psi_index,omega_index) = -spdiags(ones(psi_nodes,1),0,psi_nodes,psi_nodes);

%% SLIDING LAW
%at T cell centres
u_bed = u(T_bdy_nodes_bed);
dubed_du = ones(length(T_bdy_nodes_bed),1);

%at psi cell centres
index_up_node_bed = circshift(1:length(T_bdy_nodes_bed),1,2);
index_down_node_bed = 1:length(T_bdy_nodes_bed);
dphi_dy_bed = (-phi(T_bdy_nodes_bed(index_up_node_bed))+phi(T_bdy_nodes_bed(index_down_node_bed)))./psi_Delta_y_cell(psi_bdy_nodes_bed);
ddphidybed_dphiup =  -1./psi_Delta_y_cell(psi_bdy_nodes_bed);
index_ddphidybed_dphiup = T_bdy_nodes_bed(index_up_node_bed);
ddphidybed_dphidown =  1./psi_Delta_y_cell(psi_bdy_nodes_bed);
index_ddphidybed_dphidown = T_bdy_nodes_bed(index_down_node_bed);

%omega_bed = (gamma.*(1+gamma_pert_psigrid).*f_slide_Tbedpsigrid.*f_slide_hydro_psigrid).*(dphi_dy_bed + flux_psi_bottom); 
domegabed_dfslideT = gamma.*(1+gamma_pert_psigrid).*f_slide_hydro_psigrid.*(dphi_dy_bed + flux_psi_bottom); %then still need to multiply by derivatives wrt T_bed
domegabed_ddphidybed =  (gamma.*(1+gamma_pert_psigrid).*f_slide_Tbedpsigrid.*f_slide_hydro_psigrid);
domegabed_dfluxpsibottom =  (gamma.*(1+gamma_pert_psigrid).*f_slide_Tbedpsigrid.*f_slide_hydro_psigrid);
domegabed_dfslidehydro = (gamma.*(1+gamma_pert_psigrid).*f_slide_Tbedpsigrid).*(dphi_dy_bed + flux_psi_bottom);

%% additional condition for psi_ghost, int_width omega_bed = 0
%fout(psighost_index) = sum(omega_bed.*psi_Delta_y_cell(psi_bdy_nodes_bed)); 

fout(psighost_index,T_index) = sum(sparse(1:length(psi_bdy_nodes_bed), index_dTbedpsigriddTbedup, domegabed_dfslideT.*DfslideTbedpsigrid.*dTbedpsigrid_dTbed_up.*psi_Delta_y_cell(psi_bdy_nodes_bed),length(psi_bdy_nodes_bed),T_nodes)+...
    sparse(1:length(psi_bdy_nodes_bed),index_dTbedpsigriddTbeddown, domegabed_dfslideT.*DfslideTbedpsigrid.*dTbedpsigrid_dTbed_down.*psi_Delta_y_cell(psi_bdy_nodes_bed),length(psi_bdy_nodes_bed),T_nodes), 1);
fout(psighost_index,u_index) = sum(sparse(1:length(psi_bdy_nodes_bed), index_dTbedpsigriddTbedup, domegabed_dfslidehydro.*Df_slide_hydro_psigrid_dub.*dubedpsigrid_dubed_up.*psi_Delta_y_cell(psi_bdy_nodes_bed),length(psi_bdy_nodes_bed),T_nodes)+...
    sparse(1:length(psi_bdy_nodes_bed),index_dTbedpsigriddTbeddown, domegabed_dfslidehydro.*Df_slide_hydro_psigrid_dub.*dubedpsigrid_dubed_down.*psi_Delta_y_cell(psi_bdy_nodes_bed),length(psi_bdy_nodes_bed),T_nodes), 1);
fout(psighost_index,Pi_index) = sum(sparse(1:length(psi_bdy_nodes_bed), index_dPipsigriddPiup, domegabed_dfslidehydro.*Df_slide_hydro_psigrid_dPi.*dPipsigrid_dPi_up.*psi_Delta_y_cell(psi_bdy_nodes_bed),length(psi_bdy_nodes_bed),length(Pi))+...
    sparse(1:length(psi_bdy_nodes_bed),index_dPipsigriddPidown, domegabed_dfslidehydro.*Df_slide_hydro_psigrid_dPi.*dPipsigrid_dPi_down.*psi_Delta_y_cell(psi_bdy_nodes_bed),length(psi_bdy_nodes_bed),length(Pi)), 1);
fout(psighost_index,phi_index) = sum(sparse(1:length(psi_bdy_nodes_bed), index_ddphidybed_dphiup, psi_Delta_y_cell(psi_bdy_nodes_bed).*domegabed_ddphidybed.*(ddphidybed_dphiup), length(psi_bdy_nodes_bed),T_nodes)+...
    sparse(1:length(psi_bdy_nodes_bed), index_ddphidybed_dphidown, psi_Delta_y_cell(psi_bdy_nodes_bed).*domegabed_ddphidybed.*(ddphidybed_dphidown), length(psi_bdy_nodes_bed),T_nodes),1) ;
fout(psighost_index,psi_index) = sum(sparse(1:length(psi_bdy_nodes_bed),1:length(psi_bdy_nodes_bed), domegabed_dfluxpsibottom.*psi_Delta_y_cell(psi_bdy_nodes_bed),length(psi_bdy_nodes_bed),length(psi_bdy_nodes_bed))*dfluxpsibottom_dpsi,1);
fout(psighost_index,psighost_index) = sum(sparse(1:length(psi_bdy_nodes_bed),1:length(psi_bdy_nodes_bed), domegabed_dfluxpsibottom.*psi_Delta_y_cell(psi_bdy_nodes_bed),length(psi_bdy_nodes_bed),length(psi_bdy_nodes_bed))*...
    sparse(1:length(psi_bdy_nodes_bed),1:length(psi_ghost),dfluxpsibottom_dpsighost, length(psi_bdy_nodes_bed),length(psi_ghost)),1);

%% VORTICITY
%net_omega_horflux = psi_connect_hor * domega_dy;
domegahorflux_domega = psi_connect_hor * ddomegady_domega;
%net_omega_verflux = psi_connect_ver*(h_prev.^(-1).*domega_dz.*psi_Delta_y_partial_cell); 
domegaverflux_domega = psi_connect_ver* (spdiags(h_prev^(-1).*psi_Delta_y_partial_cell,0,psi_nedges_ver, psi_nedges_ver)*ddomegadz_domega);

%enforce Dirichlet conditions on top and bottom boundaries (psi_x=0 on
%inflow and outflow bdies)

%ice surface: omega = - (h-h_prev)/dx*1/2*(du/dy_{prev} +du/dy_{current})
du_dy_top = (u(T_bdy_nodes_top(index_up_node_top_T))-u(T_bdy_nodes_top(index_down_node_top_T)))./psi_Delta_y_cell(psi_bdy_nodes_top);
du_dy_top_prev = (u_prev(T_bdy_nodes_top(index_up_node_top_T))-u_prev(T_bdy_nodes_top(index_down_node_top_T)))./psi_Delta_y_cell(psi_bdy_nodes_top);
ddudytop_du = sparse(psi_bdy_nodes_top,T_bdy_nodes_top(index_up_node_top_T), 1./psi_Delta_y_cell(psi_bdy_nodes_top),psi_nodes,T_nodes) + ...
    sparse(psi_bdy_nodes_top,T_bdy_nodes_top(index_down_node_top_T), -1./psi_Delta_y_cell(psi_bdy_nodes_top),psi_nodes,T_nodes);
ddudytop_duup_vec = 1./psi_Delta_y_cell(psi_bdy_nodes_top);
ddudytop_dudp_vec = -1./psi_Delta_y_cell(psi_bdy_nodes_top);

%omega_top = (h-h_prev)/(2*dx).*(du_dy_top_prev + du_dy_top);   
domegatop_du =  (h-h_prev)/(2*dx).*ddudytop_du;
domegatop_duu_vec = (h-h_prev)/(2*dx).*ddudytop_duup_vec;
domegatop_dud_vec = (h-h_prev)/(2*dx).*ddudytop_dudp_vec;
domegatop_dh = sparse(psi_bdy_nodes_top, ones(length(psi_bdy_nodes_top),1), 1/(2*dx).*(du_dy_top_prev + du_dy_top), psi_nodes, 1); 
domegatop_dh_vec = 1/(2*dx).*(du_dy_top_prev + du_dy_top);

%domegadztop = (-9*omega(psi_bdy_nodes_top) + omega(psi_bdy_nodes_top + length(psi_bdy_nodes_top)) +8*omega_top)./(3*psi_Delta_z_cell(psi_bdy_nodes_top));
ddomegadztop_domega = sparse(psi_bdy_nodes_top, psi_bdy_nodes_top, -9.*psi_Delta_y_cell(psi_bdy_nodes_top)./(3*psi_Delta_z_cell(psi_bdy_nodes_top)), psi_nodes,psi_nodes) + ...
    sparse(psi_bdy_nodes_top, psi_bdy_nodes_top + length(psi_bdy_nodes_top), 1.*psi_Delta_y_cell(psi_bdy_nodes_top)./(3*psi_Delta_z_cell(psi_bdy_nodes_top)), psi_nodes,psi_nodes);
ddomegadztop_du = 8*sparse(1:length(psi_bdy_nodes_top),1:length(psi_bdy_nodes_top), 1.*psi_Delta_y_cell(psi_bdy_nodes_top)./(3*psi_Delta_z_cell(psi_bdy_nodes_top)), psi_nodes,psi_nodes)*domegatop_du;
ddomegadztop_dh = 8*sparse(1:length(psi_bdy_nodes_top),1:length(psi_bdy_nodes_top), 1.*psi_Delta_y_cell(psi_bdy_nodes_top)./(3*psi_Delta_z_cell(psi_bdy_nodes_top)), psi_nodes,psi_nodes)*domegatop_dh;

%net_omega_verflux(psi_bdy_nodes_top) = net_omega_verflux(psi_bdy_nodes_top) + h_prev.^(-1).*domegadztop.*psi_Delta_y_cell(psi_bdy_nodes_top);
domegaverflux_domega = domegaverflux_domega + h_prev.^(-1).*ddomegadztop_domega;
domegaverflux_du =  h_prev.^(-1).*ddomegadztop_du;
domegaverflux_dh =  h_prev.^(-1).*ddomegadztop_dh;

%bed: omega = omega_bed -((8 a - 9 f1 + f2)/(3 z))
%domegadz_bed =  (9*omega(psi_bdy_nodes_bed) - omega(psi_bdy_nodes_bed - length(psi_bdy_nodes_bed))-8*omega_bed)./(3*psi_Delta_z_cell(psi_bdy_nodes_bed));
domegadzbed_domegabed = -8./(3*psi_Delta_z_cell(psi_bdy_nodes_bed));

ddomegadzbed_domega = sparse(psi_bdy_nodes_bed,psi_bdy_nodes_bed, 9.*psi_Delta_y_cell(psi_bdy_nodes_bed)./(3*psi_Delta_z_cell(psi_bdy_nodes_bed)),psi_nodes,psi_nodes)+ ...
    sparse(psi_bdy_nodes_bed,psi_bdy_nodes_bed - length(psi_bdy_nodes_bed), -1.*psi_Delta_y_cell(psi_bdy_nodes_bed)./(3*psi_Delta_z_cell(psi_bdy_nodes_bed)),psi_nodes,psi_nodes);
ddomegadzbed_dpsi = sparse(1:length(psi_bdy_nodes_bed),1:length(psi_bdy_nodes_bed), domegadzbed_domegabed.*psi_Delta_y_cell(psi_bdy_nodes_bed).*domegabed_dfluxpsibottom,length(psi_bdy_nodes_bed),length(psi_bdy_nodes_bed))*dfluxpsibottom_dpsi;
ddomegadzbed_dpsighost = sparse(psi_bdy_nodes_bed,1:length(psi_ghost), domegadzbed_domegabed.*psi_Delta_y_cell(psi_bdy_nodes_bed).*domegabed_dfluxpsibottom.*dfluxpsibottom_dpsighost,psi_nodes,length(psi_ghost));
ddomegadzbed_dphi = sparse(psi_bdy_nodes_bed,index_ddphidybed_dphiup, domegadzbed_domegabed.*psi_Delta_y_cell(psi_bdy_nodes_bed).*domegabed_ddphidybed.*ddphidybed_dphiup,psi_nodes,T_nodes) + ...
    sparse(psi_bdy_nodes_bed,index_ddphidybed_dphidown, domegadzbed_domegabed.*psi_Delta_y_cell(psi_bdy_nodes_bed).*domegabed_ddphidybed.*ddphidybed_dphidown,psi_nodes,T_nodes);
ddomegadzbed_dT = sparse(psi_bdy_nodes_bed,index_dTbedpsigriddTbedup, domegadzbed_domegabed.*psi_Delta_y_cell(psi_bdy_nodes_bed).*domegabed_dfslideT.*DfslideTbedpsigrid.*dTbedpsigrid_dTbed_up,psi_nodes,T_nodes)+...
    sparse(psi_bdy_nodes_bed,index_dTbedpsigriddTbeddown, domegadzbed_domegabed.*psi_Delta_y_cell(psi_bdy_nodes_bed).*domegabed_dfslideT.*DfslideTbedpsigrid.*dTbedpsigrid_dTbed_down,psi_nodes,T_nodes);
ddomegadzbed_du = sparse(psi_bdy_nodes_bed,index_dubedpsigriddubedup, domegadzbed_domegabed.*psi_Delta_y_cell(psi_bdy_nodes_bed).*domegabed_dfslidehydro.*Df_slide_hydro_psigrid_dub.*dubedpsigrid_dubed_up,psi_nodes,T_nodes)+...
    sparse(psi_bdy_nodes_bed,index_dubedpsigriddubeddown, domegadzbed_domegabed.*psi_Delta_y_cell(psi_bdy_nodes_bed).*domegabed_dfslidehydro.*Df_slide_hydro_psigrid_dub.*dubedpsigrid_dubed_down,psi_nodes,T_nodes);
ddomegadzbed_dPi = sparse(psi_bdy_nodes_bed,index_dPipsigriddPiup, domegadzbed_domegabed.*psi_Delta_y_cell(psi_bdy_nodes_bed).*domegabed_dfslidehydro.*Df_slide_hydro_psigrid_dPi.*dPipsigrid_dPi_up,psi_nodes,length(Pi))+...
    sparse(psi_bdy_nodes_bed,index_dPipsigriddPidown, domegadzbed_domegabed.*psi_Delta_y_cell(psi_bdy_nodes_bed).*domegabed_dfslidehydro.*Df_slide_hydro_psigrid_dPi.*dPipsigrid_dPi_down,psi_nodes,length(Pi));

%net_omega_verflux(psi_bdy_nodes_bed) = net_omega_verflux(psi_bdy_nodes_bed) - h_prev.^(-1).*psi_Delta_y_cell(psi_bdy_nodes_bed).*domegadz_bed;;
domegaverflux_domega = domegaverflux_domega - h_prev.^(-1).*ddomegadzbed_domega;
domegaverflux_dpsi = [sparse(psi_nodes-length(psi_bdy_nodes_bed), psi_nodes); - h_prev.^(-1).*ddomegadzbed_dpsi];
domegaverflux_dpsighost = - h_prev.^(-1).*ddomegadzbed_dpsighost;
domegaverflux_dT = - h_prev.^(-1).*ddomegadzbed_dT;
domegaverflux_dphi = - h_prev.^(-1).*ddomegadzbed_dphi;
domegaverflux_du = domegaverflux_du - h_prev.^(-1).*ddomegadzbed_du;
domegaverflux_dPi = - h_prev.^(-1).*ddomegadzbed_dPi;

%conservation law
%div_omegafluxes =   net_omega_verflux./(h_prev.*psi_Delta_z_cell.*psi_Delta_y_cell)+ net_omega_horflux./psi_Delta_y_cell ;
fout(omega_index,omega_index) = spdiags(1./(h_prev.*psi_Delta_z_cell.*psi_Delta_y_cell),0,psi_nodes,psi_nodes)*domegaverflux_domega + ...
   spdiags(1./(psi_Delta_y_cell),0,psi_nodes,psi_nodes)*domegahorflux_domega;
fout(omega_index,psi_index) = spdiags(1./(h_prev.*psi_Delta_z_cell.*psi_Delta_y_cell),0,psi_nodes,psi_nodes)*domegaverflux_dpsi;
fout(omega_index,psighost_index) = spdiags(1./(h_prev.*psi_Delta_z_cell.*psi_Delta_y_cell),0,psi_nodes,psi_nodes)*domegaverflux_dpsighost;
fout(omega_index,phi_index) = spdiags(1./(h_prev.*psi_Delta_z_cell.*psi_Delta_y_cell),0,psi_nodes,psi_nodes)*domegaverflux_dphi;
fout(omega_index,T_index) = sparse(1:psi_nodes, 1:psi_nodes,1./(h_prev.*psi_Delta_z_cell.*psi_Delta_y_cell),psi_nodes,psi_nodes)*domegaverflux_dT;
fout(omega_index,h_index) = spdiags(1./(h_prev.*psi_Delta_z_cell.*psi_Delta_y_cell),0,psi_nodes,psi_nodes)*domegaverflux_dh;
fout(omega_index,u_index) = spdiags(1./(h_prev.*psi_Delta_z_cell.*psi_Delta_y_cell),0,psi_nodes,psi_nodes)*domegaverflux_du;
fout(omega_index,Pi_index) = spdiags(1./(h_prev.*psi_Delta_z_cell.*psi_Delta_y_cell),0,psi_nodes,psi_nodes)*domegaverflux_dPi;

%% PHI
%diffusive fluxes
%net_phi_horflux = T_connect_hor * dphi_dy;
dphihorflux_dphi = T_connect_hor * ddphidy_dphi;

%net_phi_diff_verflux = T_connect_ver * (h_prev.^(-1).*dphi_dz.*T_Delta_y_partial_cell);;
dphiverfluxdiff_dphi = T_connect_ver * (spdiags(h_prev^(-1).*T_Delta_y_partial_cell,0,T_nedges_ver, T_nedges_ver)* ddphidz_dphi);

%advective flux
%flux_phi_adv = -1/2*T_coor_hor_cell_edges_z.*(u_vert.*(h - h_prev)/dx + (h_prev - h_pprev)/dx_prev*u_vert_prev);
dfluxphiadv_du = -1/2*(h - h_prev)/dx*sparse(1:length(u_vert), 1:length(u_vert),T_coor_hor_cell_edges_z.*T_Delta_y_partial_cell,length(u_vert),length(u_vert)) *duvert_du; 
dfluxphiadv_dh = -1/2*T_coor_hor_cell_edges_z.*u_vert.*T_Delta_y_partial_cell./dx;

%net_phi_adv_verflux = T_connect_ver * (flux_phi_adv.*T_Delta_y_partial_cell);
dnetfluxphiadv_du = T_connect_ver * dfluxphiadv_du;
dnetfluxphiadv_dh = T_connect_ver *dfluxphiadv_dh;

%net_phi_verflux =  net_phi_diff_verflux + net_phi_adv_verflux;
dnetphiverflux_dphi = dphiverfluxdiff_dphi;
dnetphiverflux_du = dnetfluxphiadv_du;
dnetphiverflux_dh = dnetfluxphiadv_dh;

%S_phi = (h_av*u - h_av_prev*u_prev)./dx_int;
dSphi_du = sparse(1:T_nodes,1:T_nodes, h_av/dx_int*ones(length(T_nodes),1), T_nodes,T_nodes);
dSphi_dh = sparse(1:T_nodes,ones(T_nodes,1), dhav_dh/dx_int*u, T_nodes,1);

%conservation law 
%div_phi_fluxes = net_phi_verflux./(h_prev*T_Delta_z_cell_volume.*T_Delta_y_cell)+ net_phi_horflux./T_Delta_y_cell;
%fout(2*psi_nodes+1:2*psi_nodes+T_nodes) = div_phi_fluxes + h_prev.^(-1)*S_phi;
fout(phi_index,phi_index) = spdiags(1./(h_prev.*T_Delta_z_cell_volume.*T_Delta_y_cell),0,T_nodes,T_nodes)*dnetphiverflux_dphi + ...
   spdiags(1./(T_Delta_y_cell),0,T_nodes,T_nodes)*dphihorflux_dphi;
fout(phi_index,u_index) = spdiags(1./(h_prev.*T_Delta_z_cell_volume.*T_Delta_y_cell),0,T_nodes,T_nodes)*dnetphiverflux_du + h_prev.^(-1)*dSphi_du;
fout(phi_index,h_index) = spdiags(1./(h_prev.*T_Delta_z_cell_volume.*T_Delta_y_cell),0,T_nodes,T_nodes)*dnetphiverflux_dh + h_prev.^(-1)*dSphi_dh;

fout(phi_index(1), :) = sparse(1, phi_index(1), 1, 1, length(v_in));

%% ALONG FLOW VELOCITY, U
%diffusive fluxes
%net_u_horflux = T_connect_hor * du_dy;
duhorflux_du = T_connect_hor * ddudy_du;

net_u_verflux = T_connect_ver * (h_av.^(-1).*du_dz.*T_Delta_y_partial_cell);
duverflux_du = T_connect_ver * (spdiags(h_av^(-1).*T_Delta_y_partial_cell,0,T_nedges_ver, T_nedges_ver)* ddudz_du);
duverflux_dh = T_connect_ver * (-h_av^(-2).*T_Delta_y_partial_cell.*du_dz .*dhav_dh);

%eq. 21: prescribed flux at cell centre at the bed

flux_u_bottom  = (gamma.*(1+gamma_pert).*f_slide_Tbed.*f_slide_hydro.*u_bed).*T_Delta_y_cell(T_bdy_nodes_bed); 
dfluxubottom_du = sparse(T_bdy_nodes_bed, T_bdy_nodes_bed,gamma.*(1+gamma_pert).*f_slide_Tbed.*f_slide_hydro.*T_Delta_y_cell(T_bdy_nodes_bed), T_nodes, T_nodes) + ...
    sparse(T_bdy_nodes_bed, T_bdy_nodes_bed,gamma.*(1+gamma_pert).*f_slide_Tbed.*Df_slide_hydro_dub.*u_bed.*T_Delta_y_cell(T_bdy_nodes_bed), T_nodes, T_nodes);
dfluxubottom_dT = sparse(T_bdy_nodes_bed, T_bdy_nodes_bed, gamma.*(1+gamma_pert).*f_slide_hydro.*u_bed.*DfslideTbed.*T_Delta_y_cell(T_bdy_nodes_bed),T_nodes,T_nodes);
dfluxubottom_dPi = sparse(T_bdy_nodes_bed, 1:length(T_bdy_nodes_bed),gamma.*(1+gamma_pert).*f_slide_Tbed.*Df_slide_hydro_dPi.*u_bed.*T_Delta_y_cell(T_bdy_nodes_bed), T_nodes, length(T_bdy_nodes_bed));

net_u_verflux(T_bdy_nodes_bed) = net_u_verflux(T_bdy_nodes_bed) - flux_u_bottom;
duverflux_du = duverflux_du - dfluxubottom_du;
duverflux_dT = -dfluxubottom_dT;
duverflux_dPi = - dfluxubottom_dPi;

%source term
%S_u = (h + b - h_prev - b_prev)/(dx).* ones(T_nodes, 1);
dSu_dh = 1./(dx).* ones(T_nodes, 1);

%conservation law
%fout(2*psi_nodes+T_nodes+1:2*psi_nodes+2*T_nodes) = div_u_fluxes - S_u;
%div_u_fluxes =   net_u_verflux./(h_av.*T_Delta_z_cell.*T_Delta_y_cell)+ net_u_horflux./T_Delta_y_cell; %nb: h_int depends on h, so must be differentiated!
fout(u_index,u_index) = spdiags(1./(h_av.*T_Delta_z_cell_volume.*T_Delta_y_cell),0,T_nodes,T_nodes)*duverflux_du + ...
   spdiags(1./(T_Delta_y_cell),0,T_nodes,T_nodes)*duhorflux_du;
fout(u_index,h_index) = spdiags(1./(h_av.*T_Delta_z_cell_volume.*T_Delta_y_cell),0,T_nodes,T_nodes)*duverflux_dh - net_u_verflux./(h_av.^2.*T_Delta_z_cell_volume.*T_Delta_y_cell)*dhav_dh - dSu_dh;
fout(u_index,T_index) = spdiags(1./(h_av.*T_Delta_z_cell_volume.*T_Delta_y_cell),0,T_nodes,T_nodes)*duverflux_dT;
fout(u_index,Pi_index) = spdiags(1./(h_av.*T_Delta_z_cell_volume.*T_Delta_y_cell),0,T_nodes,T_nodes)*duverflux_dPi;

%% P'
%diffusive fluxes
%diffusive fluxes
%net_p_horflux = T_connect_hor * dp_dy;
dphorflux_dp = T_connect_hor * ddpdy_dp;

%net_p_verflux = T_connect_ver * (h_prev^(-1).*dp_dz.*T_Delta_y_partial_cell);
dpverflux_dp = T_connect_ver * (spdiags(h_prev^(-1).*T_Delta_y_partial_cell,0,T_nedges_ver, T_nedges_ver)* ddpdz_dp);

%eq.32: prescribed flux at top nodes
% domegady_top = (omega_top(index_down_node_top_psi) - omega_top(index_up_node_top_psi))./T_Delta_y_cell(T_bdy_nodes_top);
domegadytop_domegatop = 1./T_Delta_y_cell(T_bdy_nodes_top);
domegatop_duu = domegatop_duu_vec;
domegatop_dud = domegatop_dud_vec;

indexdomegady_u = T_bdy_nodes_top(index_up_node_top_T);
indexdomegady_d = T_bdy_nodes_top(index_down_node_top_T);
indexcomposite_du = indexdomegady_u(index_down_node_top_psi);
indexcomposite_dd = indexdomegady_d(index_down_node_top_psi);
indexcomposite_uu = indexdomegady_u(index_up_node_top_psi);
indexcomposite_ud = indexdomegady_d(index_up_node_top_psi);

ddomegadytop_du = sparse(T_bdy_nodes_top,indexcomposite_du,domegadytop_domegatop.*domegatop_duu, T_nodes, T_nodes) +...
    sparse(T_bdy_nodes_top,indexcomposite_dd,domegadytop_domegatop.*domegatop_dud, T_nodes, T_nodes) - ...
    (sparse(T_bdy_nodes_top,indexcomposite_uu,domegadytop_domegatop.*domegatop_duu, T_nodes, T_nodes) +...
    sparse(T_bdy_nodes_top,indexcomposite_ud,domegadytop_domegatop.*domegatop_dud, T_nodes, T_nodes));

ddomegadytop_dh = sparse(T_bdy_nodes_top,ones(length(index_down_node_top_psi),1), domegatop_dh_vec(index_down_node_top_psi)./T_Delta_y_cell(T_bdy_nodes_top), T_nodes, 1) +...
    sparse(T_bdy_nodes_top,ones(length(index_up_node_top_psi),1),-domegatop_dh_vec(index_up_node_top_psi)./T_Delta_y_cell(T_bdy_nodes_top), T_nodes, 1);

% flux_p_top = -domegady_top.*T_Delta_y_cell(T_bdy_nodes_top);
dfluxptop_du = -sparse(1:T_nodes,1:T_nodes, T_Delta_y_cell, T_nodes, T_nodes)*ddomegadytop_du;
dfluxptop_dh = -sparse(1:T_nodes,1:T_nodes, T_Delta_y_cell, T_nodes, T_nodes)*ddomegadytop_dh;

%net_p_verflux(T_bdy_nodes_top) = net_p_verflux(T_bdy_nodes_top) + flux_p_top;
dpverflux_du = dfluxptop_du;
dpverflux_dh = dfluxptop_dh;

%eq. 32: prescribed flux at cell centre, bed
%domegady_bottom = (omega_bed(index_down_node_bed_psi) - omega_bed(index_up_node_bed_psi))./T_Delta_y_cell(T_bdy_nodes_bed);

ddomegadybottom_dpsi = sparse(T_bdy_nodes_bed,dfluxpsibottom_dpsibed_vec_index(index_down_node_bed_psi),domegabed_dfluxpsibottom(index_down_node_bed_psi).* dfluxpsibottom_dpsibed_vec(index_down_node_bed_psi)./T_Delta_y_cell(T_bdy_nodes_bed), T_nodes, psi_nodes) + ...
     sparse(T_bdy_nodes_bed,dfluxpsibottom_dpsiabed_vec_index(index_down_node_bed_psi),domegabed_dfluxpsibottom(index_down_node_bed_psi).* dfluxpsibottom_dpsiabed_vec(index_down_node_bed_psi)./T_Delta_y_cell(T_bdy_nodes_bed), T_nodes, psi_nodes) + ...
     sparse(T_bdy_nodes_bed,dfluxpsibottom_dpsibed_vec_index(index_up_node_bed_psi),-domegabed_dfluxpsibottom(index_up_node_bed_psi).* dfluxpsibottom_dpsibed_vec(index_up_node_bed_psi)./T_Delta_y_cell(T_bdy_nodes_bed), T_nodes, psi_nodes) + ...
     sparse(T_bdy_nodes_bed,dfluxpsibottom_dpsiabed_vec_index(index_up_node_bed_psi),-domegabed_dfluxpsibottom(index_up_node_bed_psi).* dfluxpsibottom_dpsiabed_vec(index_up_node_bed_psi)./T_Delta_y_cell(T_bdy_nodes_bed), T_nodes, psi_nodes);

ddomegadybottom_dpsighost = sparse(T_bdy_nodes_bed,ones(length(T_bdy_nodes_bed),1),domegabed_dfluxpsibottom(index_down_node_bed_psi).*dfluxpsibottom_dpsighost(index_down_node_bed_psi)./T_Delta_y_cell(T_bdy_nodes_bed), T_nodes, 1) + ...
    sparse(T_bdy_nodes_bed,ones(length(T_bdy_nodes_bed),1),-domegabed_dfluxpsibottom(index_up_node_bed_psi).*dfluxpsibottom_dpsighost(index_up_node_bed_psi)./T_Delta_y_cell(T_bdy_nodes_bed), T_nodes, 1);

ddomegadybottom_ddphidybed_up = domegabed_ddphidybed./T_Delta_y_cell(T_bdy_nodes_bed).*ddphidybed_dphiup;
ddomegadybottom_ddphidybed_down = domegabed_ddphidybed./T_Delta_y_cell(T_bdy_nodes_bed).*ddphidybed_dphidown;

indexcomposite_ddphidybed_du = index_ddphidybed_dphiup(index_down_node_bed_psi);
indexcomposite_ddphidybed_dd = index_ddphidybed_dphidown(index_down_node_bed_psi);
indexcomposite_ddphidybed_uu = index_ddphidybed_dphiup(index_up_node_bed_psi);
indexcomposite_ddphidybed_ud = index_ddphidybed_dphidown(index_up_node_bed_psi);

ddomegadybottom_dphi = sparse(T_bdy_nodes_bed,indexcomposite_ddphidybed_du, ddomegadybottom_ddphidybed_up(index_down_node_bed_psi), T_nodes, T_nodes) +...
    sparse(T_bdy_nodes_bed,indexcomposite_ddphidybed_dd, ddomegadybottom_ddphidybed_down(index_down_node_bed_psi), T_nodes, T_nodes) - ...
    (sparse(T_bdy_nodes_bed,indexcomposite_ddphidybed_uu, ddomegadybottom_ddphidybed_up(index_up_node_bed_psi), T_nodes, T_nodes) +...
    sparse(T_bdy_nodes_bed,indexcomposite_ddphidybed_ud, ddomegadybottom_ddphidybed_down(index_up_node_bed_psi), T_nodes, T_nodes));

ddomegadybottom_dT = sparse(T_bdy_nodes_bed,index_dTbedpsigriddTbedup(index_down_node_bed_psi),domegabed_dfslideT(index_down_node_bed_psi).*DfslideTbedpsigrid(index_down_node_bed_psi).*dTbedpsigrid_dTbed_up(index_down_node_bed_psi)./T_Delta_y_cell(T_bdy_nodes_bed), T_nodes, T_nodes) + ...
    sparse(T_bdy_nodes_bed,index_dTbedpsigriddTbeddown(index_down_node_bed_psi),domegabed_dfslideT(index_down_node_bed_psi).*DfslideTbedpsigrid(index_down_node_bed_psi).*dTbedpsigrid_dTbed_down(index_down_node_bed_psi)./T_Delta_y_cell(T_bdy_nodes_bed), T_nodes, T_nodes) - ...
    (sparse(T_bdy_nodes_bed,index_dTbedpsigriddTbedup(index_up_node_bed_psi), domegabed_dfslideT(index_up_node_bed_psi).*DfslideTbedpsigrid(index_up_node_bed_psi).*dTbedpsigrid_dTbed_up(index_up_node_bed_psi)./T_Delta_y_cell(T_bdy_nodes_bed), T_nodes, T_nodes) + ...
    sparse(T_bdy_nodes_bed,index_dTbedpsigriddTbeddown(index_up_node_bed_psi),domegabed_dfslideT(index_up_node_bed_psi).*DfslideTbedpsigrid(index_up_node_bed_psi).*dTbedpsigrid_dTbed_down(index_up_node_bed_psi)./T_Delta_y_cell(T_bdy_nodes_bed), T_nodes, T_nodes));

domegaduup = domegabed_dfslidehydro.*Df_slide_hydro_psigrid_dub.*dubedpsigrid_dubed_up;
domegadudown = domegabed_dfslidehydro.*Df_slide_hydro_psigrid_dub.*dubedpsigrid_dubed_down;
ddomegadybottom_du = sparse(T_bdy_nodes_bed,index_dubedpsigriddubedup(index_down_node_bed_psi),domegaduup(index_down_node_bed_psi)./T_Delta_y_cell(T_bdy_nodes_bed), T_nodes, T_nodes) + ...
    sparse(T_bdy_nodes_bed,index_dubedpsigriddubeddown(index_down_node_bed_psi),domegadudown(index_down_node_bed_psi)./T_Delta_y_cell(T_bdy_nodes_bed), T_nodes, T_nodes) - ...
    (sparse(T_bdy_nodes_bed,index_dubedpsigriddubedup(index_up_node_bed_psi), domegaduup(index_up_node_bed_psi)./T_Delta_y_cell(T_bdy_nodes_bed), T_nodes, T_nodes) + ...
    sparse(T_bdy_nodes_bed,index_dubedpsigriddubeddown(index_up_node_bed_psi),domegadudown(index_up_node_bed_psi)./T_Delta_y_cell(T_bdy_nodes_bed), T_nodes, T_nodes));

domegadPiup = domegabed_dfslidehydro.*Df_slide_hydro_psigrid_dPi.*dPipsigrid_dPi_up;
domegadPidown = domegabed_dfslidehydro.*Df_slide_hydro_psigrid_dPi.*dPipsigrid_dPi_down;
ddomegadybottom_dPi = sparse(T_bdy_nodes_bed,index_dPipsigriddPiup(index_down_node_bed_psi),domegadPiup(index_down_node_bed_psi)./T_Delta_y_cell(T_bdy_nodes_bed), T_nodes, length(T_bdy_nodes_bed)) + ...
    sparse(T_bdy_nodes_bed,index_dPipsigriddPidown(index_down_node_bed_psi),domegadPidown(index_down_node_bed_psi)./T_Delta_y_cell(T_bdy_nodes_bed), T_nodes, length(T_bdy_nodes_bed)) - ...
    (sparse(T_bdy_nodes_bed,index_dPipsigriddPiup(index_up_node_bed_psi), domegadPiup(index_up_node_bed_psi)./T_Delta_y_cell(T_bdy_nodes_bed), T_nodes, length(T_bdy_nodes_bed)) + ...
    sparse(T_bdy_nodes_bed,index_dPipsigriddPidown(index_up_node_bed_psi),domegadPidown(index_up_node_bed_psi)./T_Delta_y_cell(T_bdy_nodes_bed), T_nodes, length(T_bdy_nodes_bed)));

%flux_p_bottom  = -domegady_bottom.*T_Delta_y_cell(T_bdy_nodes_bed); 
dfluxpbottom_dphi = -sparse(1:T_nodes,1:T_nodes,T_Delta_y_cell,T_nodes,T_nodes)*ddomegadybottom_dphi;
dfluxpbottom_dpsi = -sparse(1:T_nodes,1:T_nodes,T_Delta_y_cell,T_nodes,T_nodes)*ddomegadybottom_dpsi;
dfluxpbottom_dT = -sparse(1:T_nodes,1:T_nodes,T_Delta_y_cell,T_nodes,T_nodes)*ddomegadybottom_dT;
dfluxpbottom_dpsighost = -sparse(1:T_nodes,1:T_nodes,T_Delta_y_cell,T_nodes,T_nodes)*ddomegadybottom_dpsighost;
dfluxpbottom_du = -sparse(1:T_nodes,1:T_nodes,T_Delta_y_cell,T_nodes,T_nodes)*ddomegadybottom_du;
dfluxpbottom_dPi = -sparse(1:T_nodes,1:T_nodes,T_Delta_y_cell,T_nodes,T_nodes)*ddomegadybottom_dPi;

%net_p_verflux(T_bdy_nodes_bed) = net_p_verflux(T_bdy_nodes_bed) - flux_p_bottom;
dpverflux_dphi = -dfluxpbottom_dphi;
dpverflux_dpsi = -dfluxpbottom_dpsi;
dpverflux_dT = -dfluxpbottom_dT;
dpverflux_dpsighost = -dfluxpbottom_dpsighost;
dpverflux_du = dpverflux_du -dfluxpbottom_du;
dpverflux_dPi = -dfluxpbottom_dPi;

%dont forget dpverflux_du = dfluxptop_du; dpverflux_dh = dfluxptop_dh;
%conservation law 
%div_p_fluxes =   net_p_verflux./(h_prev.*T_Delta_z_cell.*T_Delta_y_cell)+ net_p_horflux./T_Delta_y_cell;
%fout(2*psi_nodes+2*T_nodes+1: 2*psi_nodes+3*T_nodes) = div_p_fluxes;
fout(p_index,p_index) = sparse(1:T_nodes, 1:T_nodes,1./(h_prev.*T_Delta_z_cell_volume.*T_Delta_y_cell),T_nodes,T_nodes)*dpverflux_dp+ sparse(1:T_nodes, 1:T_nodes,1./(T_Delta_y_cell),T_nodes,T_nodes)*dphorflux_dp;
fout(p_index,u_index) = sparse(1:T_nodes, 1:T_nodes,1./(h_prev.*T_Delta_z_cell_volume.*T_Delta_y_cell),T_nodes,T_nodes)*dpverflux_du;
fout(p_index,h_index) = sparse(1:T_nodes, 1:T_nodes,1./(h_prev.*T_Delta_z_cell_volume.*T_Delta_y_cell),T_nodes,T_nodes)*dpverflux_dh;
fout(p_index,phi_index) = sparse(1:T_nodes, 1:T_nodes,1./(h_prev.*T_Delta_z_cell_volume.*T_Delta_y_cell),T_nodes,T_nodes)*dpverflux_dphi;
fout(p_index,psi_index) = sparse(1:T_nodes, 1:T_nodes,1./(h_prev.*T_Delta_z_cell_volume.*T_Delta_y_cell),T_nodes,T_nodes)*dpverflux_dpsi;
fout(p_index,T_index) = sparse(1:T_nodes, 1:T_nodes,1./(h_prev.*T_Delta_z_cell_volume.*T_Delta_y_cell),T_nodes,T_nodes)*dpverflux_dT;
fout(p_index,psighost_index) = sparse(1:T_nodes, 1:T_nodes,1./(h_prev.*T_Delta_z_cell_volume.*T_Delta_y_cell),T_nodes,T_nodes)*dpverflux_dpsighost;
fout(p_index,Pi_index) = sparse(1:T_nodes, 1:T_nodes,1./(h_prev.*T_Delta_z_cell_volume.*T_Delta_y_cell),T_nodes,T_nodes)*dpverflux_dPi;

%NOTE: prescribe dirichlet condition on *one* surface node so as to fix p
%fout(2*psi_nodes+2*T_nodes+1) = p(1);

fout(p_index(1),:) = sparse(1,p_index(1),1,1,length(v_in));


%% CONSTITUTIVE RELATIONS for p(p'), v, w_eff

%notation is p -> p', pr->p
pr = p - (u-u_prev)./dx_int;  
dpr_dp = sparse(1:T_nodes, 1:T_nodes, ones(length(p),1),T_nodes,T_nodes);
dpr_du = sparse(1:T_nodes, 1:T_nodes, -ones(length(u),1)./dx_int,T_nodes,T_nodes);

%% transverse velocity, eq. 34
dpsi_dz_full = [sparse(length(psi_bdy_nodes_top),1); dpsi_dz; sparse(length(psi_bdy_nodes_bed),1)];
ddpsidzfull_dpsi = [sparse(length(psi_bdy_nodes_top), psi_nodes); ddpsidz_dpsi; sparse(length(psi_bdy_nodes_bed), psi_nodes)];

v = dphi_dy + h_prev^(-1).*dpsi_dz_full; 
dv_dpsi = h_prev^(-1).*ddpsidzfull_dpsi;
dv_dphi = ddphidy_dphi;

%% vertical velocity, w_eff, on T grid
[dpsi_dy_reindexed, ddpsidyreindexed_dpsi] = verticalvelocity(parameters,dpsi_dy,Ddiscvar.ddpsidy_dpsi);
w_eff = -dpsi_dy_reindexed + h_prev.^(-1).*dphi_dz - T_coor_hor_cell_edges_z./2.*(u_vert *(h-h_prev)/dx + u_vert_prev .*(h_prev-h_pprev)/dx_prev);
dweff_dh = - T_coor_hor_cell_edges_z./2.*(u_vert/dx);
dweff_ddphidz = h_prev.^(-1);
dweff_duvert = - T_coor_hor_cell_edges_z./2.*(h-h_prev)/dx;

if parameters.flag1d == 0
    dweff_ddpsidy = -1;
elseif parameters.flag1d == 1
    dweff_ddpsidy = 0;
end
%% HEAT EQUATION, ICE

%diffusive fluxes
%net_T_diffhorflux = T_connect_hor * dT_dy;
dTdiffhorflux_dT = T_connect_hor * ddTdy_dT;

%net_T_diffverflux = T_connect_ver*( dT_dz.*T_Delta_y_partial_cell);
dTdiffverflux_dT =  T_connect_ver * (spdiags(T_Delta_y_partial_cell,0,T_nedges_ver, T_nedges_ver)* ddTdz_dT);

%advective fluxes 
%vertical
%T_adv_flux_ver = w_eff.*T_vert;
dTadverflux_ddphidz = dweff_ddphidz.*T_vert;
dTadverflux_du = spdiags(T_Delta_y_partial_cell.*dweff_duvert.*T_vert,0, length(T_vert), length(T_vert))*duvert_du;
dTadverflux_dpsi = spdiags(T_Delta_y_partial_cell.*dweff_ddpsidy.*T_vert,0, length(T_vert), length(T_vert))*ddpsidyreindexed_dpsi;
dTadverflux_dTver =  w_eff;
dTadverflux_dh =  dweff_dh.*T_vert;

%net_T_advfluxver = T_connect_ver*(T_adv_flux_ver.*T_Delta_y_partial_cell);
dnetTadvverflux_dT = T_connect_ver * ( spdiags(T_Delta_y_partial_cell.*dTadverflux_dTver ,0,T_nedges_ver, T_nedges_ver)*dTvert_dT);
dnetTadvverflux_dh =  T_connect_ver * (T_Delta_y_partial_cell.*dTadverflux_dh);
dnetTadvverflux_dphi = T_connect_ver * ( spdiags(T_Delta_y_partial_cell.*dTadverflux_ddphidz ,0,T_nedges_ver, T_nedges_ver)*ddphidz_dphi);
dnetTadvverflux_du = T_connect_ver * dTadverflux_du;
dnetTadvverflux_dpsi = T_connect_ver * dTadverflux_dpsi;

%transverse
%T_adv_flux_transverse = v.*T_hor.*h_prev;
dTadvfluxtransvers_dT = spdiags(v.*h_prev,0,T_nodes, T_nodes)*dThor_dT;
dTadvfluxtransvers_dpsi = spdiags(h_prev.*T_hor,0,T_nodes, T_nodes)*dv_dpsi;
dTadvfluxtransvers_dphi = spdiags(h_prev.*T_hor,0,T_nodes, T_nodes)*dv_dphi;

%net_T_advfluxtransverse = T_connect_hor *T_adv_flux_transverse;
dnetTadvtransflux_dpsi = T_connect_hor *dTadvfluxtransvers_dpsi;
dnetTadvtransflux_dphi = T_connect_hor *dTadvfluxtransvers_dphi;
dnetTadvtransflux_dT = T_connect_hor *dTadvfluxtransvers_dT;

%downstream

if parameters.flag.plug == 0
    %net_T_advfluxalong = (u.*h_prev.*T -u_prev.*h_pprev.*T_prev)./dx_int;
    dnetTadvalongflux_dT = sparse(1:T_nodes, 1:T_nodes,u.*h_prev./dx_int,T_nodes,T_nodes);
    dnetTadvalongflux_du = sparse(1:T_nodes, 1:T_nodes,T.*h_prev./dx_int,T_nodes,T_nodes);
elseif parameters.flag.plug == 1
    %net_T_advfluxalong = (u_prev_T.*h_prev_T.*T -u_prev_T.*h_prev_T.*T_prev)./dx_int; 
    dnetTadvalongflux_dT = sparse(1:T_nodes, 1:T_nodes,u_prev_T.*h_prev_T./dx_int,T_nodes,T_nodes);
    dnetTadvalongflux_du = sparse(T_nodes,T_nodes);
    
end

%source term
%S_T = (2*h_prev)^(-2).*((du_dz_centre_full).^2 + (du_dz_centre_full_prev).^2 + 2.*du_dz_centre_full.*du_dz_centre_full_prev)+...
    %1/4*((du_dy_centre).^2 + (du_dy_centre_prev).^2 + 2.*du_dy_centre.*du_dy_centre_prev);
du_dz_centre_top = sparse(length(T_bdy_nodes_top),1);
ddudzcentretop_du = sparse(length(T_bdy_nodes_top),T_nodes);
du_dz_centre_bottom = h*gamma.*(1+gamma_pert).*f_slide_Tbed.*f_slide_hydro.*u_bed; 
ddudzcentrebottom_du = sparse(1:length(T_bdy_nodes_bed), T_bdy_nodes_bed, h*gamma.*(1+gamma_pert).*f_slide_Tbed.*f_slide_hydro.*dubed_du, length(T_bdy_nodes_bed), T_nodes) + ...
    sparse(1:length(T_bdy_nodes_bed), T_bdy_nodes_bed, h*gamma.*(1+gamma_pert).*f_slide_Tbed.*u_bed.*Df_slide_hydro_dub.*dubed_du, length(T_bdy_nodes_bed), T_nodes); 
ddudzcentrebottom_dh = sparse(T_bdy_nodes_bed, ones(length(T_bdy_nodes_bed),1), gamma.*(1+gamma_pert).*f_slide_Tbed.*f_slide_hydro.*u_bed, T_nodes, 1);
ddudzcentrebottom_dT = sparse(T_bdy_nodes_bed, T_bdy_nodes_bed, h*gamma.*(1+gamma_pert).*f_slide_hydro.*u_bed.*DfslideTbed, T_nodes,T_nodes);
ddudzcentrebottom_dPi = sparse(T_bdy_nodes_bed, 1:length(T_bdy_nodes_bed), h*gamma.*(1+gamma_pert).*f_slide_Tbed.*u_bed.*Df_slide_hydro_dPi, T_nodes, length(T_bdy_nodes_bed)); 

du_dz_centre_full = [du_dz_centre_top; du_dz_centre; du_dz_centre_bottom];
ddudzcentrefull_du = [ddudzcentretop_du; d_dudzcentre_du; ddudzcentrebottom_du];
ddudzcentrefull_dh = ddudzcentrebottom_dh;
ddudzcentrefull_dT = ddudzcentrebottom_dT;
ddudzcentrefull_dPi = ddudzcentrebottom_dPi;

dST_du = (sparse(1:T_nodes, 1:T_nodes,(2*h_prev)^(-2).*2*du_dz_centre_full,T_nodes,T_nodes) + sparse(1:T_nodes, 1:T_nodes,(2*h_prev)^(-2).*2*du_dz_centre_full_prev,T_nodes,T_nodes))*ddudzcentrefull_du + ...
    1/4*(sparse(1:T_nodes, 1:T_nodes,2*du_dy_centre,T_nodes,T_nodes) + sparse(1:T_nodes, 1:T_nodes,2*du_dy_centre_prev,T_nodes,T_nodes))*ddudycentre_du;
dST_dh = (sparse(1:T_nodes,1:T_nodes,(2*h_prev)^(-2).*2*du_dz_centre_full,T_nodes,T_nodes) + sparse(1:T_nodes,1:T_nodes,(2*h_prev)^(-2).*2*du_dz_centre_full_prev,T_nodes,T_nodes))*ddudzcentrefull_dh;
%dST_dT = (spdiags((2*h_prev)^(-2).*2*du_dz_centre_full,0,T_nodes) + spdiags((2*h_prev)^(-2).*2*du_dz_centre_full_prev,0,T_nodes))*ddudzcentrefull_dTbed;

%conservation law
% fout(2*psi_nodes+3*T_nodes+1: 2*psi_nodes+4*T_nodes) = ...
% Pe/h_prev.*(net_T_advfluxver./(T_Delta_z_cell.*T_Delta_y_cell)+ net_T_advfluxtransverse./T_Delta_y_cell + net_T_advfluxalong)+...
% net_T_diffverflux./(h_prev.^2*T_Delta_z_cell.*T_Delta_y_cell)+ net_T_diffhorflux./T_Delta_y_cell - alpha*S_T;

if parameters.flag.heat_full == 1 || isfield(parameters,'flag') == 0 || parameters.flag.plug == 0
    
    fout(T_index,T_index) = Pe/h_prev*(sparse(1:T_nodes,1:T_nodes,1./(T_Delta_z_cell.*T_Delta_y_cell), T_nodes,T_nodes)*dnetTadvverflux_dT + ...
        sparse(1:T_nodes,1:T_nodes,1./(T_Delta_y_cell),T_nodes,T_nodes)*dnetTadvtransflux_dT + dnetTadvalongflux_dT ) - ...
        sparse(1:T_nodes,1:T_nodes,1./(h_prev.^2*T_Delta_z_cell.*T_Delta_y_cell), T_nodes, T_nodes)*dTdiffverflux_dT -  sparse(1:T_nodes, 1:T_nodes,1./T_Delta_y_cell, T_nodes, T_nodes)*dTdiffhorflux_dT;% - alpha*dST_dT;
    fout(T_index,u_index) = Pe/h_prev*(sparse(1:T_nodes, 1:T_nodes, 1./(T_Delta_z_cell.*T_Delta_y_cell), T_nodes,T_nodes)*dnetTadvverflux_du + dnetTadvalongflux_du )- alpha*dST_du;
    fout(T_index,phi_index) = Pe/h_prev*(sparse(1:T_nodes,1:T_nodes,1./(T_Delta_z_cell.*T_Delta_y_cell), T_nodes,T_nodes)*dnetTadvverflux_dphi + sparse(1:T_nodes, 1:T_nodes, 1./(T_Delta_y_cell), T_nodes,T_nodes)*dnetTadvtransflux_dphi);
    fout(T_index,psi_index) = Pe/h_prev*( sparse(1:T_nodes,1:T_nodes,1./(T_Delta_z_cell.*T_Delta_y_cell), T_nodes,T_nodes)*dnetTadvverflux_dpsi +sparse(1:T_nodes, 1:T_nodes, 1./(T_Delta_y_cell), T_nodes,T_nodes)*dnetTadvtransflux_dpsi);
    fout(T_index,h_index) = Pe/h_prev*(sparse(1:T_nodes,1:T_nodes,1./(T_Delta_z_cell.*T_Delta_y_cell), T_nodes,T_nodes)*dnetTadvverflux_dh)- alpha*dST_dh;
elseif parameters.flag.heat_full == 0 || parameters.flag.plug == 1
    fout(T_index,T_index) = Pe/h_prev*(dnetTadvalongflux_dT ) - ...
        sparse(1:T_nodes,1:T_nodes,1./(h_prev.^2*T_Delta_z_cell.*T_Delta_y_cell), T_nodes, T_nodes)*dTdiffverflux_dT -  sparse(1:T_nodes, 1:T_nodes,1./T_Delta_y_cell, T_nodes, T_nodes)*dTdiffhorflux_dT;% - alpha*dST_dT;
end

%surface and bed boundary condition
%fout(2*psi_nodes+3*T_nodes+T_bdy_nodes_top) = T(T_bdy_nodes_top) - T_s;
fout(T_index(T_bdy_nodes_top),:) = sparse(length(T_bdy_nodes_top),size(fout,2));
fout(T_index(T_bdy_nodes_top),T_index) = sparse(1:length(T_bdy_nodes_top),T_bdy_nodes_top,ones(length(T_bdy_nodes_top),1), length(T_bdy_nodes_top),T_nodes);

% %fout(2*psi_nodes+3*T_nodes+T_bdy_nodes_bed) = T(T_bdy_nodes_bed) - T_bed;
fout(T_index(T_bdy_nodes_bed),:) = sparse(length(T_bdy_nodes_bed),length(v_in));

%% HEAT EQUATION, BED
%net_Tb_diffhorflux = Tb_connect_hor * (-dTb_dy);
dTbdiffhorflux_dTb = Tb_connect_hor * (-ddTbdy_dTb);

%net_Tb_diffverflux = Tb_connect_ver *( -dTb_dz.*Tb_Delta_y_partial_cell);
dTbdiffverflux_dTb = Tb_connect_ver * (-spdiags(Tb_Delta_y_partial_cell,0,Tb_nedges_ver, Tb_nedges_ver)*ddTbdz_dTb);

%neumann conditions at top and bottom
%bedflux_below = - (T(T_bdy_nodes_bed) - Tb(Tb_bdy_nodes_top)).*Tb_Delta_y_cell(Tb_bdy_nodes_top)./(Tb_Delta_z_cell(Tb_bdy_nodes_top));   
dbedfluxbelow_dTb = Tb_Delta_y_cell(Tb_bdy_nodes_top)./(Tb_Delta_z_cell(Tb_bdy_nodes_top));
dbedfluxbelow_dT = -Tb_Delta_y_cell(Tb_bdy_nodes_top)./(Tb_Delta_z_cell(Tb_bdy_nodes_top));   

%net_Tb_diffverflux(Tb_bdy_nodes_top) = net_Tb_diffverflux(Tb_bdy_nodes_top) + bedflux_below;
dTbdiffverflux_dTb(Tb_bdy_nodes_top,:) = dTbdiffverflux_dTb(Tb_bdy_nodes_top,:) + sparse(1:length(Tb_bdy_nodes_top), Tb_bdy_nodes_top, dbedfluxbelow_dTb,length(Tb_bdy_nodes_top),Tb_nodes);
dTbdiffverflux_dT = sparse(1:length(Tb_bdy_nodes_top), T_bdy_nodes_bed, dbedfluxbelow_dT,Tb_nodes,T_nodes);

%conservation law
%fout(2*psi_nodes+4*T_nodes+1: 2*psi_nodes+4*T_nodes+Tb_nodes) = net_Tb_diffverflux./(Tb_Delta_z_cell.*Tb_Delta_y_cell)+ net_Tb_diffhorflux./Tb_Delta_y_cell;
fout(Tb_index,Tb_index) = sparse(1:Tb_nodes, 1:Tb_nodes,1./(Tb_Delta_z_cell.*Tb_Delta_y_cell),Tb_nodes,Tb_nodes)*dTbdiffverflux_dTb +...
    sparse(1:Tb_nodes, 1:Tb_nodes,1./(Tb_Delta_y_cell),Tb_nodes,Tb_nodes)*dTbdiffhorflux_dTb;
fout(Tb_index,T_index) = sparse(1:Tb_nodes, 1:Tb_nodes,1./(Tb_Delta_z_cell.*Tb_Delta_y_cell),Tb_nodes,Tb_nodes)*dTbdiffverflux_dT;

%% DRAINAGE
%defined modified effective pressure (this is an auxiliary variable, living
%at T_nodes)

%enforce frozen/temperate bed with indicator function from previous
%timestep, eq. 42-43
index_temperate_prev = find(I_prev == 1); 
index_frozen_prev = find(I_prev == 0); 

%nb: the two equations below all together take n_rows =
%length(T_bdy_nodes_bed);
% fout(Pi_index(index_frozen_prev)) = Pi(index_frozen_prev);
fout(Pi_index(index_frozen_prev), Pi_index(index_frozen_prev)) = sparse(1:length(index_frozen_prev), 1:length(index_frozen_prev), ones(length(index_frozen_prev),1),length(index_frozen_prev),length(index_frozen_prev));
% fout(T_index(T_bdy_nodes_bed(index_temperate_prev))) = T(T_bdy_nodes_bed(index_temperate_prev));
fout(T_index(T_bdy_nodes_bed(index_temperate_prev)), T_index(T_bdy_nodes_bed(index_temperate_prev))) = sparse(1:length(T_bdy_nodes_bed(index_temperate_prev)), 1:length(T_bdy_nodes_bed(index_temperate_prev)), ones(length(index_temperate_prev),1),length(index_temperate_prev),length(index_temperate_prev));

%bed state indicator at current time step, eq. 47
I = sparse(length(T_bdy_nodes_bed),1);
index_temperate_T = find(T_bed>0);
index_incipient = find(T_bed==0);
index_temperate_P = index_incipient(Pi(index_incipient)>0);

I(index_temperate_T) = ones(length(index_temperate_T),1);
I(index_temperate_P) = ones(length(index_temperate_P),1);

%normal stress at the base, eq. 44
%second order accurate, one-sided interpolation formula for second
%derivative of phi at the base

if isfield(parameters.flag, 'sigma') == 0 ||  strcmp(parameters.flag.sigma, 'FV') == 1
    ddphi_dz_bed = 2./(2*(T_Delta_z_cell(T_bdy_nodes_bed)).^2).*(-2*phi(T_bdy_nodes_bed-T_nodes_hor(end))+ phi(T_bdy_nodes_bed-2*T_nodes_hor(end)) +phi(T_bdy_nodes_bed));
    
    dddphidzbed_dphi = sparse(1:length(T_bdy_nodes_bed),T_bdy_nodes_bed-T_nodes_hor(end),2./(2*(T_Delta_z_cell(T_bdy_nodes_bed)).^2).*(-2), length(T_bdy_nodes_bed),T_nodes)+...
        sparse(1:length(T_bdy_nodes_bed),T_bdy_nodes_bed-2*T_nodes_hor(end),2./(2*(T_Delta_z_cell(T_bdy_nodes_bed)).^2), length(T_bdy_nodes_bed),T_nodes)+...
        sparse(1:length(T_bdy_nodes_bed),T_bdy_nodes_bed,2./(2*(T_Delta_z_cell(T_bdy_nodes_bed)).^2), length(T_bdy_nodes_bed),T_nodes);
    ddphidzbed_dT = sparse(length(T_bdy_nodes_bed),T_nodes);
    ddphidzbed_du = sparse(length(T_bdy_nodes_bed),T_nodes);
    ddphidzbed_dPi = sparse(length(T_bdy_nodes_bed),length(T_bdy_nodes_bed));
    ddphidzbed_dh = sparse(length(T_bdy_nodes_bed),1);

elseif strcmp(parameters.flag.sigma, 'spectral') == 1
    fout_phizz = conv_phizz(length(T_bdy_nodes_bed), parameters.grid.N.extra.bd_y, parameters.grid.N.coor_nodes.y(T_bdy_nodes_bed), h_prev);
    Gf_conv = fout_phizz.Gf_conv;
    Gfx_conv = fout_phizz.Gfx_conv;
    
    %define f, f_x and u_x; h should be h_i
    f = (gamma.*(1+gamma_pert).*f_slide_Tbed.*f_slide_hydro.*u_bed);
    df_du = sparse(1:length(T_bdy_nodes_bed),T_bdy_nodes_bed,(gamma.*(1+gamma_pert).*f_slide_Tbed.*f_slide_hydro), length(T_bdy_nodes_bed), T_nodes) + ...
        sparse(1:length(T_bdy_nodes_bed),T_bdy_nodes_bed,(gamma.*(1+gamma_pert).*f_slide_Tbed.*u_bed).*Df_slide_hydro_dub, length(T_bdy_nodes_bed), T_nodes);
    df_dT = sparse(1:length(T_bdy_nodes_bed),T_bdy_nodes_bed,(gamma.*(1+gamma_pert).*f_slide_hydro.* DfslideTbed.*u_bed), length(T_bdy_nodes_bed), T_nodes);
    df_dPi = sparse(1:length(T_bdy_nodes_bed),1:length(T_bdy_nodes_bed),(gamma.*(1+gamma_pert).*f_slide_Tbed.*u_bed).*Df_slide_hydro_dPi, length(T_bdy_nodes_bed), length(T_bdy_nodes_bed));
    
    f_x = (f-f_prev)./ dx_prev;
    u_x = (u_bed-u_prev(T_bdy_nodes_bed))./ dx;
    h_x = (h-h_prev)./dx;
    
    ddphi_dz_bed = Gf_conv*(f.*h_x) + Gfx_conv*f_x + 1/parameters.grid.N.extra.bd_y*sum(u_x.*T_Delta_y_cell(T_bdy_nodes_bed));
    
    ddphidz_df = Gf_conv;
    ddphidz_dfx = Gfx_conv;
    
    ddphidzbed_dT = ddphidz_df*df_dT + (ddphidz_dfx./dx_prev)*df_dT;
    ddphidzbed_dPi = ddphidz_df*df_dPi + (ddphidz_dfx./dx_prev)*df_dPi;
    ddphidzbed_du = ddphidz_df*df_du + (ddphidz_dfx./dx_prev)*df_du +  sparse(1:length(T_bdy_nodes_bed), T_bdy_nodes_bed,T_Delta_y_cell(T_bdy_nodes_bed)./(parameters.grid.N.extra.bd_y*dx),length(T_bdy_nodes_bed), T_nodes);
    dddphidzbed_dphi = sparse(length(T_bdy_nodes_bed),T_nodes);
    ddphidzbed_dh = Gf_conv*(f./dx) ;
end

index_fluxpsidown = circshift(1:length(psi_bdy_nodes_bed),-1,2);
index_fluxpsiup = 1:length(psi_bdy_nodes_bed);

dfluxpsibottom_dpsi_array = sparse(1:length(psi_bdy_nodes_bed),dfluxpsibottom_dpsibed_vec_index,dfluxpsibottom_dpsibed_vec, length(psi_bdy_nodes_bed),psi_nodes ) + ...
    sparse(1:length(psi_bdy_nodes_bed),dfluxpsibottom_dpsiabed_vec_index,dfluxpsibottom_dpsiabed_vec, length(psi_bdy_nodes_bed),psi_nodes );
dfluxpsibottom_dpsighost_array = dfluxpsibottom_dpsighost;

ddpsi_dz_dy = (flux_psi_bottom(index_fluxpsidown) - flux_psi_bottom(index_fluxpsiup))./psi_Delta_y_cell(psi_bdy_nodes_bed);
ddpsydzdy_dpsi = spdiags(1./psi_Delta_y_cell(psi_bdy_nodes_bed),0,length(psi_bdy_nodes_bed))*(dfluxpsibottom_dpsi_array(index_fluxpsidown,:) - dfluxpsibottom_dpsi_array(index_fluxpsiup,:));
ddpsydzdy_dpsighost = sparse(1:length(psi_bdy_nodes_bed),1:length(psi_bdy_nodes_bed), 1./psi_Delta_y_cell(psi_bdy_nodes_bed),length(psi_bdy_nodes_bed),length(psi_bdy_nodes_bed))...
    *(dfluxpsibottom_dpsighost_array(index_fluxpsidown,:) - dfluxpsibottom_dpsighost_array(index_fluxpsiup,:));

sigma = pr(T_bdy_nodes_bed) + 2*h_prev^(-1).*ddpsi_dz_dy - 2*h_prev^(-2).*ddphi_dz_bed;

dsigma_dp = dpr_dp(T_bdy_nodes_bed,:);
dsigma_du = dpr_du(T_bdy_nodes_bed,:)- 2*h_prev^(-2).* ddphidzbed_du;
dsigma_dT = - 2*h_prev^(-2).* ddphidzbed_dT;
dsigma_dPi = - 2*h_prev^(-2).* ddphidzbed_dPi;
dsigma_dphi = - 2*h_prev^(-2).*dddphidzbed_dphi;
dsigma_dpsi = 2*h_prev^(-1).*ddpsydzdy_dpsi;
dsigma_dpsighost =  2*h_prev^(-1).*ddpsydzdy_dpsighost;
dsigma_dh = - 2*h_prev^(-2).*ddphidzbed_dh;

%meltwater fluxes
%permeabilities
[kappa_fout, kappa_Dfout] = permeabilities(Pi, parameters);
kappa = kappa_fout.kappa;
dkappa_dPi = kappa_Dfout.kappa;

kappa_2 = kappa_fout.kappa_2;
dkappa2_dPi = kappa_Dfout.kappa2;

%qx = kappa.*(h-h_prev+r^(-1)*(b-b_prev))/dx;
dqx_dPi = -dkappa_dPi.*(h-h_prev+r^(-1)*(b-b_prev))./dx;
dqx_dh =  -kappa.*(1)/dx;

index_sigma_edge_down = circshift(1:length(sigma),-1,2);
index_sigma_edge_up = 1:length(sigma);

if (isfield(parameters.flag, 'flux_y') == 1 && strcmp(parameters.flag.flux_y, 'upwind') == 1)
    dsigma_dy = parameters.v_in_prev.dsigma_dy_prev; 
    ddsigmady_dp = sparse(length(T_bdy_nodes_bed),T_nodes); 
    ddsigmady_du = sparse(length(T_bdy_nodes_bed),T_nodes); 
    ddsigmady_dT = sparse(length(T_bdy_nodes_bed),T_nodes); 
    ddsigmady_dPi = sparse(length(T_bdy_nodes_bed),length(T_bdy_nodes_bed)); 
    ddsigmady_dphi = sparse(length(T_bdy_nodes_bed),T_nodes); 
    ddsigmady_dpsi = sparse(length(T_bdy_nodes_bed),psi_nodes); 
    ddsigmady_dpsighost = sparse(length(T_bdy_nodes_bed),1); 
    ddsigmady_dh = sparse(length(T_bdy_nodes_bed),1); 
else
    dsigma_dy = (sigma(index_sigma_edge_down)-sigma(index_sigma_edge_up))./T_Delta_y_cell(T_bdy_nodes_bed);
    ddsigmady_dp = spdiags(1./T_Delta_y_cell(T_bdy_nodes_bed),0,length(T_bdy_nodes_bed))*(dsigma_dp(index_sigma_edge_down,:) - dsigma_dp(index_sigma_edge_up,:));
    ddsigmady_dh = spdiags(1./T_Delta_y_cell(T_bdy_nodes_bed),0,length(T_bdy_nodes_bed))*(dsigma_dh(index_sigma_edge_down,:) - dsigma_dh(index_sigma_edge_up,:));
    ddsigmady_du = spdiags(1./T_Delta_y_cell(T_bdy_nodes_bed),0,length(T_bdy_nodes_bed))*(dsigma_du(index_sigma_edge_down,:) - dsigma_du(index_sigma_edge_up,:));
    ddsigmady_dT = spdiags(1./T_Delta_y_cell(T_bdy_nodes_bed),0,length(T_bdy_nodes_bed))*(dsigma_dT(index_sigma_edge_down,:) - dsigma_dT(index_sigma_edge_up,:));
    ddsigmady_dPi = spdiags(1./T_Delta_y_cell(T_bdy_nodes_bed),0,length(T_bdy_nodes_bed))*(dsigma_dPi(index_sigma_edge_down,:) - dsigma_dPi(index_sigma_edge_up,:));
    ddsigmady_dphi = spdiags(1./T_Delta_y_cell(T_bdy_nodes_bed),0,length(T_bdy_nodes_bed))*(dsigma_dphi(index_sigma_edge_down,:) - dsigma_dphi(index_sigma_edge_up,:));
    ddsigmady_dpsi = spdiags(1./T_Delta_y_cell(T_bdy_nodes_bed),0,length(T_bdy_nodes_bed))*(dsigma_dpsi(index_sigma_edge_down,:) - dsigma_dpsi(index_sigma_edge_up,:));
    ddsigmady_dpsighost = spdiags(1./T_Delta_y_cell(T_bdy_nodes_bed),0,length(T_bdy_nodes_bed))*(dsigma_dpsighost(index_sigma_edge_down,:) - dsigma_dpsighost(index_sigma_edge_up,:));

end

dPi_dy = (Pi(index_sigma_edge_down)-Pi(index_sigma_edge_up))./T_Delta_y_cell(T_bdy_nodes_bed);
ddPidy_dPi = sparse(1:length(T_bdy_nodes_bed),index_sigma_edge_down,1./T_Delta_y_cell(T_bdy_nodes_bed), length(T_bdy_nodes_bed),length(T_bdy_nodes_bed))+...
    sparse(1:length(T_bdy_nodes_bed),index_sigma_edge_up,-1./T_Delta_y_cell(T_bdy_nodes_bed), length(T_bdy_nodes_bed),length(T_bdy_nodes_bed));

%qy =  -I_prev_edge.*(kappa_edge.*dsigma_dy + beta.*kappa_2_edge.*dPi_dy);
I_prev_edge = 1/2*(I_prev(index_sigma_edge_up) + I_prev(index_sigma_edge_down));
I_prev_edge(I_prev_edge<1) = 0;

if isfield(parameters.flag, 'flux_y') == 0
    
    Pi_edge = Pi(index_sigma_edge_up);
    dPi_edge_dPi_up = ones(length(Pi),1);
    dPi_edge_dPi_down = zeros(length(Pi),1);
    
    [kappa_edge_fout, kappa_edge_Dfout] = permeabilities(Pi_edge, parameters);
    kappa_edge = kappa_edge_fout.kappa;
    kappa_2_edge = kappa_edge_fout.kappa_2;
    dkappa_edge_dPiedge = kappa_edge_Dfout.kappa;
    dkappa2_edge_dPiedge = kappa_edge_Dfout.kappa2;
    
    dkappa_edge_dPi_up = dkappa_edge_dPiedge.*dPi_edge_dPi_up;
    dkappa_edge_dPi_down = dkappa_edge_dPiedge.*dPi_edge_dPi_down;
    
    dkappa2_edge_dPi_up = dkappa2_edge_dPiedge.*dPi_edge_dPi_up;
    dkappa2_edge_dPi_down = dkappa2_edge_dPiedge.*dPi_edge_dPi_down;
    
elseif strcmp(parameters.flag.flux_y, 'centered') == 1
    Pi_edge = 1/2*(Pi(index_sigma_edge_up) + Pi(index_sigma_edge_down));
    dPi_edge_dPi_up = 1/2*ones(length(Pi),1);
    dPi_edge_dPi_down = 1/2*ones(length(Pi),1);
    
    [kappa_edge_fout, kappa_edge_Dfout] = permeabilities(Pi_edge, parameters);
    kappa_edge = kappa_edge_fout.kappa;
    kappa_2_edge = kappa_edge_fout.kappa_2;
    dkappa_edge_dPiedge = kappa_edge_Dfout.kappa;
    dkappa2_edge_dPiedge = kappa_edge_Dfout.kappa2;
    
    dkappa_edge_dPi_up = dkappa_edge_dPiedge.*dPi_edge_dPi_up;
    dkappa_edge_dPi_down = dkappa_edge_dPiedge.*dPi_edge_dPi_down;
    
    dkappa2_edge_dPi_up = dkappa2_edge_dPiedge.*dPi_edge_dPi_up;
    dkappa2_edge_dPi_down = dkappa2_edge_dPiedge.*dPi_edge_dPi_down;
    
    
elseif strcmp(parameters.flag.flux_y, 'centered_cutoff') == 1
    Pi_edge = 1/2*(Pi(index_sigma_edge_up) + Pi(index_sigma_edge_down));
    dPi_edge_dPi_up = 1/2*ones(length(Pi),1);
    dPi_edge_dPi_down = 1/2*ones(length(Pi),1);
    
    [kappa_edge_fout, kappa_edge_Dfout] = permeabilities(Pi_edge, parameters);
    index_p = find(kappa_edge_fout.kappa>0);
    
    kappa_edge = max(kappa_edge_fout.kappa,0);
    dkappa_edge_dPiedge = zeros(length(kappa_edge),1);
    dkappa_edge_dPiedge(index_p)= kappa_edge_Dfout.kappa(index_p);
    
    dkappa_edge_dPi_up = dkappa_edge_dPiedge.*dPi_edge_dPi_up;
    dkappa_edge_dPi_down = dkappa_edge_dPiedge.*dPi_edge_dPi_down;
    
    kappa_2_edge = kappa_edge_fout.kappa_2;
    dkappa2_edge_dPiedge = kappa_edge_Dfout.kappa2;
    
    dkappa2_edge_dPi_up = dkappa2_edge_dPiedge.*dPi_edge_dPi_up;
    dkappa2_edge_dPi_down = dkappa2_edge_dPiedge.*dPi_edge_dPi_down;
    
elseif strcmp(parameters.flag.flux_y, 'centered_harmonic') == 1
    
    Pi_up = Pi(index_sigma_edge_up);
    Pi_down = Pi(index_sigma_edge_down);
    [kappa_up_fout, kappa_up_Dfout] = permeabilities(Pi_up, parameters);
    [kappa_down_fout, kappa_down_Dfout] = permeabilities(Pi_down, parameters);
    
    kappa_up = max(kappa_up_fout.kappa,0);
    index_z_up  = find(kappa_up==0);
    dkappaup_dPi_up = kappa_up_Dfout.kappa;
    dkappaup_dPi_up(index_z_up) = 0;
    
    kappa_down =  max(kappa_down_fout.kappa,0);
    index_z_down  = find(kappa_down==0);
    dkappadown_dPi_down = kappa_down_Dfout.kappa;
    dkappadown_dPi_down(index_z_down) = 0;

    kappa_edge = 2.* kappa_up.*kappa_down./(kappa_up + kappa_down);
    kappa_edge(kappa_up.*kappa_down==0) = 0;
    
    dkappa_edge_dkappaup = (2.*kappa_down.*(kappa_up + kappa_down) - 2.* kappa_up.*kappa_down )./(kappa_up + kappa_down).^2;
    dkappa_edge_dkappaup(2.*kappa_down.*(kappa_up + kappa_down) - 2.* kappa_up.*kappa_down == 0) = 0;
    
    dkappa_edge_dkappadown = (2.*kappa_up.*(kappa_up + kappa_down) - 2.* kappa_up.*kappa_down )./(kappa_up + kappa_down).^2;
    dkappa_edge_dkappadown(2.*kappa_up.*(kappa_up + kappa_down) - 2.* kappa_up.*kappa_down == 0) = 0;
    
    dkappa_edge_dPi_up = dkappa_edge_dkappaup.*dkappaup_dPi_up;
    dkappa_edge_dPi_down = dkappa_edge_dkappadown.*dkappadown_dPi_down;
    
    kappa_2_edge = kappa_up_fout.kappa_2;
    dkappa2_edge_dPi_up = kappa_up_Dfout.kappa2.*ones(length(Pi_up),1);
    dkappa2_edge_dPi_down = kappa_down_Dfout.kappa2.*ones(length(Pi_down),1);
    
elseif strcmp(parameters.flag.flux_y, 'upwind') == 1
    
    adv_speed = -parameters.v_in_prev.dsigma_dy_prev; 
    index_p = find(adv_speed>=0);
    index_n = find(adv_speed<0);
    
    index_Pi_edge = zeros(length(index_p),1);
    index_Pi_edge(index_p) = index_sigma_edge_up(index_p);
    index_Pi_edge(index_n) = index_sigma_edge_down(index_n);
    
    Pi_edge = Pi(index_Pi_edge);
    dPi_edge_dPi_up = zeros(length(Pi),1);
    dPi_edge_dPi_up(index_p) = 1;
    dPi_edge_dPi_down = zeros(length(Pi),1);
    dPi_edge_dPi_down(index_n) = 1;
    
    [kappa_edge_fout, kappa_edge_Dfout] = permeabilities(Pi_edge, parameters);
    
    kappa_edge = kappa_edge_fout.kappa;
    kappa_2_edge = kappa_edge_fout.kappa_2;
    dkappa_edge_dPiedge = kappa_edge_Dfout.kappa;
    dkappa2_edge_dPiedge = kappa_edge_Dfout.kappa2;
    
    dkappa_edge_dPi_up = dkappa_edge_dPiedge.*dPi_edge_dPi_up;
    dkappa_edge_dPi_down = dkappa_edge_dPiedge.*dPi_edge_dPi_down;
    
    dkappa2_edge_dPi_up = dkappa2_edge_dPiedge.*dPi_edge_dPi_up;
    dkappa2_edge_dPi_down = dkappa2_edge_dPiedge.*dPi_edge_dPi_down;
    
    if isempty(index_p) == 1
         dkappa2_edge_dPi_up = sparse(length(Pi),1);
         dkappa_edge_dPi_up = sparse(length(Pi),1);
    elseif isempty(index_n) == 1
         dkappa2_edge_dPi_down = sparse(length(Pi),1);
         dkappa_edge_dPi_down = sparse(length(Pi),1);
    end
      
end


if isfield(parameters.flag, 'fluxy_sigma') == 0
    ddqy_dPi = -sparse(1:length(T_bdy_nodes_bed), 1:length(T_bdy_nodes_bed), I_prev_edge,length(T_bdy_nodes_bed), length(T_bdy_nodes_bed))*...
        (sparse(1:length(T_bdy_nodes_bed),index_sigma_edge_up, dkappa_edge_dPi_up.*dsigma_dy+ beta.*dkappa2_edge_dPi_up.*dPi_dy, length(T_bdy_nodes_bed), length(T_bdy_nodes_bed)) +...
        sparse(1:length(T_bdy_nodes_bed),index_sigma_edge_down, dkappa_edge_dPi_down.*dsigma_dy+ beta.*dkappa2_edge_dPi_down.*dPi_dy, length(T_bdy_nodes_bed), length(T_bdy_nodes_bed))+...
        sparse(1:length(T_bdy_nodes_bed), 1:length(T_bdy_nodes_bed), beta.*kappa_2_edge, length(T_bdy_nodes_bed), length(T_bdy_nodes_bed))*ddPidy_dPi + ...
        sparse(1:length(T_bdy_nodes_bed), 1:length(T_bdy_nodes_bed), kappa_edge, length(T_bdy_nodes_bed), length(T_bdy_nodes_bed))*ddsigmady_dPi);
    ddqy_dT = -sparse(1:length(T_bdy_nodes_bed), 1:length(T_bdy_nodes_bed),I_prev_edge.*kappa_edge,length(T_bdy_nodes_bed), length(T_bdy_nodes_bed))*ddsigmady_dT;
    ddqy_dp = -sparse(1:length(T_bdy_nodes_bed), 1:length(T_bdy_nodes_bed),I_prev_edge.*kappa_edge,length(T_bdy_nodes_bed), length(T_bdy_nodes_bed))*ddsigmady_dp;
    ddqy_du = -sparse(1:length(T_bdy_nodes_bed), 1:length(T_bdy_nodes_bed), I_prev_edge.*kappa_edge,length(T_bdy_nodes_bed), length(T_bdy_nodes_bed))*ddsigmady_du;
    ddqy_dpsi = -sparse(1:length(T_bdy_nodes_bed), 1:length(T_bdy_nodes_bed), I_prev_edge.*kappa_edge,length(T_bdy_nodes_bed), length(T_bdy_nodes_bed))*ddsigmady_dpsi;
    ddqy_dphi = -sparse(1:length(T_bdy_nodes_bed), 1:length(T_bdy_nodes_bed), I_prev_edge.*kappa_edge,length(T_bdy_nodes_bed), length(T_bdy_nodes_bed))*ddsigmady_dphi;
    ddqy_dpsighost = -sparse(1:length(T_bdy_nodes_bed), 1:length(T_bdy_nodes_bed), I_prev_edge.*kappa_edge,length(T_bdy_nodes_bed), length(T_bdy_nodes_bed))*ddsigmady_dpsighost;
    ddqy_dh = -sparse(1:length(T_bdy_nodes_bed), 1:length(T_bdy_nodes_bed), I_prev_edge.*kappa_edge,length(T_bdy_nodes_bed), length(T_bdy_nodes_bed))*ddsigmady_dh;
     
    
elseif strcmp(parameters.flag.fluxy_sigma, 'full')== 1
    ddqy_dPi = -sparse(1:length(T_bdy_nodes_bed), 1:length(T_bdy_nodes_bed), I_prev_edge,length(T_bdy_nodes_bed), length(T_bdy_nodes_bed))*...
        (sparse(1:length(T_bdy_nodes_bed),index_sigma_edge_up, dkappa_edge_dPi_up.*dsigma_dy+ beta.*dkappa2_edge_dPi_up.*dPi_dy, length(T_bdy_nodes_bed), length(T_bdy_nodes_bed)) +...
        sparse(1:length(T_bdy_nodes_bed),index_sigma_edge_down, dkappa_edge_dPi_down.*dsigma_dy+ beta.*dkappa2_edge_dPi_down.*dPi_dy, length(T_bdy_nodes_bed), length(T_bdy_nodes_bed))+...
        sparse(1:length(T_bdy_nodes_bed), 1:length(T_bdy_nodes_bed), beta.*kappa_2_edge, length(T_bdy_nodes_bed), length(T_bdy_nodes_bed))*ddPidy_dPi + ...
        sparse(1:length(T_bdy_nodes_bed), 1:length(T_bdy_nodes_bed), kappa_edge, length(T_bdy_nodes_bed), length(T_bdy_nodes_bed))*ddsigmady_dPi);
    ddqy_dT = -sparse(1:length(T_bdy_nodes_bed), 1:length(T_bdy_nodes_bed),I_prev_edge.*kappa_edge,length(T_bdy_nodes_bed), length(T_bdy_nodes_bed))*ddsigmady_dT;
    ddqy_dp = -sparse(1:length(T_bdy_nodes_bed), 1:length(T_bdy_nodes_bed),I_prev_edge.*kappa_edge,length(T_bdy_nodes_bed), length(T_bdy_nodes_bed))*ddsigmady_dp;
    ddqy_du = -sparse(1:length(T_bdy_nodes_bed), 1:length(T_bdy_nodes_bed), I_prev_edge.*kappa_edge,length(T_bdy_nodes_bed), length(T_bdy_nodes_bed))*ddsigmady_du;
    ddqy_dpsi = -sparse(1:length(T_bdy_nodes_bed), 1:length(T_bdy_nodes_bed), I_prev_edge.*kappa_edge,length(T_bdy_nodes_bed), length(T_bdy_nodes_bed))*ddsigmady_dpsi;
    ddqy_dphi = -sparse(1:length(T_bdy_nodes_bed), 1:length(T_bdy_nodes_bed), I_prev_edge.*kappa_edge,length(T_bdy_nodes_bed), length(T_bdy_nodes_bed))*ddsigmady_dphi;
    ddqy_dpsighost = -sparse(1:length(T_bdy_nodes_bed), 1:length(T_bdy_nodes_bed), I_prev_edge.*kappa_edge,length(T_bdy_nodes_bed), length(T_bdy_nodes_bed))*ddsigmady_dpsighost;
    ddqy_dh = -sparse(1:length(T_bdy_nodes_bed), 1:length(T_bdy_nodes_bed), I_prev_edge.*kappa_edge,length(T_bdy_nodes_bed), length(T_bdy_nodes_bed))*ddsigmady_dh;
    
elseif strcmp(parameters.flag.fluxy_sigma, 'Pionly')== 1
    ddqy_dPi = -sparse(1:length(T_bdy_nodes_bed), 1:length(T_bdy_nodes_bed), I_prev_edge,length(T_bdy_nodes_bed), length(T_bdy_nodes_bed))*...
        (sparse(1:length(T_bdy_nodes_bed),index_sigma_edge_up, beta.*dkappa2_edge_dPi_up.*dPi_dy, length(T_bdy_nodes_bed), length(T_bdy_nodes_bed)) +...
        sparse(1:length(T_bdy_nodes_bed),index_sigma_edge_down, beta.*dkappa2_edge_dPi_down.*dPi_dy, length(T_bdy_nodes_bed), length(T_bdy_nodes_bed))+...
        sparse(1:length(T_bdy_nodes_bed), 1:length(T_bdy_nodes_bed), beta.*kappa_2_edge, length(T_bdy_nodes_bed), length(T_bdy_nodes_bed))*ddPidy_dPi);
    ddqy_dT = sparse(length(T_bdy_nodes_bed), T_nodes);
    ddqy_dp = sparse(length(T_bdy_nodes_bed), T_nodes);
    ddqy_du = sparse(length(T_bdy_nodes_bed), T_nodes);
    ddqy_dpsi = sparse(length(T_bdy_nodes_bed), psi_nodes);
    ddqy_dphi = sparse(length(T_bdy_nodes_bed), T_nodes);
    ddqy_dpsighost = sparse(length(T_bdy_nodes_bed), 1);
    ddqy_dh = sparse(length(T_bdy_nodes_bed), 1);
    
end

%% BASAL ENERGY BUDGET
%this is effectively a constraint on basal temperature, and should be
%treated as an additional equation in the unknown T_bed

%m = -(bedflux_above - bedflux_below) + alpha*(gamma.*(1+gamma_pert).*f_slide_Tbed.*f_slide_hydro).*u(T_bdy_nodes_bed).^(2);

%bedflux_above = - (T(T_bdy_nodes_bed-T_nodes_hor(end))- T(T_bdy_nodes_bed))./(h_prev*T_Delta_z_cell(T_bdy_nodes_bed));   
dbedfluxabove_dT = sparse(1:length(T_bdy_nodes_bed),T_bdy_nodes_bed-T_nodes_hor(end), -1./(h_prev*T_Delta_z_cell(T_bdy_nodes_bed)), length(T_bdy_nodes_bed), T_nodes)+ ...
    sparse(1:length(T_bdy_nodes_bed),T_bdy_nodes_bed, 1./(h_prev*T_Delta_z_cell(T_bdy_nodes_bed)), length(T_bdy_nodes_bed), T_nodes);

%bedflux_below = - (T(T_bdy_nodes_bed) - Tb(Tb_bdy_nodes_top))./(Tb_Delta_z_cell(Tb_bdy_nodes_top));   
dbedfluxbelow_dTb = sparse(1:length(Tb_bdy_nodes_top),Tb_bdy_nodes_top, 1./(Tb_Delta_z_cell(Tb_bdy_nodes_top)), length(Tb_bdy_nodes_top), Tb_nodes);
dbedfluxbelow_dT = sparse(1:length(Tb_bdy_nodes_top),T_bdy_nodes_bed,-1./(Tb_Delta_z_cell(Tb_bdy_nodes_top)), length(Tb_bdy_nodes_top), T_nodes);

dfriction_dT = sparse(1:length(T_bdy_nodes_bed),T_bdy_nodes_bed,  alpha*(gamma.*(1+gamma_pert)).*f_slide_hydro.*u(T_bdy_nodes_bed).^(2).*DfslideTbed,length(T_bdy_nodes_bed),T_nodes);
    dfriction_du = sparse(1:length(T_bdy_nodes_bed),T_bdy_nodes_bed,alpha*(gamma.*(1+gamma_pert).*f_slide_hydro.*f_slide_Tbed).*2.*u(T_bdy_nodes_bed) + ...
        alpha*(gamma.*(1+gamma_pert).*f_slide_Tbed.*Df_slide_hydro_dub).*u(T_bdy_nodes_bed).^(2),length(T_bdy_nodes_bed),T_nodes);
    dfriction_dPi = sparse(1:length(T_bdy_nodes_bed),1:length(T_bdy_nodes_bed),  alpha*(gamma.*(1+gamma_pert).*f_slide_Tbed.*Df_slide_hydro_dPi).*u(T_bdy_nodes_bed).^(2),length(T_bdy_nodes_bed),length(T_bdy_nodes_bed));
    dfriction_dh = 0;

dm_du = dfriction_du;
dm_dTb = -(-dbedfluxbelow_dTb);
dm_dT = -(dbedfluxabove_dT - dbedfluxbelow_dT) + dfriction_dT;
dm_dPi = dfriction_dPi;
dm_dh = dfriction_dh;

index_sigma_down = 1:length(sigma);
index_sigma_up = circshift(1:length(sigma),1,2);
%div_qy = (qy(index_sigma_down)-qy(index_sigma_up))./(T_Delta_y_cell(T_bdy_nodes_bed));
ddivqy_dPi = sparse(1:length(T_bdy_nodes_bed),1:length(T_bdy_nodes_bed),1./(T_Delta_y_cell(T_bdy_nodes_bed)),length(T_bdy_nodes_bed),length(T_bdy_nodes_bed))*(ddqy_dPi(index_sigma_down,:) - ddqy_dPi(index_sigma_up,:));
ddivqy_dT = sparse(1:length(T_bdy_nodes_bed),1:length(T_bdy_nodes_bed),1./(T_Delta_y_cell(T_bdy_nodes_bed)),length(T_bdy_nodes_bed),length(T_bdy_nodes_bed))*(ddqy_dT(index_sigma_down,:) - ddqy_dT(index_sigma_up,:));
ddivqy_dp = sparse(1:length(T_bdy_nodes_bed),1:length(T_bdy_nodes_bed),1./(T_Delta_y_cell(T_bdy_nodes_bed)),length(T_bdy_nodes_bed),length(T_bdy_nodes_bed))*(ddqy_dp(index_sigma_down,:) - ddqy_dp(index_sigma_up,:));
ddivqy_du = sparse(1:length(T_bdy_nodes_bed),1:length(T_bdy_nodes_bed),1./(T_Delta_y_cell(T_bdy_nodes_bed)),length(T_bdy_nodes_bed),length(T_bdy_nodes_bed))*(ddqy_du(index_sigma_down,:) - ddqy_du(index_sigma_up,:));
ddivqy_dpsi = sparse(1:length(T_bdy_nodes_bed),1:length(T_bdy_nodes_bed),1./(T_Delta_y_cell(T_bdy_nodes_bed)),length(T_bdy_nodes_bed),length(T_bdy_nodes_bed))*(ddqy_dpsi(index_sigma_down,:) - ddqy_dpsi(index_sigma_up,:));
ddivqy_dphi = sparse(1:length(T_bdy_nodes_bed),1:length(T_bdy_nodes_bed),1./(T_Delta_y_cell(T_bdy_nodes_bed)),length(T_bdy_nodes_bed),length(T_bdy_nodes_bed))*(ddqy_dphi(index_sigma_down,:) - ddqy_dphi(index_sigma_up,:));
ddivqy_dpsighost = sparse(1:length(T_bdy_nodes_bed),1:length(T_bdy_nodes_bed),1./(T_Delta_y_cell(T_bdy_nodes_bed)),length(T_bdy_nodes_bed),length(T_bdy_nodes_bed))*(ddqy_dpsighost(index_sigma_down,:) - ddqy_dpsighost(index_sigma_up,:));
ddivqy_dh = sparse(1:length(T_bdy_nodes_bed),1:length(T_bdy_nodes_bed),1./(T_Delta_y_cell(T_bdy_nodes_bed)),length(T_bdy_nodes_bed),length(T_bdy_nodes_bed))*(ddqy_dh(index_sigma_down,:) - ddqy_dh(index_sigma_up,:));

%conservation law
%enth_eq = (qx-qx_prev)/dx_int + div_qy - m;
dentheq_dPi = (dx_int)^(-1).*sparse(1:length(T_bdy_nodes_bed), 1:length(T_bdy_nodes_bed), I_prev,length(T_bdy_nodes_bed), length(T_bdy_nodes_bed))...
    *sparse(1:length(T_bdy_nodes_bed),1:length(T_bdy_nodes_bed),dqx_dPi,length(T_bdy_nodes_bed),length(T_bdy_nodes_bed)) + ddivqy_dPi -dm_dPi;
dentheq_dh = (dx_int)^(-1).*sparse(1:length(T_bdy_nodes_bed), 1:length(T_bdy_nodes_bed), I_prev,length(T_bdy_nodes_bed), length(T_bdy_nodes_bed))*dqx_dh + ddivqy_dh - dm_dh;
dentheq_dp = ddivqy_dp;
dentheq_du =  - dm_du + ddivqy_du;
dentheq_dpsi = ddivqy_dpsi;
dentheq_dphi = ddivqy_dphi;
dentheq_dpsighost = ddivqy_dpsighost;
dentheq_dT = ddivqy_dT -dm_dT;
dentheq_dTb = -dm_dTb;

%operator split
%constrain Pi where bed was temperate at previous timestep
%fout(Pi_index(index_temp_prev)) = enth_eq(index_temperate_prev);
fout(Pi_index(index_temperate_prev), Pi_index) = dentheq_dPi(index_temperate_prev,:);
fout(Pi_index(index_temperate_prev), h_index) = dentheq_dh(index_temperate_prev,:);
fout(Pi_index(index_temperate_prev), p_index) = dentheq_dp(index_temperate_prev,:);
fout(Pi_index(index_temperate_prev), u_index) = dentheq_du(index_temperate_prev,:);
fout(Pi_index(index_temperate_prev), psi_index) = dentheq_dpsi(index_temperate_prev,:);
fout(Pi_index(index_temperate_prev), phi_index) = dentheq_dphi(index_temperate_prev,:);
fout(Pi_index(index_temperate_prev), psighost_index) = dentheq_dpsighost(index_temperate_prev,:);
fout(Pi_index(index_temperate_prev), T_index) = dentheq_dT(index_temperate_prev,:);
fout(Pi_index(index_temperate_prev), Tb_index) = dentheq_dTb(index_temperate_prev,:);

%constrain T_bed where bed was frozen at current timestep
%fout(T_index(T_bdy_nodes_bed(index_frozen_prev))) = enth_eq(index_frozen_prev);
fout(T_index(T_bdy_nodes_bed(index_frozen_prev)), Pi_index) = dentheq_dPi(index_frozen_prev,:);
fout(T_index(T_bdy_nodes_bed(index_frozen_prev)), h_index) = dentheq_dh(index_frozen_prev,:);
fout(T_index(T_bdy_nodes_bed(index_frozen_prev)), p_index) = dentheq_dp(index_frozen_prev,:);
fout(T_index(T_bdy_nodes_bed(index_frozen_prev)), u_index) = dentheq_du(index_frozen_prev,:);
fout(T_index(T_bdy_nodes_bed(index_frozen_prev)), psi_index) = dentheq_dpsi(index_frozen_prev,:);
fout(T_index(T_bdy_nodes_bed(index_frozen_prev)), phi_index) = dentheq_dphi(index_frozen_prev,:);
fout(T_index(T_bdy_nodes_bed(index_frozen_prev)), psighost_index) = dentheq_dpsighost(index_frozen_prev,:);
fout(T_index(T_bdy_nodes_bed(index_frozen_prev)), T_index) = dentheq_dT(index_frozen_prev,:);
fout(T_index(T_bdy_nodes_bed(index_frozen_prev)), Tb_index) = dentheq_dTb(index_frozen_prev,:);
