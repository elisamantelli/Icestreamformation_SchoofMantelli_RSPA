function [fout, faux] = network_timestep_v5_mech_Tfixed_jacobian_v2(v_in,parameters)
%Jacobian of network_timestep_v5_mech_Tfixed_v2.,
%Tested against numerical jacobian. Elisa Mantelli, April 2020

%parameters 
gamma = parameters.gamma;

%psi grid
psi_nodes = parameters.grid.psi.n_nodes.tot;                                   %number of nodes
psi_nodes_hor = parameters.grid.psi.n_nodes.hor;                               %number of nodes
psi_nedges_ver = parameters.grid.psi.n_edges.vert;

psi_bdy_nodes_top = parameters.grid.psi.bdy_nodes.top;                         %list of nodes adjacent to the top of the box;
psi_bdy_nodes_bed = parameters.grid.psi.bdy_nodes.bed;

psi_Delta_z_cell = parameters.grid.psi.Delta_z_cell;                           %length of cells, ver (list)
psi_Delta_y_cell = parameters.grid.psi.Delta_y_cell;                           %length of cells, hor (list)
psi_Delta_y_partial_cell = parameters.grid.psi.delta_y_vflux;                  %length of hor cell esges crossed by vertical network edges, hor (list)

index_up_node_top_psi = 1:length(psi_bdy_nodes_top);
index_down_node_top_psi = circshift(1:length(psi_bdy_nodes_top), -1,2);

index_up_node_bed_psi = 1:length(psi_bdy_nodes_bed);
index_down_node_bed_psi = circshift(1:length(psi_bdy_nodes_bed), -1,2);

psi_connect_ver =  parameters.grid.psi.connect_ver;
psi_connect_hor = parameters.grid.psi.connect_hor;

%T grid
T_nodes = parameters.grid.N.n_nodes.tot;                                   %number of nodes
T_nedges_ver = parameters.grid.N.n_edges.vert;

T_bdy_nodes_top = parameters.grid.N.bdy_nodes.top;                         %list of nodes adjacent to the top of the box;
T_bdy_nodes_bed = parameters.grid.N.bdy_nodes.bottom;

T_Delta_z_cell_volume = parameters.grid.N.Delta_z_cell_volume;             %length of cells, ver (list)
T_Delta_y_cell = parameters.grid.N.Delta_y_cell;                           %length of cells, hor (list)
T_Delta_y_partial_cell = parameters.grid.N.delta_y_vflux;                  %length of hor cell esges crossed by vertical network edges, hor (list)
T_coor_hor_cell_edges_z = parameters.grid.N.coor_hor_cell_edges.z;

index_up_node_top_T = circshift(1:length(T_bdy_nodes_top),1);
index_down_node_top_T = 1:length(T_bdy_nodes_top);

T_connect_ver =  parameters.grid.N.connect_ver;
T_connect_hor = parameters.grid.N.connect_hor;

Tb_nodes = parameters.grid.Tb.n_nodes.tot; 

psi_index = 1:psi_nodes;
omega_index = psi_nodes+1: 2*psi_nodes;
phi_index = 2*psi_nodes+1:2*psi_nodes+T_nodes;
u_index = 2*psi_nodes+T_nodes+1:2*psi_nodes+2*T_nodes;
p_index = 2*psi_nodes+2*T_nodes+1: 2*psi_nodes+3*T_nodes;
psighost_index = 2*psi_nodes+3*T_nodes+1;

psi = v_in(psi_index);
phi = v_in(phi_index);
u = v_in(u_index);
h = parameters.h;
psi_ghost = v_in(psighost_index);

%unpack variable at previous timestep
h_prev = parameters.v_in_prev.h_prev;
u_prev = parameters.v_in_prev.u_prev;

%initialize output
fout = sparse(length(v_in),length(v_in));

%%  PRELIMINARIES 

%DISCRETE DEP. VARIABLES
v_in_disc = [v_in(1:p_index(end)); sparse(T_nodes+Tb_nodes,1); psi_ghost; h; 0; sparse(length(T_bdy_nodes_bed),1)];
[discvar, Ddiscvar] = discretisation_hanging_full_v4(v_in_disc, parameters);

ddpsidz_dpsi = Ddiscvar.ddpsidz_dpsi;
ddpsidy_dpsi = Ddiscvar.ddpsidy_dpsi;

ddomegadz_domega = Ddiscvar.ddomegadz_domega;
ddomegady_domega = Ddiscvar.ddomegady_domega;

u_vert = discvar.u_vert;
duvert_du = Ddiscvar.duvert_du;

ddudy_du = Ddiscvar.ddudy_du;

du_dz = discvar.du_dz; 
ddudz_du = Ddiscvar.ddudz_du;

ddphidy_dphi = Ddiscvar.ddphidy_dphi;
ddphidz_dphi = Ddiscvar.ddphidz_dphi;

ddpdy_dp = Ddiscvar.ddpdy_dp;
ddpdz_dp = Ddiscvar.ddpdz_dp;

%dT_dy = discvar.dT_dy;

%construct regularized bedwater content
T_bed = parameters.T_bed;
gamma_pert = parameters.gamma_pert;

Tbed_upnode = circshift(1:length(T_bdy_nodes_bed),1,2);
Tbed_downnode = 1:length(T_bdy_nodes_bed);
gamma_pert_psigrid = (gamma_pert(Tbed_upnode)+gamma_pert(Tbed_downnode))./2;
T_bed_psigrid = (T_bed(Tbed_upnode)+T_bed(Tbed_downnode))./2;


%TIMESTEPPING
x_current = parameters.timestep.x;
x_prev = x_current - (1/2*parameters.timestep.dx + 1/2*parameters.timestep.dx_prev);

dx_int = parameters.timestep.dx_prev; % x_{i+1/2} - x_{i-1/2}
dx = x_current -x_prev; % x_{i+1} - x_{i}

h_av = (h*parameters.timestep.dx_prev+h_prev*parameters.timestep.dx)./(parameters.timestep.dx_prev + parameters.timestep.dx); %eq. 16
dhav_dh = parameters.timestep.dx_prev./(parameters.timestep.dx_prev + parameters.timestep.dx)*ones(length(h),1);

%construct regularized bedwater content
if parameters.flag_Tdep == 1
    [f_slide_Tbed, DfslideTbed] = regularizedfriction_temperature(T_bed, parameters);
    [f_slide_Tbedpsigrid, DfslideTbedpsigrid] = regularizedfriction_temperature(T_bed_psigrid, parameters);
else
    f_slide_Tbed = ones(length(T_bdy_nodes_bed),1);
    DfslideTbed = zeros(length(T_bdy_nodes_bed),1);
    
    f_slide_Tbedpsigrid = ones(length(psi_bdy_nodes_bed),1);
    DfslideTbedpsigrid = zeros(length(psi_bdy_nodes_bed),1);
end

 %% STREAM FUNCTION, eq. 23 with bcs 28, 30
%net_psi_horflux = psi_connect_hor*dpsi_dy;
dpsihorflux_dpsi = psi_connect_hor*ddpsidy_dpsi;

%net_psi_verflux = psi_connect_ver*(h_prev^(-1).*dpsi_dz.*psi_Delta_y_partial_cell);
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

%at psi cell centres
index_up_node_bed = circshift(1:length(T_bdy_nodes_bed),1);
index_down_node_bed = 1:length(T_bdy_nodes_bed);
dphi_dy_bed = -(phi(T_bdy_nodes_bed(index_up_node_bed))-phi(T_bdy_nodes_bed(index_down_node_bed)))./psi_Delta_y_cell(psi_bdy_nodes_bed);
ddphidybed_dphiup =  -1./psi_Delta_y_cell(psi_bdy_nodes_bed);
index_ddphidybed_dphiup = T_bdy_nodes_bed(index_up_node_bed);
ddphidybed_dphidown =  1./psi_Delta_y_cell(psi_bdy_nodes_bed);
index_ddphidybed_dphidown = T_bdy_nodes_bed(index_down_node_bed);

domegabed_dfslide = gamma.*(1+gamma_pert_psigrid).*(dphi_dy_bed + flux_psi_bottom); %then still need to multiply by derivatives wrt T_bed
domegabed_ddphidybed =  (gamma.*(1+gamma_pert_psigrid).*f_slide_Tbedpsigrid);
domegabed_dfluxpsibottom =  (gamma.*(1+gamma_pert_psigrid).*f_slide_Tbedpsigrid);

%% additional condition for psi_ghost, int_width omega_bed = 0
%fout(psighost_index) = sum(omega_bed.*psi_Delta_y_cell(psi_bdy_nodes_bed)); 

fout(psighost_index,phi_index) = sum(sparse(1:length(psi_bdy_nodes_bed), index_ddphidybed_dphiup, psi_Delta_y_cell(psi_bdy_nodes_bed).*domegabed_ddphidybed.*(ddphidybed_dphiup), length(psi_bdy_nodes_bed),T_nodes)+...
    sparse(1:length(psi_bdy_nodes_bed), index_ddphidybed_dphidown, psi_Delta_y_cell(psi_bdy_nodes_bed).*domegabed_ddphidybed.*(ddphidybed_dphidown), length(psi_bdy_nodes_bed),T_nodes),1) ;
fout(psighost_index,psi_index) = sum(sparse(1:length(psi_bdy_nodes_bed),1:length(psi_bdy_nodes_bed), domegabed_dfluxpsibottom.*psi_Delta_y_cell(psi_bdy_nodes_bed),length(psi_bdy_nodes_bed),length(psi_bdy_nodes_bed))*dfluxpsibottom_dpsi,1);
fout(psighost_index,psighost_index) = sum(sparse(1:length(psi_bdy_nodes_bed),1:length(psi_bdy_nodes_bed), domegabed_dfluxpsibottom.*psi_Delta_y_cell(psi_bdy_nodes_bed),length(psi_bdy_nodes_bed),length(psi_bdy_nodes_bed))*...
    sparse(1:length(psi_bdy_nodes_bed),1:length(psi_ghost),dfluxpsibottom_dpsighost, length(psi_bdy_nodes_bed),length(psi_ghost)),1);

%% VORTICITY, eq. 24 with bcs 29, 31
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

%omega_top = -(h-h_prev)/(2*dx).*(du_dy_top_prev + du_dy_top);   
domegatop_du =  -(h-h_prev)/(2*dx).*ddudytop_du;
domegatop_duu_vec = -(h-h_prev)/(2*dx).*ddudytop_duup_vec;
domegatop_dud_vec = -(h-h_prev)/(2*dx).*ddudytop_dudp_vec;
domegatop_dh = sparse(psi_bdy_nodes_top, ones(length(psi_bdy_nodes_top),1), -1/(2*dx).*(du_dy_top_prev + du_dy_top), psi_nodes, 1); 
domegatop_dh_vec = -1/(2*dx).*(du_dy_top_prev + du_dy_top);

%domegadztop = (-9*omega(psi_bdy_nodes_top) + omega(psi_bdy_nodes_top + length(psi_bdy_nodes_top)) +8*omega_top)./(3*psi_Delta_z_cell(psi_bdy_nodes_top));
ddomegadztop_domega = sparse(psi_bdy_nodes_top, psi_bdy_nodes_top, -9.*psi_Delta_y_cell(psi_bdy_nodes_top)./(3*psi_Delta_z_cell(psi_bdy_nodes_top)), psi_nodes,psi_nodes) + ...
    sparse(psi_bdy_nodes_top, psi_bdy_nodes_top + length(psi_bdy_nodes_top), 1.*psi_Delta_y_cell(psi_bdy_nodes_top)./(3*psi_Delta_z_cell(psi_bdy_nodes_top)), psi_nodes,psi_nodes);
ddomegadztop_du = 8*sparse(1:length(psi_bdy_nodes_top),1:length(psi_bdy_nodes_top), 1.*psi_Delta_y_cell(psi_bdy_nodes_top)./(3*psi_Delta_z_cell(psi_bdy_nodes_top)), psi_nodes,psi_nodes)*domegatop_du;
ddomegadztop_dh = 8*sparse(1:length(psi_bdy_nodes_top),1:length(psi_bdy_nodes_top), 1.*psi_Delta_y_cell(psi_bdy_nodes_top)./(3*psi_Delta_z_cell(psi_bdy_nodes_top)), psi_nodes,psi_nodes)*domegatop_dh;

%net_omega_verflux(psi_bdy_nodes_top) = net_omega_verflux(psi_bdy_nodes_top) + h_prev.^(-1).*domegadztop;
domegaverflux_domega = domegaverflux_domega + h_prev.^(-1).*ddomegadztop_domega;
domegaverflux_du =  h_prev.^(-1).*ddomegadztop_du;
domegaverflux_dh =  h_prev.^(-1).*ddomegadztop_dh;

%bed: omega = omega_bed -((8 a - 9 f1 + f2)/(3 z))
%domegadz_bed =  (9*omega(psi_bdy_nodes_bed) - omega(psi_bdy_nodes_bed - length(psi_bdy_nodes_bed))-8*omega_bed)./(3*psi_Delta_z_cell(psi_bdy_nodes_bed));
domegadzbed_domegabed = -8./(3*psi_Delta_z_cell(psi_bdy_nodes_bed));

ddomegadzbed_domega = sparse(psi_bdy_nodes_bed,psi_bdy_nodes_bed, 9.*psi_Delta_y_cell(psi_bdy_nodes_bed)./(3*psi_Delta_z_cell(psi_bdy_nodes_bed)),psi_nodes,psi_nodes)+ ...
    sparse(psi_bdy_nodes_bed,psi_bdy_nodes_bed - length(psi_bdy_nodes_bed), -1.*psi_Delta_y_cell(psi_bdy_nodes_bed)./(3*psi_Delta_z_cell(psi_bdy_nodes_bed)),psi_nodes,psi_nodes);
ddomegadzbed_dpsi = sparse(1:length(psi_bdy_nodes_bed),1:length(psi_bdy_nodes_bed), domegadzbed_domegabed.*domegabed_dfluxpsibottom.*psi_Delta_y_cell(psi_bdy_nodes_bed),length(psi_bdy_nodes_bed),length(psi_bdy_nodes_bed))*dfluxpsibottom_dpsi;
ddomegadzbed_dpsighost = sparse(psi_bdy_nodes_bed,1:length(psi_ghost), domegadzbed_domegabed.*domegabed_dfluxpsibottom.*dfluxpsibottom_dpsighost.*psi_Delta_y_cell(psi_bdy_nodes_bed),psi_nodes,length(psi_ghost));
ddomegadzbed_dphi = sparse(psi_bdy_nodes_bed,index_ddphidybed_dphiup, domegadzbed_domegabed.*domegabed_ddphidybed.*ddphidybed_dphiup.*psi_Delta_y_cell(psi_bdy_nodes_bed),psi_nodes,T_nodes) + ...
    sparse(psi_bdy_nodes_bed,index_ddphidybed_dphidown, domegadzbed_domegabed.*domegabed_ddphidybed.*ddphidybed_dphidown.*psi_Delta_y_cell(psi_bdy_nodes_bed),psi_nodes,T_nodes);

%net_omega_verflux(psi_bdy_nodes_bed) = net_omega_verflux(psi_bdy_nodes_bed) - h_prev.^(-1).*domegadz_bed;
domegaverflux_domega = domegaverflux_domega - h_prev.^(-1).*ddomegadzbed_domega;
domegaverflux_dpsi = [sparse(psi_nodes-length(psi_bdy_nodes_bed), psi_nodes); - h_prev.^(-1).*ddomegadzbed_dpsi];
domegaverflux_dpsighost = - h_prev.^(-1).*ddomegadzbed_dpsighost;
domegaverflux_dphi = - h_prev.^(-1).*ddomegadzbed_dphi;

%conservation law
%div_omegafluxes =   net_omega_verflux./(h_prev.*psi_Delta_z_cell.*psi_Delta_y_cell)+ net_omega_horflux./psi_Delta_y_cell ;
fout(omega_index,omega_index) = spdiags(1./(h_prev.*psi_Delta_z_cell.*psi_Delta_y_cell),0,psi_nodes,psi_nodes)*domegaverflux_domega + ...
   spdiags(1./(psi_Delta_y_cell),0,psi_nodes,psi_nodes)*domegahorflux_domega;
fout(omega_index,psi_index) = spdiags(1./(h_prev.*psi_Delta_z_cell.*psi_Delta_y_cell),0,psi_nodes,psi_nodes)*domegaverflux_dpsi;
fout(omega_index,psighost_index) = spdiags(1./(h_prev.*psi_Delta_z_cell.*psi_Delta_y_cell),0,psi_nodes,psi_nodes)*domegaverflux_dpsighost;
fout(omega_index,phi_index) = spdiags(1./(h_prev.*psi_Delta_z_cell.*psi_Delta_y_cell),0,psi_nodes,psi_nodes)*domegaverflux_dphi;
fout(omega_index,u_index) = spdiags(1./(h_prev.*psi_Delta_z_cell.*psi_Delta_y_cell),0,psi_nodes,psi_nodes)*domegaverflux_du;

%% PHI, eq. 22 with bcs 26, 27

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
%fout(phi_index,h_index) = spdiags(1./(h_prev.*T_Delta_z_cell_volume.*T_Delta_y_cell),0,T_nodes,T_nodes)*dnetphiverflux_dh + h_prev.^(-1)*dSphi_dh;

fout(phi_index(1), :) = sparse(1, phi_index(1), 1, 1, length(v_in));
%% ALONG FLOW VELOCITY, U, eq. 19 with 20-21
%diffusive fluxes
%net_u_horflux = T_connect_hor * du_dy;
duhorflux_du = T_connect_hor * ddudy_du;

%net_u_verflux = T_connect_ver * (h_av.^(-1).*du_dz.*T_Delta_y_partial_cell);
duverflux_du = T_connect_ver * (spdiags(h_prev^(-1).*T_Delta_y_partial_cell,0,T_nedges_ver, T_nedges_ver)* ddudz_du);
duverflux_dh = T_connect_ver * spdiags(-h_prev^(-2).*T_Delta_y_partial_cell.*du_dz .*dhav_dh,0,T_nedges_ver, T_nedges_ver);

%eq. 21: prescribed flux at cell centre at the bed
% flux_u_bottom  = (gamma.*f_slide_Tbed.*u_bed).*T_Delta_y_cell(T_bdy_nodes_bed); 
% net_u_verflux(T_bdy_nodes_bed) = net_u_verflux(T_bdy_nodes_bed) - flux_u_bottom;

%flux_u_bottom  = (gamma.*(1+gamma_pert).*f_slide_Tbed.*u_bed).*T_Delta_y_cell(T_bdy_nodes_bed); 
dfluxubottom_du = sparse(T_bdy_nodes_bed, T_bdy_nodes_bed,gamma.*(1+gamma_pert).*f_slide_Tbed.*T_Delta_y_cell(T_bdy_nodes_bed), T_nodes, T_nodes);

%net_u_verflux(T_bdy_nodes_bed) = net_u_verflux(T_bdy_nodes_bed) - flux_u_bottom;
duverflux_du = duverflux_du - dfluxubottom_du;

%source term
%S_u = (h + b - h_prev - b_prev)/(dx).* ones(T_nodes, 1);
dSu_dh = 1./(dx).* ones(T_nodes, 1);

%conservation law
%fout(2*psi_nodes+T_nodes+1:2*psi_nodes+2*T_nodes) = div_u_fluxes - S_u;
%div_u_fluxes =   net_u_verflux./(h_av.*T_Delta_z_cell.*T_Delta_y_cell)+ net_u_horflux./T_Delta_y_cell; %nb: h_int depends on h, so must be differentiated!
fout(u_index,u_index) = spdiags(1./(h_av.*T_Delta_z_cell_volume.*T_Delta_y_cell),0,T_nodes,T_nodes)*duverflux_du + ...
   spdiags(1./(T_Delta_y_cell),0,T_nodes,T_nodes)*duhorflux_du;
%fout(u_index,h_index) = spdiags(1./(h_av.*T_Delta_z_cell_volume.*T_Delta_y_cell),0,T_nodes,T_nodes)*duverflux_dh - net_u_verflux./(h_av.^2.*T_Delta_z_cell_volume.*T_Delta_y_cell)*dhav_dh - dSu_dh;

%% P', eq. 25 with 32
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
% dfluxpsibottom_dpsibed_vec = 9./(3*psi_Delta_z_cell(psi_bdy_nodes_bed).*h_prev);
% dfluxpsibottom_dpsibed_vec_index = psi_bdy_nodes_bed;
% dfluxpsibottom_dpsiabed_vec =-1./(3*psi_Delta_z_cell(psi_bdy_nodes_bed).*h_prev);
% dfluxpsibottom_dpsiabed_vec_index = psi_bdy_nodes_bed-psi_nodes_hor(end);

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

%flux_p_bottom  = -domegady_bottom.*T_Delta_y_cell(T_bdy_nodes_bed); 
dfluxpbottom_dphi = -sparse(1:T_nodes,1:T_nodes,T_Delta_y_cell,T_nodes,T_nodes)*ddomegadybottom_dphi;
dfluxpbottom_dpsi = -sparse(1:T_nodes,1:T_nodes,T_Delta_y_cell,T_nodes,T_nodes)*ddomegadybottom_dpsi;
dfluxpbottom_dpsighost = -sparse(1:T_nodes,1:T_nodes,T_Delta_y_cell,T_nodes,T_nodes)*ddomegadybottom_dpsighost;

%net_p_verflux(T_bdy_nodes_bed) = net_p_verflux(T_bdy_nodes_bed) - flux_p_bottom;
dpverflux_dphi = -dfluxpbottom_dphi;
dpverflux_dpsi = -dfluxpbottom_dpsi;
dpverflux_dpsighost = -dfluxpbottom_dpsighost;

%dont forget dpverflux_du = dfluxptop_du; dpverflux_dh = dfluxptop_dh;
%conservation law 
%div_p_fluxes =   net_p_verflux./(h_prev.*T_Delta_z_cell.*T_Delta_y_cell)+ net_p_horflux./T_Delta_y_cell;
%fout(2*psi_nodes+2*T_nodes+1: 2*psi_nodes+3*T_nodes) = div_p_fluxes;
fout(p_index,p_index) = sparse(1:T_nodes, 1:T_nodes,1./(h_prev.*T_Delta_z_cell_volume.*T_Delta_y_cell),T_nodes,T_nodes)*dpverflux_dp+ sparse(1:T_nodes, 1:T_nodes,1./(T_Delta_y_cell),T_nodes,T_nodes)*dphorflux_dp;
fout(p_index,u_index) = sparse(1:T_nodes, 1:T_nodes,1./(h_prev.*T_Delta_z_cell_volume.*T_Delta_y_cell),T_nodes,T_nodes)*dpverflux_du;
%fout(p_index,h_index) = sparse(1:T_nodes, 1:T_nodes,1./(h_prev.*T_Delta_z_cell_volume.*T_Delta_y_cell),T_nodes,T_nodes)*dpverflux_dh;
fout(p_index,phi_index) = sparse(1:T_nodes, 1:T_nodes,1./(h_prev.*T_Delta_z_cell_volume.*T_Delta_y_cell),T_nodes,T_nodes)*dpverflux_dphi;
fout(p_index,psi_index) = sparse(1:T_nodes, 1:T_nodes,1./(h_prev.*T_Delta_z_cell_volume.*T_Delta_y_cell),T_nodes,T_nodes)*dpverflux_dpsi;
fout(p_index,psighost_index) = sparse(1:T_nodes, 1:T_nodes,1./(h_prev.*T_Delta_z_cell_volume.*T_Delta_y_cell),T_nodes,T_nodes)*dpverflux_dpsighost;

%NOTE: prescribe dirichlet condition on *one* surface node so as to fix p
%fout(2*psi_nodes+2*T_nodes+1) = p(1);
fout(p_index(1),:) = sparse(1,p_index(1),1,1,length(v_in));
































