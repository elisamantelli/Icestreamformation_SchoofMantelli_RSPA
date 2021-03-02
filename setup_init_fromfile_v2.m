function fout = setup_init_fromfile_v2(parameters,filename,index)

% setup_init_fromfile_v2.m provides the initial condition 

st = load(filename);

parameters_1D = st.parameters;
u = st.fout.v_in(st.index.u,index-1);
phi = st.fout.v_in(st.index.phi,index);
T = st.fout.v_in(st.index.T,index);

%caluate u_vert, du_dz_centre_full for this timestep
[~, faux_1D] = network_timestep_v5(st.fout.v_in(:,index-1),parameters_1D);
u_vert = faux_1D.u_vert;
dudzcentre_full = faux_1D.du_dz_centre_full;

%interpolate onto 2D grid: u_interp, phi_interp,u_vert_interp,
%dudzcentre_full_interp
Y_2D_uneven = parameters.grid.N.coor_nodes.y;
Z_2D_uneven = parameters.grid.N.coor_nodes.z;

%Y_1D = parameters_1D.grid.N.coor_nodes.y;
Z_1D = parameters_1D.grid.N.coor_nodes.z;
Z_1D_edges = parameters_1D.grid.N.coor_hor_cell_edges.z;
Y_hor = linspace(min(Y_2D_uneven),max(Y_2D_uneven),length(parameters.grid.N.bdy_nodes.bottom));
Y_hor_edges = linspace(min(parameters.grid.N.coor_hor_cell_edges.y),max(parameters.grid.N.coor_hor_cell_edges.y),length(parameters.grid.N.bdy_nodes.bottom));

Z_2D_even = repmat(Z_1D,[1,length(parameters.grid.N.bdy_nodes.bottom)]);
Z_2D_edges_even = [repmat(Z_1D_edges,[1,length(parameters.grid.N.bdy_nodes.bottom)]); zeros(1,length(parameters.grid.N.bdy_nodes.bottom))];
Y_2D_even = repmat(Y_hor,[size(Z_2D_even,1),1]);
Y_2D_edges_even = repmat(Y_hor_edges,[size(Z_2D_edges_even,1),1]);

u_array = repmat(u,[1,size(Z_2D_even,2)]);
phi_array = repmat(phi,[1,size(Z_2D_even,2)]);
uvert_array = repmat([u_vert; u(end)],[1,size(Z_2D_even,2)]);
dudzcentrefull_array = repmat(dudzcentre_full,[1,size(Z_2D_even,2)]);
T_array = repmat(T,[1,size(Z_2D_even,2)]);

u_interp = interp2(Y_2D_even,Z_2D_even,full(u_array),Y_2D_uneven,Z_2D_uneven);
T_interp = interp2(Y_2D_even,Z_2D_even,full(T_array),Y_2D_uneven,Z_2D_uneven);
phi_interp = interp2(Y_2D_even,Z_2D_even,full(phi_array),Y_2D_uneven,Z_2D_uneven);
u_vert_interp = interp2(Y_2D_edges_even,Z_2D_edges_even,full(uvert_array),parameters.grid.N.coor_hor_cell_edges.y,parameters.grid.N.coor_hor_cell_edges.z);
dudzcentre_full_interp = interp2(Y_2D_even,Z_2D_even,full(dudzcentrefull_array),Y_2D_uneven,Z_2D_uneven);

%% set up first iteration with perturbed bed temperature

if strcmp(parameters.flag.pert, 'rand') == 1
    rng('default') % For reproducibility
    T_pert = T_interp + parameters.amplitude*normrnd(0,abs(T(end)/2),[parameters.grid.N.n_nodes.tot, 1]);
    T_pert(parameters.grid.N.bdy_nodes.top) = parameters.T_s;
    parameters.T_bed =  T_pert(parameters.grid.N.bdy_nodes.bottom);
    parameters.T_bed(end) =  parameters.T_bed(1);
    parameters.gamma_pert = zeros(length(parameters.grid.N.coor_nodes.y(parameters.grid.N.bdy_nodes.bottom)), 1);
elseif strcmp(parameters.flag.pert, 'rand_gamma') == 1
    rng(1);
    T_pert = T_interp;
    parameters.gamma_pert = parameters.amplitude_gamma*normrnd(0,abs(parameters.gamma/3),[length(parameters.grid.N.coor_nodes.y(parameters.grid.N.bdy_nodes.bottom)), 1]);
    parameters.T_bed = T_pert(parameters.grid.N.bdy_nodes.bottom);
elseif strcmp(parameters.flag.pert, 'mono') == 1
    T_pert = T_interp;
    parameters.T_bed =  T(end) +parameters.amplitude*sin(2*pi.*parameters.n_wavelengths*parameters.grid.N.coor_nodes.y(parameters.grid.N.bdy_nodes.bottom)./(parameters.grid.N.extra.bd_y));
    parameters.gamma_pert = zeros(length(parameters.grid.N.coor_nodes.y(parameters.grid.N.bdy_nodes.bottom)), 1);
end

%set up the v_in vector for newton iteration of the mechanical problem at
%prescribed bed temperature

%set variables at previous time step
parameters.v_in_prev.h_prev = st.fout.v_in(st.index.h,index-1);
parameters.v_in_prev.u_prev = u_interp;
%parameters.v_in_prev.Q_prev = st.fout.v_in(st.index.Q,index)*parameters.grid.N.extra.bd_y/parameters_1D.grid.N.extra.bd_y; 
parameters.v_in_prev.h_pprev = st.fout.v_in(st.index.h,index-2);
parameters.v_in_prev.h_av_prev = 1/2*(parameters.v_in_prev.h_prev + parameters.v_in_prev.h_pprev);
parameters.v_in_prev.u_vert_prev = u_vert_interp; 
parameters.v_in_prev.du_dz_centre_full = dudzcentre_full_interp; 


parameters.h = st.fout.v_in(st.index.h,index);
parameters.Q = st.fout.v_in(st.index.Q,index)*parameters.grid.N.extra.bd_y/parameters_1D.grid.N.extra.bd_y; 
u_in = u_interp;
p_in = zeros(parameters.grid.N.n_nodes.tot,1); 
phi_in = phi_interp;
omega_in = zeros(parameters.grid.psi.n_nodes.tot,1);
psi_in = zeros(parameters.grid.psi.n_nodes.tot,1);
psi_ghost = 0;

v_in_pert_mech = [psi_in; omega_in; phi_in; u_in; p_in; psi_ghost];

srch.tolF = 2e-03;
srch.verbose = 1;
srch.itmax = 10;

[fout_Tfixed,~,faux_Tfixed] = Newton_v2(@network_timestep_v5_mech_Tfixed_v2,@network_timestep_v5_mech_Tfixed_jacobian_v2,v_in_pert_mech,parameters,srch);

% % PLOTTING
% psi = fout_Tfixed(1:parameters.grid.psi.n_nodes.tot);
% phi = fout_Tfixed(2*parameters.grid.psi.n_nodes.tot+1:2*parameters.grid.psi.n_nodes.tot+parameters.grid.N.n_nodes.tot);
% w = faux_Tfixed.w;
% v = faux_Tfixed.v;
% dpsi_dy = faux_Tfixed.dpsi_dy;
% 
% y = parameters.grid.psi.coor_ver_cell_edges.y;
% z = parameters.grid.psi.coor_ver_cell_edges.z;
% [dpsi_dy_even,ygrid,zgrid] =gridfit(y,z,dpsi_dy,200,200);
% figure; subplot(144)
% contourf(ygrid,zgrid,dpsi_dy_even,200);
% xlabel('y');zlabel('z'); title('dpsi/dy');
% 
% y = parameters.grid.N.coor_nodes.y;
% z = parameters.grid.N.coor_nodes.z;
% [psi_even,ygrid,zgrid] =gridfit(y,z,phi,200,200);
%  subplot(141)
% contourf(ygrid,zgrid,psi_even,200);
% xlabel('y');zlabel('z'); title('phi');
% 
% y = parameters.grid.psi.coor_hor_cell_edges.y;
% z = parameters.grid.psi.coor_hor_cell_edges.z;
% [phi_even,ygrid,zgrid] =gridfit(y,z,v,200,200);
% subplot(142)
% contourf(ygrid,zgrid,phi_even,200);
% xlabel('y');zlabel('z'); title('dpsi/dz on T grid');
% 
% y = parameters.grid.N.coor_hor_cell_edges.y;
% z = parameters.grid.N.coor_hor_cell_edges.z;
% [w_even,ygrid,zgrid] =gridfit(y,z,w,200,200);
% subplot(143)
% contourf(ygrid,zgrid,w_even,200);
% xlabel('y');zlabel('z'); title('-dpsi/dy  on T grid');

%% construct temperature field in the bed

%bed
[Y_bed,Z_bed] = meshgrid(linspace(min(parameters.grid.Tb.coor_nodes.y),max(parameters.grid.Tb.coor_nodes.y),length(parameters.grid.Tb.bdy_nodes.top)), ...
    linspace(min(parameters.grid.Tb.coor_nodes.z),max(parameters.grid.Tb.coor_nodes.z),150));
T_bed_array_init = repmat(parameters.T_bed.',[size(Z_bed,1),1]);
if strcmp(parameters.flag.pert, 'rand_gamma') == 1
    T_bed_array = T_bed_array_init - parameters.nu.*Z_bed;
else
    T_bed_array = T_bed_array_init.*exp(Z_bed/0.1) - parameters.nu.*Z_bed;
end

Yb = parameters.grid.Tb.coor_nodes.y;
Zb = parameters.grid.Tb.coor_nodes.z;
if parameters.flag1d == 0
    Tb_pert = interp2(Y_bed,Z_bed,T_bed_array,Yb,Zb);
elseif parameters.flag1d == 1
    Tb_pert = interp2(Y_bed,Z_bed,T_bed_array,Yb,Zb);
    %Tb_pert = interp1(Z_bed,T_bed_array,Zb);
end

 
%% input to timestepping for the fully coupled problem

T_nodes = parameters.grid.N.n_nodes.tot;                                   %number of nodes 
psi_nodes = parameters.grid.psi.n_nodes.tot; 

%set variables at previous time step
parameters.v_in_prev.h_prev = parameters.h;
parameters.v_in_prev.u_prev = fout_Tfixed(2*psi_nodes+T_nodes+1:2*psi_nodes+2*T_nodes);
parameters.v_in_prev.Q_prev = st.fout.v_in(st.index.Q,index)*parameters.grid.N.extra.bd_y/parameters_1D.grid.N.extra.bd_y;
%b_prev = parameters.v_in_prev.b_prev;
parameters.v_in_prev.qx_prev = sparse(length(parameters.grid.N.bdy_nodes.bottom),1);
parameters.v_in_prev.T_prev = T_pert; 
%parameters.v_in_prev.T_prev(parameters.grid.N.bdy_nodes.bottom) = max(faux.T_bed) + 0.1*sin(2*pi.*parameters.grid.N.coor_nodes.y(parameters.grid.N.bdy_nodes.bottom)./parameters.grid.N.extra.bd_y);
parameters.v_in_prev.h_pprev = st.fout.v_in(st.index.h,index-1);
parameters.v_in_prev.h_av_prev = 1/2*(parameters.v_in_prev.h_prev + parameters.v_in_prev.h_pprev);
parameters.v_in_prev.du_dz_centre_full_prev = faux_Tfixed.du_dz_centre_full;
parameters.v_in_prev.du_dy_centre_prev = faux_Tfixed.du_dy_centre;
parameters.v_in_prev.u_vert_prev = faux_Tfixed.u_vert;
parameters.v_in_prev.I_prev = sparse(length(parameters.grid.N.bdy_nodes.bottom),1);
parameters.v_in_prev.f_prev = faux_Tfixed.f;
%dx_prev =  parameters.timestep.dx_prev;

%set up initial condition/guess
h_in = parameters.v_in_prev.h_prev;
Q_in = parameters.v_in_prev.Q_prev;
u_in = fout_Tfixed(2*psi_nodes+T_nodes+1:2*psi_nodes+2*T_nodes);
p_in = fout_Tfixed(2*psi_nodes+2*T_nodes+1: 2*psi_nodes+3*T_nodes);
T_in = T_pert; 
Tb_in = Tb_pert;
phi_in = fout_Tfixed(2*psi_nodes+1:2*psi_nodes+T_nodes);
omega_in = fout_Tfixed(psi_nodes+1: 2*psi_nodes);
psi_in = fout_Tfixed(1:psi_nodes);
psi_ghost = fout_Tfixed(2*psi_nodes+3*T_nodes+1);
Pi_in = zeros(length(parameters.grid.N.bdy_nodes.bottom),1);

v_in = [psi_in; omega_in; phi_in; u_in; p_in; T_in; Tb_in; psi_ghost; h_in; Q_in; Pi_in];

fout.parameters = parameters;
fout.v_in = v_in;