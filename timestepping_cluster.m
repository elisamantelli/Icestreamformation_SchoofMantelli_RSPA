
function timestepping_cluster(parameters, v_in, varargin)

if nargin == 3
    save_dir = char(varargin);
    current_dir = cd;
else
    save_dir = cd;
    current_dir = cd;
end

parameters.filename = ['4_divide3D_friction_' parameters.flag.hydrology '_latdrainage_' parameters.flag.fluxy_sigma '_pertampl_' num2str(parameters.amplitude) '_pertmode_' parameters.flag.pert '_width_' num2str(parameters.grid.N.extra.bd_y)...
    '_hd_' num2str(parameters.h_init) '_alpha_' num2str(parameters.alpha) '_gamma_' num2str(parameters.gamma) '_nu_' num2str(parameters.nu) '_Pe_' num2str(parameters.Pe) '_dx_' num2str(parameters.timestep.dx) '_beta_' num2str(parameters.beta) ...
    '_Pi0_' num2str(parameters.reg.Pi_0) '_Pis_' num2str(parameters.reg.Pi_s) '_mu_' num2str(parameters.mu) ];
% parameters.filename = ['divide2D_friction_' parameters.flag.hydrology '_delta_' num2str(parameters.reg.epsilon_f) '_hd_' num2str(parameters.h_init) '_alpha_' ...
%     num2str(parameters.alpha) '_gamma_' num2str(parameters.gamma) '_nu_' num2str(parameters.nu) '_Pe_' num2str(parameters.Pe) '_dx_' num2str(parameters.timestep.dx) '_slope_' num2str(parameters.bed.b1) '_ratiohor_' num2str(parameters.ratio_hor)];
%

%indexing of v_in
T_nodes = parameters.grid.N.n_nodes.tot;   
Tb_nodes = parameters.grid.Tb.n_nodes.tot;
psi_nodes = parameters.grid.psi.n_nodes.tot;
T_bdy_nodes_bed = parameters.grid.N.bdy_nodes.bottom;

index.psi = 1:psi_nodes;
index.omega = psi_nodes+1: 2*psi_nodes;
index.phi = 2*psi_nodes+1:2*psi_nodes+T_nodes;
index.u = 2*psi_nodes+T_nodes+1:2*psi_nodes+2*T_nodes;
index.p = 2*psi_nodes+2*T_nodes+1: 2*psi_nodes+3*T_nodes;
index.T = 2*psi_nodes+3*T_nodes+1: 2*psi_nodes+4*T_nodes;
index.Tb = 2*psi_nodes+4*T_nodes+1: 2*psi_nodes+4*T_nodes+Tb_nodes;
index.h = 2*psi_nodes+4*T_nodes+Tb_nodes+2;
index.Q = 2*psi_nodes+4*T_nodes+Tb_nodes+3;
index.psighost = 2*psi_nodes+4*T_nodes+Tb_nodes+1;
index.Pi =2*psi_nodes+4*T_nodes+Tb_nodes+4: 2*psi_nodes+4*T_nodes+Tb_nodes+3+length(T_bdy_nodes_bed);

%parameters for Newton iteration
srch.tolF= 1e-02;
srch.verbose= 1;
srch.itmax = 15;
srch.toldelta = 1e-05;

step = 1;

%pre-allocation and saving the first step
f_v_in = sparse(length(v_in), parameters.stepmax);
f_a = sparse(1,parameters.stepmax);
f_x = sparse(1,parameters.stepmax);
f_gamma = sparse(length(T_bdy_nodes_bed),parameters.stepmax);

f_v_in(:,1) = v_in;
f_a(1) = parameters.a;
f_x(1) = parameters.timestep.x;
if isfield(parameters,'gamma_pert')==1
    f_gamma(:,1) = parameters.gamma_pert;
elseif isfield(parameters,'gamma_pert')==0
    f_gamma(:,1) = f_gamma(:,1);
end

while step <= parameters.stepmax-1 %&& parameters.v_in_prev.h_prev >=0
    
    parameters.timestep.x = parameters.timestep.x + 1/2*parameters.timestep.dx + 1/2*parameters.timestep.dx_prev;
    display(['x = ' num2str(parameters.timestep.x)])
    
    v_in_input = v_in;
    
    %preconditioning
    v_in_input(index.h) = parameters.v_in_prev.h_prev + (parameters.v_in_prev.h_prev - parameters.v_in_prev.h_pprev)*(parameters.timestep.dx + parameters.timestep.dx_prev)/parameters.timestep.dx_prev;
    
    %solve
    tic
    [vout,error_flag, faux] = Newton_v2(@network_timestep_v5,@network_timestep_v5_jacobian,v_in_input,parameters,srch);
    toc
    
    %update timestep
    
    parameters.timestep.dx_prev = parameters.timestep.dx;
    parameters.timestep.dx_pprev = parameters.timestep.dx_prev;
    
    parameters.v_in_prev.h_prev = vout(index.h);
    parameters.v_in_prev.u_prev = vout(index.u);
    parameters.v_in_prev.Q_prev = vout(index.Q);
    parameters.v_in_prev.qx_prev = faux.qx;
    parameters.v_in_prev.T_prev = vout(index.T);
    parameters.v_in_prev.h_pprev = faux.h_prev;
    parameters.v_in_prev.h_av_prev = faux.h_av;
    parameters.v_in_prev.du_dz_centre_full_prev = faux.du_dz_centre_full;
    parameters.v_in_prev.du_dy_centre_prev = faux.du_dy_centre;
    parameters.v_in_prev.u_vert_prev = faux.u_vert;
    parameters.v_in_prev.I_prev = faux.I;
    parameters.v_in_prev.dsigma_dy_prev = faux.dsigma_dy;
    parameters.v_in_prev.f_prev= faux.f;
    
    v_in = vout;
    f_v_in(:,step+1) = vout;
    f_a(step+1) = parameters.a;
    f_x(step+1) = parameters.timestep.x;
    if isfield(parameters,'gamma_pert')==1
        f_gamma(:,step+1) = parameters.gamma_pert;
    elseif isfield(parameters,'gamma_pert')==0
        f_gamma(:,step+1) = f_gamma(:,step+1);
    end
    
    
    
    %update parameters
    parameters.a = parameters.a + parameters.da;
    %perturbed friction coeff
    if strcmp(parameters.flag.pert, 'mono') == 1
        parameters.gamma_pert = ones(length(parameters.grid.N.coor_nodes.y(parameters.grid.N.bdy_nodes.bottom)), 1);
    elseif strcmp(parameters.flag.pert, 'rand') == 1
        parameters.gamma_pert = zeros(length(parameters.grid.N.coor_nodes.y(parameters.grid.N.bdy_nodes.bottom)), 1);
    elseif strcmp(parameters.flag.pert, 'rand_gamma') == 1
        parameters.gamma_pert = parameters.amplitude_gamma*normrnd(0,abs(parameters.gamma/3),[length(parameters.grid.N.coor_nodes.y(parameters.grid.N.bdy_nodes.bottom)), 1]);
    else
        parameters.gamma_pert = zeros(length(parameters.grid.N.coor_nodes.y(parameters.grid.N.bdy_nodes.bottom)), 1);
    end
    
    
    %saving
    
    
    save_step = step/30;
    
    
    if floor(save_step) == save_step 
        
        fout.a = f_a(1:step+1);
        fout.x = f_x(1:step+1);
        fout.v_in = f_v_in(:,1:step+1);
        fout.gamma = f_gamma(:,1:step+1);
        
        cd(save_dir)
        save([parameters.filename '.mat'],'fout','parameters', 'index', '-v7.3')
        cd(current_dir)
        
    elseif error_flag == 1
        fout.a = f_a(1:step);
        fout.x = f_x(1:step);
        fout.v_in = f_v_in(:,1:step);
        fout.gamma = f_gamma(:,1:step);
        
        cd(save_dir)
        save([parameters.filename '.mat'],'fout','parameters', 'index', '-v7.3')
        cd(current_dir)
    end
    
    if error_flag
        warning('Iteration failure: convergence not achieved')
        quit
    end
    
    step = step+1;
    display(['step # ' num2str(step-1) ' of ' num2str(parameters.stepmax)])
    
end