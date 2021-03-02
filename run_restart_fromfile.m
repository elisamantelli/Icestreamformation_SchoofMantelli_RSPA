N = 8;
maxNumCompThreads(N)
filename = 'divide2D_friction_coulomb_delta_0.03_hd_2_alpha_1_gamma_3_nu_0.5_Pe_1_dx_0.005.mat';
st = load(filename);
j = 180;

parameters = st.parameters; 
parameters.filename = filename;

%set perturbation parameters
parameters.flag.pert = 'rand_gamma';
% parameters.h_init = 2; 
% parameters.h = parameters.h_init;
parameters.iter.x_current = parameters.timestep.x_init;
parameters.timestep.dx = 0.005; %(i+1)
parameters.flag_Tdep = 1; %turn on (1) and off(0) T_dep sliding
parameters.amplitude_gamma = parameters.reg.epsilon_f/100;    %amplitude of bed temperature perturbation; set to zero to compute steady states
parameters.n_wavelengths = 1; %number of wavelengths within assigneddomain size

%regularization parameters 

parameters.timestep.x_init = st.fout.x(j);
parameters.timestep.x = parameters.timestep.x_init; % update at each timestep
parameters.tplot = parameters.timestep.dx;
parameters.stepmax = 5000;


%grid generation
parameters.n_nodes_transverse = 100;             %#nodes in the transverse direction in the uppermost layer. Remaining grid parameters are set in the grid routine. Must be even!
parameters.ratio_hor = 7;  %dy = dz*ratio
parameters.ratio_vert = 2;
parameters.flag1d = 0;
if parameters.flag1d == 0
    parameters.grid = fv_grid_transverse_v6(parameters);
elseif parameters.flag1d == 1
    parameters.grid = fv_grid_transverse_1d_v4(parameters);
end



fout_init = setup_init_fromfile_v2(parameters,filename,j);
v_in = fout_init.v_in;
parameters = fout_init.parameters;

timestepping_cluster(parameters, v_in)

quit;