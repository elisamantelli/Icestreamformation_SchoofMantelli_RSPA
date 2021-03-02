 function fout = fv_grid_transverse_v6(parameters)
 
%2-dimensional (z and y) grid for three dimensional calculations. 
% Implements layered grid structure in the vertical, uses periodic boundaries laterally,
% and is set up for thrid order polynomial interpolation at hanging nodes 

% Elisa Mantelli, Dec 2020
 
nodesy = parameters.n_nodes_transverse;
ratio_ver = parameters.ratio_vert;
ratio_hor = parameters.ratio_hor;

%Psi grid
extra.psi.n_layers = 4; %2
extra.psi.bd_z.ice = 1;
extra.psi.bd_z.bed = -3;

extra.psi.layer_height = [0.5; 0.25; 0.125; 0.125];%[0.5; 0.5];%

extra.psi.n_v = ratio_ver*[4; 4; 4; 8];%2*[10; 10; 10; 15];%[4; 4; 4; 6];%[10; 20];%

dz = extra.psi.layer_height./extra.psi.n_v;
dy = [ratio_hor*dz(1:3); ratio_hor*0.5*dz(3)];                             %THIS CHANGES DEPENDING ON THE NUMBER OF LAYERS
extra.psi.bd_y = dy(1)*nodesy;                                             %THIS MUST BE AN EVEN NUMBER
extra.psi.n_y = extra.psi.bd_y./dy;

fout.psi = grid_psi(extra.psi,dy,dz);
fout.psi.extra = extra.psi;

%N, T_ice grid
extra.N.n_layers = extra.psi.n_layers; 
extra.N.bd_z.ice = extra.psi.bd_z.ice;
extra.N.bd_z.bed = extra.psi.bd_z.bed;

extra.N.layer_height = extra.psi.layer_height;
extra.N.n_v = extra.psi.n_v;

Ndz = extra.N.layer_height./extra.N.n_v;
extra.N.n_v(end) = extra.N.n_v(end) + 1;                                   % add bed (z=0) as a node to the lowermost layer in T grid
Ndy = dy;                                                                  %THIS CHANGES DEPENDING ON THE NUMBER OF LAYERS
extra.N.bd_y = extra.psi.bd_y;
extra.N.n_y = extra.N.bd_y./Ndy;       

fout.N = grid_T(extra.N, Ndy,Ndz,1);
fout.N.extra = extra.N;

%T_bed grid
extra.Tb.n_layers = extra.psi.n_layers + 1;                                %note that this accounts both for the layer across the bed and the extra layer for z<-1
extra.Tb.bd_z.ice = 0-Ndz(end);                                            %this is the coordinate of the top node
extra.Tb.bd_z.bed = extra.psi.bd_z.bed;

extra.Tb.layer_height = [flip(extra.psi.layer_height(1:end)); abs(extra.Tb.bd_z.bed)- extra.psi.bd_z.ice];
extra.Tb.layer_height(1) = extra.Tb.layer_height(1)-Ndz(end);
extra.Tb.n_v = [flip(extra.psi.n_v(1:end));floor(extra.Tb.layer_height(end)/(Ndz(1)))];
extra.Tb.n_v(1) = extra.Tb.n_v(1)-1;

Tbdz = extra.Tb.layer_height./extra.Tb.n_v;
%extra.Tb.n_v(end) = extra.Tb.n_v(end) + 1;                                 % add bed as a node to the lowermost layer
Tbdy = [flip([2*(Ndy(1)); Ndy])]; %[1.5*Tdz(1:2); 1.5*Tdz(1); 2*(1.5*Tdz(1))];% %this is the size of cells %THIS CHANGES DEPENDING ON THE NUMBER OF LAYERS
extra.Tb.bd_y = extra.psi.bd_y;
extra.Tb.n_y = extra.Tb.bd_y./Tbdy;                                        % this accounts for the fact that T nodes are on the corners of psi cells, and therefore on domain boundaries

fout.Tb = grid_T(extra.Tb, Tbdy,Tbdz,0);
fout.Tb.extra = extra.Tb;
 end

function fout = grid_psi(grid, dy,dz)
% save/define key grid parameters
fout.n_nodes.vert = grid.n_v;                                              %list n_layer-by-one
fout.n_nodes.hor = grid.n_y;                                               %list n_layer-by-one
fout.n_nodes.tot = sum(fout.n_nodes.hor.*fout.n_nodes.vert);
fout.n_nodes.tot_virtual = sum(fout.n_nodes.hor.*fout.n_nodes.vert) + sum(2*fout.n_nodes.hor(2:end));
fout.n_edges.hor = sum((fout.n_nodes.hor).*fout.n_nodes.vert);
fout.n_edges.vert = 0;                                                     %build along with the grid

%list of nodes where boundary conditions apply
fout.bdy_nodes.top = (1:fout.n_nodes.hor(1))';
fout.bdy_nodes.bed = (fout.n_nodes.tot - fout.n_nodes.hor(end)+1:fout.n_nodes.tot)';

fout.bdy_nodes.inflow = [];                                                %build along with the grid
fout.bdy_nodes.outflow = [];                                               %build along with the grid

%build the grid a layer at a time, starting from the top

%initialization
fout.coor_nodes.y = [];
fout.coor_nodes.z = [];
fout.Delta_y_cell = [];
fout.Delta_z_cell = [];
fout.delta_y_vflux = [];
fout.coor_hor_cell_edges.z = [];
fout.coor_hor_cell_edges.y = [];
fout.coor_ver_cell_edges.z = [];
fout.coor_ver_cell_edges.y = [];
fout.up_node.hor = [];
fout.up_node.vert = [];
fout.down_node.hor = [];
fout.down_node.vert = [];

%these now have length of #vertical network edges in the real network. Must
%have instead length = number of virtual nodes; must also contain indices
%of real network, so the only thing to update possibly is the ordering
fout.coor_node_virtual.y = [];  %y coordinate of the virtual node - this is where the interpolation is performed -OK
fout.coor_node_virtual.z = [];

%lateral interpolation
fout.coor_node_h_y.up = [];   %coordinate of the real node horizontally upstream of current virtual node
fout.coor_node_h_y.down = []; %coordinate of the real node horizontally downstream of current virtual node
fout.index_node_hdown = [];   %index of the REAL node downstream (for regular -non-interface- nodes) this will be the index of the real node co-located with the virtual one
fout.index_node_hup = [];
fout.coor_node_h_y.upup = [];
fout.coor_node_h_y.downdown = [];
fout.index_node_hdowndown = [];
fout.index_node_hupup = [];

%vertical interpolation - ADD COORDINATES!!
fout.down_node_virtual.vert = [];
fout.up_node_virtual.vert = [];
fout.downdown_node_virtual.vert = [];
fout.upup_node_virtual.vert = [];

fout.coor_down_node_virtual_vert = [];
fout.coor_up_node_virtual_vert = [];
fout.coor_downdown_node_virtual_vert = [];
fout.coor_upup_node_virtual_vert = [];


for ss = 1:grid.n_layers
    
    index_la = sum(grid.n_y(1:ss-1).*grid.n_v(1:ss-1));
    if ss == 1
        index_la_virtual = index_la;
    elseif ss>1
        index_la_virtual = sum(grid.n_y(1:ss-1).*grid.n_v(1:ss-1)) + sum(4*grid.n_y(1:ss-1));
    end
    n_nodes_layer = grid.n_y(ss)*grid.n_v(ss);
    n_edges_layer_ver = grid.n_y(ss)*(grid.n_v(ss)-1);
    n_nodes_layer_virtual = n_nodes_layer;
    index_nodes = index_la+(1:n_nodes_layer).';
    index_nodes_virtual = index_la_virtual+(1:n_nodes_layer).';
    index_nodes_layer = (1:n_nodes_layer).';
    index_nodes_layer_virtual = (1:n_nodes_layer_virtual).';
    %inflow and outflow boundary
    fout.bdy_nodes.inflow = [fout.bdy_nodes.inflow; index_la + (1:grid.n_y(ss):n_nodes_layer).'];
    fout.bdy_nodes.outflow = [fout.bdy_nodes.outflow; index_la+(grid.n_y(ss):grid.n_y(ss):n_nodes_layer).']; 
    
    %cell centres locations (list n_nodes_layer-by-one)
    z_top = 1-sum(grid.layer_height(1:ss-1));
    z_bottom = 1-sum(grid.layer_height(1:ss));
    
    z_nodes_line = (z_top-dz(ss)/2:-dz(ss):z_bottom).';
    coor_nodes.z = reshape(repmat(z_nodes_line,[1,grid.n_y(ss)]).', n_nodes_layer, []);
    
    y_nodes_line = (dy(ss)/2:dy(ss):grid.bd_y).';
    coor_nodes.y = reshape(repmat(y_nodes_line,[grid.n_v(ss),1]).', n_nodes_layer, []);
    coor_node_virtual.y = coor_nodes.y;
    coor_node_virtual.z = coor_nodes.z;
    
    %cell size (n_nodes_layer-by-one)
    Delta_y_cell = dy(ss)*ones(n_nodes_layer,1);
    Delta_z_cell = dz(ss)*ones(n_nodes_layer,1);
    
    % length of horizontal cell edges for computation of vertical fluxes (only interior vertical edges)
    delta_y_vflux = dy(ss)*ones((grid.n_v(ss)-1)*grid.n_y(ss),1);
    
    %coordinate of horizontal cell edges, including only the edges internal
    %to the layer. Edges on the bottom boundary are
    %added in the for loop below. 
    %NOTE: node numbering is from upper left corner to bottom
    %right corner along rows; same for edges) 
    
    n_edges_int = (grid.n_v(ss)-1)*grid.n_y(ss);
    coor_hor_cell_edges.z = (z_top-dz(ss):-dz(ss):z_bottom+dz(ss)).';
    coor_hor_cell_edges.z = reshape(repmat(coor_hor_cell_edges.z,[1,grid.n_y(ss)]).', n_edges_int,[]);
   
    coor_hor_cell_edges.y = y_nodes_line;
    coor_hor_cell_edges.y = reshape(repmat(coor_hor_cell_edges.y.',[grid.n_v(ss)-1,1]).',n_edges_int,[]);
    
    coor_ver_cell_edges.y = y_nodes_line +dy(ss)/2;
    coor_ver_cell_edges.y = reshape(repmat(coor_ver_cell_edges.y.',[grid.n_v(ss),1]).',n_nodes_layer,[]);
    
    coor_ver_cell_edges.z = coor_nodes.z;
    
    %upnodes/downnodes of network edges inside the layer, excluding top, bottom
    %boundaries. Assumes a PERIODIC GRID IN Y
    %NOTE: network edges numbering is from upper left corner to bottom
    %right corner along rows, both for horizontal and vertical edges
    
    %hor network edges (up_node and down_node list up node and down node
    %for each edge; numbering along rows, edge #1 connects node#1 to node#2)

    up = index_nodes_layer;
    down = reshape((circshift(reshape(index_nodes_layer, [grid.n_y(ss) grid.n_v(ss)]).',-1,2)).',n_nodes_layer,1);
    up_node.hor = index_nodes(up); 
    down_node.hor = index_nodes(down);
    
    %ver network edges (excluding edges connecting adjacent layers)
    up = index_nodes_layer; up(1:grid.n_y(ss))=[];
    down = index_nodes_layer; down(n_nodes_layer-grid.n_y(ss)+1:n_nodes_layer)=[];
    up_node.ver = index_nodes(up); 
    down_node.ver = index_nodes(down); 
    
    up = index_nodes_layer_virtual; up(1:grid.n_y(ss))=[];
    down = index_nodes_layer_virtual; down(n_nodes_layer_virtual-grid.n_y(ss)+1:n_nodes_layer_virtual)=[];
    
    up_node_virtual.ver = index_nodes_virtual(up); 
    down_node_virtual.ver = index_nodes_virtual(down); 
    upup_node_virtual.ver = sparse(length(up),1);
    downdown_node_virtual.ver = sparse(length(up),1);
    coor_down_node_virtual_vert = sparse(length(up),1);
    coor_up_node_virtual_vert = sparse(length(up),1);
    coor_downdown_node_virtual_vert = sparse(length(up),1);
    coor_upup_node_virtual_vert = sparse(length(up),1);
    
    %virtual network
    
    %this is the actual down_node, so length must be n_nodes_layer_virtual.
    %Also, these contain the index of the up/down/.. nodes indexed in terms
    %of the actual network
    index_node_hdown = index_nodes;
    coor_node_h_y.down = coor_node_virtual.y;
    
    %one further downstream
    index_node_hdowndown = reshape((circshift(reshape(index_node_hdown, [grid.n_y(ss) grid.n_v(ss)]).',-1,2)).',n_nodes_layer_virtual,1);
    coor_intermediate = (circshift(reshape(coor_node_h_y.down, [grid.n_y(ss) grid.n_v(ss)]).',-1,2));
    coor_intermediate(:,end) = coor_intermediate(:,end-1)+ dy(ss);
    coor_node_h_y.downdown = reshape(coor_intermediate.',n_nodes_layer_virtual,1);
    
    %one upstream
    index_node_hup = reshape((circshift(reshape(index_node_hdown, [grid.n_y(ss) grid.n_v(ss)]).',1,2)).',n_nodes_layer_virtual,1);
    coor_intermediate = (circshift(reshape(coor_node_h_y.down, [grid.n_y(ss) grid.n_v(ss)]).',1,2));
    coor_intermediate(:,1) = coor_intermediate(:,2)- dy(ss);
    coor_node_h_y.up = reshape(coor_intermediate.',n_nodes_layer_virtual,1);
    
    %two upstream
    index_node_hupup = reshape((circshift(reshape(index_node_hdown, [grid.n_y(ss) grid.n_v(ss)]).',2,2)).',n_nodes_layer_virtual,1);
    coor_intermediate = (circshift(reshape(coor_node_h_y.down, [grid.n_y(ss) grid.n_v(ss)]).',2,2));
    coor_intermediate(:,2) = coor_intermediate(:,3)- dy(ss);
    coor_intermediate(:,1) = coor_intermediate(:,2)- dy(ss);
    coor_node_h_y.upup = reshape(coor_intermediate.',n_nodes_layer_virtual,1);
    
    %add ver network edges that connect lower to upper layer
    if ss>1
        %down nodes
        %vertical, actual network
        down_node_interface = reshape(repmat((index_la-grid.n_y(ss-1)+1:index_la),[2 1]),[ ],1);
        down_node.ver = [down_node_interface;down_node.ver];
        
        %edges of the vertical virtual network
        downdown_node_virtual_interface = (index_la_virtual-2*grid.n_y(ss)+1:index_la_virtual-grid.n_y(ss)).';
        downdown_node_virtual.ver = [downdown_node_virtual_interface; downdown_node_virtual.ver]; 
        
        down_node_virtual_interface = (index_la_virtual-grid.n_y(ss)+1:index_la_virtual).';
        down_node_virtual.ver = [down_node_virtual_interface; down_node_virtual.ver]; 
        
        % up nodes (must be in the same number as down_node_interface)
        up_node_interface = (index_la+1 : index_la+grid.n_y(ss)).';
        up_node.ver = [up_node_interface; up_node.ver];
        
        up_node_virtual_interface = (index_la_virtual+1:index_la_virtual+grid.n_y(ss)).';
        up_node_virtual.ver = [up_node_virtual_interface; up_node_virtual.ver]; 
        
        upup_node_virtual_interface = (index_la_virtual+grid.n_y(ss)+1:index_la_virtual+2*grid.n_y(ss)).';
        upup_node_virtual.ver = [upup_node_virtual_interface; upup_node_virtual.ver]; 
        
        %add location (x,z) of hor cell edges and virtual node networks
        coor_hor_cell_edges.y = [y_nodes_line ; coor_hor_cell_edges.y];
        coor_hor_cell_edges.z = [z_top*ones(grid.n_y(ss),1);coor_hor_cell_edges.z];
        
        coor_node_virtual.y = [y_nodes_line ; y_nodes_line ; coor_node_virtual.y];
        coor_node_virtual.z = [z_top + 3/2*dz(ss-1)*ones(grid.n_y(ss),1); z_top+ dz(ss-1)/2*ones(grid.n_y(ss),1); coor_node_virtual.z];
        coor_down_node_virtual_vert = [ z_top+ dz(ss-1)/2*ones(grid.n_y(ss),1) ;coor_down_node_virtual_vert];
        coor_up_node_virtual_vert = [z_top - dz(ss)/2*ones(grid.n_y(ss),1);coor_up_node_virtual_vert];
        coor_downdown_node_virtual_vert = [z_top + 3/2*dz(ss-1)*ones(grid.n_y(ss),1) ;coor_downdown_node_virtual_vert];
        coor_upup_node_virtual_vert = [z_top - 3*dz(ss)/2*ones(grid.n_y(ss),1); coor_upup_node_virtual_vert];
       
        %add surface crossed by vertical network edge
        delta_y_vflux = [Delta_y_cell(1:grid.n_y(ss)); delta_y_vflux];
        
        %define nodes and nodes coordinates for connected by the hor
        %network edges intersected by the virtual verical netweo edge.
        %Provides nodes for a third order polynomial interpolation
        %laterally - LATERAL INTERPOLATION
        
        %one node downstream
        index_node_hdown_interface_1 = reshape(repmat((index_la-grid.n_y(ss-1)+2:index_la),[2 1]),[ ],1);
        index_node_hdown_interface_2 = reshape(repmat((index_la-2*grid.n_y(ss-1)+2:index_la - grid.n_y(ss-1)),[2 1]),[ ],1);
        index_node_hdown = [index_la-2*grid.n_y(ss-1)+1; index_node_hdown_interface_2; index_la-2*grid.n_y(ss-1)+1;...
            index_la-grid.n_y(ss-1)+1; index_node_hdown_interface_1; index_la-grid.n_y(ss-1)+1;...
            index_node_hdown];
        coor_node_h_y_down_interface = reshape(repmat((fout.coor_nodes.y(index_la-grid.n_y(ss-1)+2:index_la)).',[2 1]),[ ],1);
        coor_node_h_y.down = [fout.coor_nodes.y(index_la-grid.n_y(ss-1)+1); coor_node_h_y_down_interface; coor_node_h_y_down_interface(end) + dy(ss-1);...
            fout.coor_nodes.y(index_la-grid.n_y(ss-1)+1); coor_node_h_y_down_interface; coor_node_h_y_down_interface(end) + dy(ss-1);...
            coor_node_h_y.down];
        
        %two nodes downstream
        index_node_hdown_interface_1 = [index_la-grid.n_y(ss-1)+1; index_node_hdown_interface_1; index_la-grid.n_y(ss-1)+1];
        index_node_hdowndown_interface_1 = index_node_hdown_interface_1 + 1;
        index_node_hdowndown_interface_1(end-2:end-1) = index_la-grid.n_y(ss-1)+1; 
        index_node_hdown_interface_2 = [index_la-2*grid.n_y(ss-1)+1; index_node_hdown_interface_2; index_la-2*grid.n_y(ss-1)+1];
        index_node_hdowndown_interface_2 = index_node_hdown_interface_2 + 1;
        index_node_hdowndown_interface_2(end-2:end-1) = index_la-2*grid.n_y(ss-1)+1; 
        
        index_node_hdowndown = [index_node_hdowndown_interface_2; index_node_hdowndown_interface_1; index_node_hdowndown];
        
        coor_node_h_y_down_interface = [fout.coor_nodes.y(index_la-grid.n_y(ss-1)+1); coor_node_h_y_down_interface; coor_node_h_y_down_interface(end) + dy(ss-1)];
        coor_node_h_y_downdown_interface = coor_node_h_y_down_interface + dy(ss-1);
        coor_node_h_y.downdown = [coor_node_h_y_downdown_interface; coor_node_h_y_downdown_interface; coor_node_h_y.downdown];
        
        %one node upstream
        index_node_hup_interface_1 = reshape(repmat((index_la-grid.n_y(ss-1)+1:index_la-1),[2 1]),[ ],1);
        index_node_hup_interface_2 = reshape(repmat((index_la-2*grid.n_y(ss-1)+1:index_la-grid.n_y(ss-1)-1),[2 1]),[ ],1);
        index_node_hup = [index_la-grid.n_y(ss-1);index_node_hup_interface_2; index_la-grid.n_y(ss-1); ...
            index_la; index_node_hup_interface_1;index_la;...
            index_node_hup];
        
        coor_node_h_y_up_interface = reshape(repmat((fout.coor_nodes.y(index_la-grid.n_y(ss-1)+1:index_la-1)).',[2 1]),[ ],1);
        coor_node_h_y.up = [coor_node_h_y_up_interface(1)-dy(ss-1); coor_node_h_y_up_interface; fout.coor_nodes.y(index_la);...
            coor_node_h_y_up_interface(1)-dy(ss-1); coor_node_h_y_up_interface; fout.coor_nodes.y(index_la); ...
            coor_node_h_y.up];
        
        %two nodes upstream
        index_node_hup_interface_1 = [index_la; index_node_hup_interface_1;index_la];
        index_node_hupup_interface_1 = index_node_hup_interface_1 - 1;
        index_node_hupup_interface_1(2:3) = index_la;
        index_node_hup_interface_2 = [index_la-grid.n_y(ss-1);index_node_hup_interface_2; index_la-grid.n_y(ss-1)];
        index_node_hupup_interface_2 = index_node_hup_interface_2 - 1;
        index_node_hupup_interface_2(2:3) = index_la-grid.n_y(ss-1);
        
        index_node_hupup = [index_node_hupup_interface_2; index_node_hupup_interface_1; index_node_hupup];
        
        coor_node_h_y_up_interface = [coor_node_h_y_up_interface(1)-dy(ss-1); coor_node_h_y_up_interface; fout.coor_nodes.y(index_la)];
        coor_node_h_y_upup_interface = coor_node_h_y_up_interface - dy(ss-1);
        coor_node_h_y.upup = [coor_node_h_y_upup_interface; coor_node_h_y_upup_interface; coor_node_h_y.upup ];
        
    end
    
    %must add all layers together
    fout.coor_nodes.y = [fout.coor_nodes.y; coor_nodes.y];
    fout.coor_nodes.z = [fout.coor_nodes.z; coor_nodes.z];
    fout.Delta_y_cell = [fout.Delta_y_cell; Delta_y_cell];
    fout.Delta_z_cell = [fout.Delta_z_cell; Delta_z_cell];
    fout.delta_y_vflux = [fout.delta_y_vflux; delta_y_vflux ];
    
    fout.coor_hor_cell_edges.z = [fout.coor_hor_cell_edges.z; coor_hor_cell_edges.z];
    fout.coor_hor_cell_edges.y = [fout.coor_hor_cell_edges.y; coor_hor_cell_edges.y];
    fout.coor_ver_cell_edges.z = [fout.coor_ver_cell_edges.z; coor_ver_cell_edges.z];
    fout.coor_ver_cell_edges.y = [fout.coor_ver_cell_edges.y; coor_ver_cell_edges.y];
    
    fout.up_node.hor = [fout.up_node.hor; up_node.hor];
    fout.up_node.vert = [fout.up_node.vert; up_node.ver];
    fout.down_node.hor = [fout.down_node.hor; down_node.hor];
    fout.down_node.vert = [fout.down_node.vert; down_node.ver];
    
    fout.down_node_virtual.vert = [fout.down_node_virtual.vert; down_node_virtual.ver];
    fout.up_node_virtual.vert = [fout.up_node_virtual.vert; up_node_virtual.ver];
    fout.downdown_node_virtual.vert = [fout.downdown_node_virtual.vert; downdown_node_virtual.ver];
    fout.upup_node_virtual.vert = [fout.upup_node_virtual.vert; upup_node_virtual.ver];
   
    fout.coor_node_h_y.up = [fout.coor_node_h_y.up; coor_node_h_y.up];
    fout.coor_node_h_y.down = [fout.coor_node_h_y.down; coor_node_h_y.down];
    fout.index_node_hdown = [fout.index_node_hdown; index_node_hdown];
    fout.index_node_hup = [fout.index_node_hup; index_node_hup];
    fout.coor_node_h_y.upup = [fout.coor_node_h_y.upup; coor_node_h_y.upup ];
    fout.coor_node_h_y.downdown = [fout.coor_node_h_y.downdown; coor_node_h_y.downdown];
    fout.index_node_hdowndown = [fout.index_node_hdowndown;index_node_hdowndown];
    fout.index_node_hupup = [fout.index_node_hupup; index_node_hupup];
    
    fout.coor_node_virtual.y = [fout.coor_node_virtual.y; coor_node_virtual.y];
    fout.coor_down_node_virtual_vert = [fout.coor_down_node_virtual_vert; coor_down_node_virtual_vert];
    fout.coor_up_node_virtual_vert = [fout.coor_up_node_virtual_vert; coor_up_node_virtual_vert];
    fout.coor_downdown_node_virtual_vert = [fout.coor_downdown_node_virtual_vert; coor_downdown_node_virtual_vert];
    fout.coor_upup_node_virtual_vert = [fout.coor_upup_node_virtual_vert; coor_upup_node_virtual_vert];
end
fout.n_edges.vert = length(fout.up_node.vert);
% build connect matrix for vertical fluxes
connect_ver = sparse(fout.n_nodes.tot, fout.n_edges.vert);
for k = 1:fout.n_edges.vert
    index_up = fout.up_node.vert(k) ; 
    index_down = fout.down_node.vert(k); 
    connect_ver = connect_ver + sparse(index_up,k,1,fout.n_nodes.tot, fout.n_edges.vert);
    connect_ver = connect_ver + sparse(index_down,k,-1,fout.n_nodes.tot, fout.n_edges.vert);
end
fout.connect_ver = connect_ver;

fout.n_edges.hor = length(fout.up_node.hor);
% build connect matrix for horizontal fluxes
connect_hor = sparse(fout.n_nodes.tot, fout.n_edges.hor);
for k = 1:fout.n_edges.hor
    index_up = fout.up_node.hor(k) ; 
    index_down = fout.down_node.hor(k); 
    connect_hor = connect_hor + sparse(index_up,k,1,fout.n_nodes.tot, fout.n_edges.hor);
    connect_hor = connect_hor + sparse(index_down,k,-1,fout.n_nodes.tot, fout.n_edges.hor);
end
fout.connect_hor = connect_hor;
end

function fout = grid_T(grid, dy,dz,flag)
%set flag to 0 for T_bed grid, set to 1 for T grid both in the ice and in
%the bed or T,N grid in the ice

% save/define key grid parameters
fout.n_nodes.vert = grid.n_v;                                              %list n_layer-by-one
fout.n_nodes.hor = grid.n_y;                                               %list n_layer-by-one
fout.n_nodes.tot = sum(fout.n_nodes.hor.*fout.n_nodes.vert);

if flag==1
    fout.n_nodes.tot_virtual = sum(fout.n_nodes.hor.*fout.n_nodes.vert) + sum(2*fout.n_nodes.hor(2:end));
elseif flag==0
    fout.n_nodes.tot_virtual = sum(fout.n_nodes.hor.*fout.n_nodes.vert) + sum(2*fout.n_nodes.hor(1:end-1));
end
fout.n_edges.hor = sum((fout.n_nodes.hor).*fout.n_nodes.vert);             %THIS IS FOR A PERIODIC GRID!!
fout.n_edges.vert = 0;                                                     %build along with the grid

%list of nodes where boundary conditions apply
fout.bdy_nodes.top = (1:fout.n_nodes.hor(1))';
fout.bdy_nodes.bottom = (fout.n_nodes.tot - fout.n_nodes.hor(end)+1:fout.n_nodes.tot)';
if flag == 1 && length(dz)>1
    fout.bdy_nodes.bed = sum(fout.n_nodes.vert(1:length(fout.n_nodes.vert)/2-1).*fout.n_nodes.hor(1:length(fout.n_nodes.vert)/2-1)) +...
        1/2*(fout.n_nodes.vert(length(fout.n_nodes.vert)/2).*fout.n_nodes.hor(length(fout.n_nodes.vert)/2)) +1:...
        sum(fout.n_nodes.vert(1:length(fout.n_nodes.vert)/2-1).*fout.n_nodes.hor(1:length(fout.n_nodes.vert)/2-1)) +...
        1/2*(fout.n_nodes.vert(length(fout.n_nodes.vert)/2).*fout.n_nodes.hor(length(fout.n_nodes.vert)/2)) + ...
        fout.n_nodes.hor(length(fout.n_nodes.vert)/2);         %%DOUBLE CHECK!!
else
    fout.bdy_nodes.bed = fout.bdy_nodes.bottom;
end
fout.bdy_nodes.inflow = [];                                                %build along with the grid
fout.bdy_nodes.outflow = [];                                               %build along with the grid

%build the grid a layer at a time, starting from the top

%initialization
fout.coor_nodes.y = [];
fout.coor_nodes.z = [];
fout.coor_node_virtual.y = [];  
fout.coor_node_virtual.z = [];
fout.Delta_y_cell = [];
fout.Delta_y_edge = [];
fout.Delta_z_cell = [];
fout.Delta_z_cell_volume = [];  %the difference between Delta_z_cell_volume and Delta_z_cell is that the former is the right cell vertical height for T(ice) and all the flow variables,
                                %because it assumes half of the height for
                                %cells at the bottom and at the top.
                                %Delta_z_cell must be used for the vertical
                                %averaging, and considers full cells at the
                                %top and bottom.
fout.delta_y_vflux = [];

fout.coor_hor_cell_edges.z = [];
fout.coor_hor_cell_edges.y = [];
fout.up_node.hor = [];
fout.up_node.vert = [];
fout.up_node.adv_ver = [];
fout.up_node.adv_hor = [];
fout.down_node.hor = [];
fout.down_node.vert = [];
fout.down_node.adv_ver = [];
fout.down_node.adv_hor = [];
fout.down_node_virtual.vert = [];
fout.up_node_virtual.vert = [];

%lateral interpolation
fout.coor_node_h_y.down = [];
fout.coor_node_h_y.downdown = [];
fout.index_node_h.down = [];
fout.index_node_h.up = [];
fout.index_node_h.downdown = [];
fout.index_node_h.upup = [];
fout.coor_node_h_y.upup = [];
fout.coor_node_h_y.up = [];

%horizontal grid
fout.up_edge.hor = [];
fout.down_edge.hor = [];

for ss = 1:grid.n_layers
    
    index_la = sum(grid.n_y(1:ss-1).*grid.n_v(1:ss-1));                    %index of last (bottom right) node in layer above current in absolute terms, that is, in the global indexing system for real nodes
    if ss == 1
        index_la_virtual = index_la;                                       %index of last (bottom right) VIRTUAL node in layer above current in absolute terms, that is, in the global indexing system for virtual nodes
    elseif ss>1
        if dz(2)<dz(1)
            index_la_virtual = sum(grid.n_y(1:ss-1).*grid.n_v(1:ss-1)) + sum(2*grid.n_y(2:ss));
        else
            index_la_virtual = sum(grid.n_y(1:ss-1).*grid.n_v(1:ss-1)) + sum(2*grid.n_y(1:ss-1));
        end
    end
    n_nodes_layer = grid.n_y(ss)*grid.n_v(ss);
    n_edges_layer_ver = grid.n_y(ss)*(grid.n_v(ss)-1);
    n_nodes_layer_virtual = n_nodes_layer;
    index_nodes = index_la+(1:n_nodes_layer).';
    index_nodes_virtual = index_la_virtual+(1:n_nodes_layer).';
    index_nodes_layer = (1:n_nodes_layer).';
    index_nodes_layer_virtual = (1:n_nodes_layer_virtual).';
    
    %number of ver network edges in layers above current, after horizontal
    %averaging of ver velocities
    if ss >1
        %network for u_vert_av nodes
        index_ve_a = sum(grid.n_y(2:ss-1).*grid.n_v(2:ss-1))+ grid.n_y(1).*(grid.n_v(1) -1)  ; %includes interface
        n_edges_layer_a = n_nodes_layer; %includes interface
        
        %horizontal advection averaging, network for u_vert (.._ve_nodes) and
        %for u_vert_av (.._ve_edges)
        index_ve_nodes = grid.n_y(1).*(grid.n_v(1) -1) + 3*grid.n_y(1) + sum(grid.n_y(2:ss-1).*(grid.n_v(2:ss-1)-1) + 3*grid.n_y(2:ss-1)) ;%includes lower interface ; this index labels nodes of the u_vert network *before* lateral averaging
        index_ve_edges = grid.n_y(1).*(grid.n_v(1) -1) + 2*grid.n_y(1) + sum(grid.n_y(2:ss-1).*grid.n_v(2:ss-1)+ 2* grid.n_y(2:ss-1));%includes lower interface; this index labels nodes of the u_vert network *after* lateral averaging
        n_edges_layer = grid.n_y(ss)*(grid.n_v(ss)-1); %only interior part of the layer
    elseif ss == 1
        index_ve_a = 0;
        n_edges_layer_a = grid.n_y(ss)*(grid.n_v(ss)-1); 
        
        index_ve_nodes = 0;
        index_ve_edges = 0;
        n_edges_layer = grid.n_y(ss)*(grid.n_v(ss)-1);       
    end
    
     index_vedges_a = index_ve_a + (1:n_edges_layer_a).';
     index_vedges_layer_a = (1:n_edges_layer_a).';
     
     index_vedges_edges = index_ve_edges + (1:n_edges_layer).';
     index_vedges_nodes = index_ve_nodes + (1:n_edges_layer).';
     index_vedges_edges_layer = (1:n_edges_layer).';
     
    
    %inflow and outflow boundary
    fout.bdy_nodes.inflow = [fout.bdy_nodes.inflow; index_la + (1:grid.n_y(ss):n_nodes_layer).'];
    fout.bdy_nodes.outflow = [fout.bdy_nodes.outflow; index_la+(grid.n_y(ss):grid.n_y(ss):n_nodes_layer).']; 
    
    %cell centres locations (list n_nodes_layer-by-one)
    if flag == 1
        z_top = grid.bd_z.ice-sum(grid.layer_height(1:ss-1));
        z_bottom = grid.bd_z.ice-sum(grid.layer_height(1:ss));
    elseif flag == 0
        if ss== 1
            z_top = grid.bd_z.ice;
        elseif ss>1
            z_top = grid.bd_z.ice + dz(1)/2-sum(grid.layer_height(1:ss-1)) - dz(ss)/2;   %this is the coordinate of the topmost node for each layer
        end
        z_bottom = z_top + dz(ss)/2 - grid.layer_height(ss);             %this is the coordinate of the lowermost node for each layer
    end
    
    if ss<length(dz) && flag == 1
        z_nodes_line = (z_top:-dz(ss):z_bottom+dz(ss)).';
    elseif ss == length(dz) && flag == 1
        z_nodes_line = (z_top:-dz(ss):z_bottom).';
    elseif flag == 0
        z_nodes_line = (z_top:-dz(ss):z_bottom).';
    end
    
    coor_nodes.z = reshape(repmat(z_nodes_line,[1,grid.n_y(ss)]).', n_nodes_layer, []);
    
    y_nodes_line = (dy(ss):dy(ss):grid.bd_y).';
    coor_nodes.y = reshape(repmat(y_nodes_line,[grid.n_v(ss),1]).', n_nodes_layer, []);
    coor_node_virtual.y = coor_nodes.y;
    coor_node_virtual.z = coor_nodes.z;
    
    delta_y = dy(ss)*ones(n_nodes_layer,1);   
    delta_y_edge = dy(ss)*ones(n_nodes_layer,1);

    Delta_z_cell = dz(ss)*ones(n_nodes_layer,1);
    Delta_z_cell_volume = dz(ss)*ones(n_nodes_layer,1);
    if ss==1 && flag==1
        Delta_z_cell_volume(1: grid.n_y(ss))= dz(ss)/2*ones(grid.n_y(ss),1);
    elseif ss==length(dy) && flag==1
        Delta_z_cell_volume(end-grid.n_y(ss)+1: end)= dz(ss)/2*ones(grid.n_y(ss),1);
    end
    
    if ss>1 && flag == 1
    Delta_z_cell(1:grid.n_y(ss)) = 1/2*dz(ss) + 1/2*dz(ss-1);
    Delta_z_cell_volume(1:grid.n_y(ss)) = 1/2*dz(ss) + 1/2*dz(ss-1);
    end
    
    % length of horizontal cell edges for computation of vertical fluxes (only interior vertical edges
    %n.b.: the row of nodes that lies on a layer boundary belongs to the lower layer)

    delta_y_vflux = dy(ss)*ones((grid.n_v(ss)-1)*grid.n_y(ss),1);
    
    %coordinate of horizontal cell edges (same as network vertical edges), including only the edges internal
    %to the layer. Edges on the bottom boundary are
    %added in the for loop below. 
    %NOTE: node numbering is from upper left corner to bottom
    %right corner along rows; same for edges) 
    
    n_edges_int = (grid.n_v(ss)-1)*grid.n_y(ss);
    if ss<length(dz) && flag == 1
        coor_hor_cell_edges.z = (z_top-(1-1/2)*dz(ss):-dz(ss):z_bottom+3*dz(ss)/2).';
    elseif ss==length(dz) && flag == 1
        coor_hor_cell_edges.z = (z_top-(1-1/2)*dz(ss):-dz(ss):z_bottom+dz(ss)/2).';
    elseif flag == 0 
        coor_hor_cell_edges.z = (z_top-(1-1/2)*dz(ss):-dz(ss):z_bottom+dz(ss)).';
    end
    coor_hor_cell_edges.z = reshape(repmat(coor_hor_cell_edges.z,[1,grid.n_y(ss)]).', n_edges_int,[]);
   
    coor_hor_cell_edges.y = y_nodes_line;
    coor_hor_cell_edges.y = reshape(repmat(coor_hor_cell_edges.y.',[grid.n_v(ss)-1,1]).',n_edges_int,[]);
    
    
    %upnodes/downnodes of network edges inside the layer, excluding top, and bottom
    %boundaries. Assumes a PERIODIC GRID IN Y
    
    %NOTE: network edges numbering is from upper left corner to bottom
    %right corner along rows, both for horizontal and vertical edges
    
    %hor network edges (up_node and down_node list up node and down node
    %for each edge; numbering along rows, edge #1 connects node#1 to node#2)

    up = reshape((circshift(reshape(index_nodes_layer, [grid.n_y(ss) grid.n_v(ss)]).',1,2)).',n_nodes_layer,1);
    down = index_nodes_layer;
    up_node.hor = index_nodes(up); 
    down_node.hor = index_nodes(down);
    
    up_edge.hor = index_nodes(down);
    down_edge.hor = reshape(circshift(reshape(index_nodes, [grid.n_y(ss) grid.n_v(ss)]).',-1,2).',n_nodes_layer,1);
    
    %ver network edges (excluding edges connecting adjacent layers)
    up = index_nodes_layer; up(1:grid.n_y(ss))=[];
    down = index_nodes_layer; down(n_nodes_layer-grid.n_y(ss)+1:n_nodes_layer)=[];
    up_node.ver = index_nodes(up); 
    down_node.ver = index_nodes(down);
    
    up = index_nodes_layer_virtual; up(1:grid.n_y(ss))=[];
    down = index_nodes_layer_virtual; down(n_nodes_layer_virtual-grid.n_y(ss)+1:n_nodes_layer_virtual)=[];
    up_node_virtual.ver = index_nodes_virtual(up); 
    down_node_virtual.ver = index_nodes_virtual(down); 
    
    % up_nodes and down_nodes for computation of the vertical advective
    % flux after horizontal averaging
    up = index_vedges_layer_a; up(1:grid.n_y(ss))=[];
    down = index_vedges_layer_a; down(n_edges_layer_a-grid.n_y(ss)+1:n_edges_layer_a)=[];
    up_node_adv_ver_av =  index_vedges_a(up); 
    down_node_adv_ver_av = index_vedges_a(down);  
    
    % up_nodes and down_nodes for computation of the horizontally averaged
    % velocities
    
    up = index_vedges_edges_layer; 
    down = index_vedges_edges_layer; 
    up_node_adv_hor =  index_vedges_nodes(up); 
    down_node_adv_hor = index_vedges_nodes(down);  
    
    % virtual node network, nodes within layers
  
    index_node_hdown = index_nodes;
    coor_node_h_y.down = coor_nodes.y;
    
    %one further downstream
    index_node_hdowndown = reshape((circshift(reshape(index_node_hdown, [grid.n_y(ss) grid.n_v(ss)]).',-1,2)).',n_nodes_layer_virtual,1);
    coor_intermediate = (circshift(reshape(coor_node_h_y.down, [grid.n_y(ss) grid.n_v(ss)]).',-1,2));
    coor_intermediate(:,end) = coor_intermediate(:,end-1)+ dy(ss);
    coor_node_h_y.downdown = reshape(coor_intermediate.',n_nodes_layer_virtual,1);
    
    %one upstream
    index_node_hup = reshape((circshift(reshape(index_node_hdown, [grid.n_y(ss) grid.n_v(ss)]).',1,2)).',n_nodes_layer_virtual,1);
    coor_intermediate = (circshift(reshape(coor_node_h_y.down, [grid.n_y(ss) grid.n_v(ss)]).',1,2));
    coor_intermediate(:,1) = coor_intermediate(:,2)- dy(ss);
    coor_node_h_y.up = reshape(coor_intermediate.',n_nodes_layer_virtual,1);
    
    %two upstream
    index_node_hupup = reshape((circshift(reshape(index_node_hdown, [grid.n_y(ss) grid.n_v(ss)]).',2,2)).',n_nodes_layer_virtual,1);
    coor_intermediate = (circshift(reshape(coor_node_h_y.down, [grid.n_y(ss) grid.n_v(ss)]).',2,2));
    if grid.n_y(ss)>=4
        coor_intermediate(:,2) = coor_intermediate(:,3)- dy(ss);
        coor_intermediate(:,1) = coor_intermediate(:,2)- dy(ss);
    end
    coor_node_h_y.upup = reshape(coor_intermediate.',n_nodes_layer_virtual,1);
        
    
    %% add ver network edges that connect lower to upper layer
    if ss>1 && (dy(ss)<dy(ss-1)) %refinement
        %down nodes (T nodes always remain on domain boundaries, so upper cells are split three ways, except for first and last cell)
        down_node_interface = reshape(repmat((index_la-grid.n_y(ss-1)+1:index_la),[3 1]),[ ],1);
        down_node.ver = [circshift(down_node_interface, 1,1);down_node.ver];
        
        index_interface = 1: 3*grid.n_y(ss-1);
        index_single_virtual = 3:3:length(index_interface);
        index_extra_virtual = sort([1:3:length(index_interface), 2:3:length(index_interface)]);
        down_node_virtual_interface = sparse(index_extra_virtual,ones(length(index_extra_virtual),1), (index_la_virtual-2*grid.n_y(ss)+1:index_la_virtual-grid.n_y(ss)).',length(index_interface),1) + ...
            sparse(index_single_virtual,ones(length(index_single_virtual),1), (index_la_virtual-5*grid.n_y(ss-1)+1:index_la_virtual-4*grid.n_y(ss-1)).',length(index_interface),1);
        down_node_virtual.ver = [down_node_virtual_interface; down_node_virtual.ver];
        
        down_node_adv_ver_interface = (index_ve_a-grid.n_y(ss-1)+1:index_ve_a).';
        down_node_adv_ver_av = [down_node_adv_ver_interface; down_node_adv_ver_av];
        
        index_single_hor = 2:2:grid.n_y(ss);
        index_split_up = 1:2:grid.n_y(ss);
        index_split_down = 1:2:grid.n_y(ss);
        
        down_node_adv_hor_interface = sparse(index_single_hor.', ones(length(index_single_hor),1), index_ve_nodes - 3*grid.n_y(ss-1) + [3:3:3*grid.n_y(ss-1)], grid.n_y(ss),1) +...
            sparse(index_split_down.', ones(length(index_split_down),1), index_ve_nodes - 3*grid.n_y(ss-1) + [2:3:3*grid.n_y(ss-1)], grid.n_y(ss),1)  ;
        
        up_node_adv_hor_interface = sparse(index_single_hor.', ones(length(index_single_hor),1), index_ve_nodes- 3*grid.n_y(ss-1) + [3:3:3*grid.n_y(ss-1)], grid.n_y(ss),1) +...
            sparse(index_split_up.', ones(length(index_split_up),1), index_ve_nodes - 3*grid.n_y(ss-1) + [1:3:3*grid.n_y(ss-1)], grid.n_y(ss),1)  ;
        down_node_adv_hor = [down_node_adv_hor_interface; down_node_adv_hor];
        up_node_adv_hor = [up_node_adv_hor_interface; up_node_adv_hor];
        
        % up nodes (must be in the same number as down_node_interface)
        
        %regular network
        index_single = 3:3:length(down_node_interface);
        up_node_interface = sparse(index_single.', ones(length(index_single),1), index_la+(2:2:grid.n_y(ss)), length(down_node_interface),1);
        
        index_split_up = 1:3:length(down_node_interface);
        up_node_interface = up_node_interface + sparse(index_split_up.', ones(length(index_split_up),1),  index_la+(1:2:grid.n_y(ss)), length(down_node_interface),1);
        
        index_split_down = 2:3:length(down_node_interface);
        up_node_interface = up_node_interface + sparse(index_split_down.', ones(length(index_split_down),1),  index_la+(1:2:grid.n_y(ss)), length(down_node_interface),1);
        
        up_node.ver = [full(up_node_interface); up_node.ver];
        
        up_node_adv_ver_interface = (index_ve_a +2:2: index_ve_a +grid.n_y(ss)).';
        up_node_adv_ver_av = [up_node_adv_ver_interface; up_node_adv_ver_av];
        
        %virtual network:
        index_interface = 1: 3*grid.n_y(ss-1);
        index_single_virtual = 3:3:length(index_interface);
        index_extra_virtual = sort([1:3:length(index_interface), 2:3:length(index_interface)]);
        up_node_virtual_interface = sparse(index_extra_virtual,ones(length(index_extra_virtual),1), (index_la_virtual-1*grid.n_y(ss)+1:index_la_virtual).',length(index_interface),1) + ...
            sparse(index_single_virtual,ones(length(index_single_virtual),1), (index_la_virtual+2:2:index_la_virtual+grid.n_y(ss)).',length(index_interface),1);
        up_node_virtual.ver = [up_node_virtual_interface; up_node_virtual.ver];  
        
        %add location (x,z) of hor cell edges
        coor_interface_y = sparse(index_single.', ones(length(index_single),1), (dy(ss-1):dy(ss-1):grid.bd_y).', length(down_node_interface),1)+...
            sparse(index_split_up.', ones(length(index_split_up),1), 3/4*dy(ss):2*dy(ss):grid.bd_y, length(down_node_interface),1)+...
            sparse(index_split_down.', ones(length(index_split_down),1), 5/4*dy(ss):2*dy(ss):grid.bd_y, length(down_node_interface),1);
        
        coor_hor_cell_edges.y = [full(coor_interface_y); coor_hor_cell_edges.y];
        coor_hor_cell_edges.z = [(z_top+dz(ss-1)/2)*ones(length(up_node_interface),1);coor_hor_cell_edges.z];
        
        %add surface crossed by vertical network edge
        delta_y_vflux_interface = sparse(index_single.', ones(length(index_single),1), dy(ss)*ones(length(index_single),1), length(down_node_interface),1)+...
            sparse(index_split_up.', ones(length(index_split_up),1), 1/2*dy(ss)*ones(length(index_split_up),1), length(down_node_interface),1)+...
            sparse(index_split_down.', ones(length(index_split_down),1), 1/2*dy(ss)*ones(length(index_split_down),1), length(down_node_interface),1);
        
        delta_y_vflux = [full(delta_y_vflux_interface); delta_y_vflux];
        
        %horizontal network edge intersected by vertical network edge in the
        %downstream direction (list n_ver_network_edges-by-one)   
        
        index_split_up_virtual = 1:2:grid.n_y(ss);
        index_split_down_virtual = 2:2:grid.n_y(ss);
        %down node
        index_node_hddown_interface = reshape(repmat((index_la-grid.n_y(ss-1)+1:index_la),[2 1]),[ ],1);
        index_node_hudown_interface = sparse(index_split_up_virtual.', ones(length(index_split_up_virtual),1), index_la+(1:2:grid.n_y(ss)), grid.n_y(ss),1)+...
            sparse(index_split_down_virtual.', ones(length(index_split_down_virtual),1), index_la+(2:2:grid.n_y(ss)), grid.n_y(ss),1);
        index_node_hdown = [index_node_hddown_interface; index_node_hudown_interface; index_node_hdown];
        
        coor_node_h_y_d_down_interface = reshape(repmat((fout.coor_nodes.y(index_la-grid.n_y(ss-1)+1:index_la)).',[2 1]),[ ],1);
        coor_node_h_y_d_up_interface = sparse(index_split_up_virtual.', ones(length(index_split_up_virtual),1), (dy(ss):2*dy(ss):grid.bd_y).' , grid.n_y(ss),1)+...
            sparse(index_split_down_virtual.', ones(length(index_split_down_virtual),1), (2*dy(ss):2*dy(ss):grid.bd_y).' , grid.n_y(ss),1);
        coor_node_h_y.down = [coor_node_h_y_d_down_interface; coor_node_h_y_d_up_interface; coor_node_h_y.down];
        
        %upnode
        index_node_hdup_interface = sparse(index_split_up_virtual.', ones(length(index_split_up_virtual),1), circshift(index_la-grid.n_y(ss-1)+(1:grid.n_y(ss-1)),1,2) ,  grid.n_y(ss),1)+...
            sparse(index_split_down_virtual.', ones(length(index_split_down_virtual),1), circshift(index_la-grid.n_y(ss-1)+(1:grid.n_y(ss-1)),1),grid.n_y(ss),1);
        index_node_huup_interface = sparse(index_split_up_virtual.', ones(length(index_split_up_virtual),1), circshift(index_la+(2:2:grid.n_y(ss)),1,2) ,grid.n_y(ss),1)+...
            sparse(index_split_down_virtual.', ones(length(index_split_down_virtual),1), index_la+(1:2:grid.n_y(ss)),grid.n_y(ss),1);
        index_node_hup = [index_node_hdup_interface; index_node_huup_interface; index_node_hup];
        
        coor_node_h_y_u_down_interface = coor_node_h_y_d_down_interface - dy(ss-1);
        coor_node_h_y_u_up_interface = coor_node_h_y_d_up_interface -dy(ss);
        coor_node_h_y.up = [coor_node_h_y_u_down_interface; coor_node_h_y_u_up_interface; coor_node_h_y.up];
        
        %one node further up
        index_node_hup_interface = [index_node_hdup_interface.'; index_node_huup_interface.'];
        index_node_hupup_interface = reshape([circshift(index_node_hup_interface(1,:), 2, 2);circshift(index_node_hup_interface(2,:), 1, 2)].', [], 1);
        index_node_hupup = [index_node_hupup_interface; index_node_hupup];
        coor_node_hup_interface =  [coor_node_h_y_u_down_interface.'; coor_node_h_y_u_up_interface.'];
        coor_node_hupup_interface = reshape([coor_node_hup_interface(1,:) - dy(ss-1); coor_node_hup_interface(2,:) - dy(ss)].',[],1);
        coor_node_h_y.upup = [coor_node_hupup_interface; coor_node_h_y.upup];
        
        %one node further down
        index_node_hdown_interface = [index_node_hddown_interface.'; index_node_hudown_interface.'];
        index_node_hdowndown_interface = reshape([circshift(index_node_hdown_interface(1,:), -2, 2); circshift(index_node_hdown_interface(2,:), -1, 2)].', [], 1);
        index_node_hdowndown = [index_node_hdowndown_interface; index_node_hdowndown];
        coor_node_hdown_interface = [coor_node_h_y_d_down_interface.'; coor_node_h_y_d_up_interface.'];
        coor_node_hdowndown_interface = reshape([coor_node_hdown_interface(1,:) +  dy(ss-1); coor_node_hdown_interface(2,:) +  dy(ss)].', [],1);
        coor_node_h_y.downdown = [coor_node_hdowndown_interface; coor_node_h_y.downdown];
        
        %coordinates of virtual nodes
        y_virtual_interface = sort([3/4*dy(ss): 2*dy(ss):grid.bd_y,5/4*dy(ss): 2*dy(ss):grid.bd_y]).';
        coor_node_virtual.y = [y_virtual_interface; y_virtual_interface; coor_node_virtual.y];
        coor_node_virtual.z = [(z_top+ dz(ss-1))*ones(grid.n_y(ss),1);z_top*ones(grid.n_y(ss),1); coor_node_virtual.z];
        
    elseif ss>1 && (dy(ss)>dy(ss-1)) %coarsening as we move down along the grid
        %up nodes
        up_node_interface = reshape(repmat((index_la+1:index_la+grid.n_y(ss)),[3 1]),[ ],1);  %OK
        %up_node_interface(1) = []; up_node_interface(end) = [];
        up_node.ver = [circshift(up_node_interface, 1,1);up_node.ver];
        
        % down nodes (must be in the same number as down_node_interface)
        
        index_single = 3:3:length(up_node_interface);
        down_node_interface = sparse(index_single.', ones(length(index_single),1), index_la-grid.n_y(ss-1)+(2:2:grid.n_y(ss-1)), length(up_node_interface),1);
        
        index_split_up = 1:3:length(up_node_interface);
        down_node_interface = down_node_interface + sparse(index_split_up.', ones(length(index_split_up),1),  index_la-grid.n_y(ss-1)+(1:2:grid.n_y(ss-1)), length(up_node_interface),1);
        
        index_split_down = 2:3:length(up_node_interface);
        down_node_interface = down_node_interface + sparse(index_split_down.', ones(length(index_split_down),1),  index_la-grid.n_y(ss-1)+(1:2:grid.n_y(ss-1)), length(up_node_interface),1);
        
        down_node.ver = [full(down_node_interface); down_node.ver];
        
        index_interface = 1: 3*grid.n_y(ss);
        index_single_virtual = 3:3:length(index_interface);
        index_extra_virtual = sort([1:3:length(index_interface), 2:3:length(index_interface)]);
        down_node_virtual_interface = sparse(index_extra_virtual,ones(length(index_extra_virtual),1), (index_la_virtual-2*grid.n_y(ss-1)+1:index_la_virtual-1*grid.n_y(ss-1)).',length(index_interface),1) + ...
            sparse(index_single_virtual,ones(length(index_single_virtual),1), (index_la_virtual-3*grid.n_y(ss-1)+2:2:index_la_virtual-2*grid.n_y(ss-1)).',length(index_interface),1);
        down_node_virtual.ver = [down_node_virtual_interface; down_node_virtual.ver];
        
        
        %% CORRECTION, 4/1
%         up_node_virtual_interface = sparse(index_extra_virtual,ones(length(index_extra_virtual),1), circshift((index_la_virtual-1*grid.n_y(ss-1)+1:index_la_virtual).',1,1),length(index_interface),1) + ...
%            sparse(index_single_virtual,ones(length(index_single_virtual),1), (index_la_virtual+1:index_la_virtual+grid.n_y(ss)).',length(index_interface),1);
       up_node_virtual_interface = sparse(index_extra_virtual,ones(length(index_extra_virtual),1), (index_la_virtual-1*grid.n_y(ss-1)+1:index_la_virtual).',length(index_interface),1) + ...
            sparse(index_single_virtual,ones(length(index_single_virtual),1), (index_la_virtual+1:index_la_virtual+grid.n_y(ss)).',length(index_interface),1);
         %%   
        up_node_virtual.ver = [up_node_virtual_interface; up_node_virtual.ver];  
        
        %add location (x,z) of hor cell edges
        coor_interface_y = sparse(index_single.', ones(length(index_single),1), (2*dy(ss-1):2*dy(ss-1):grid.bd_y).', length(up_node_interface),1)+...
            sparse(index_split_up.', ones(length(index_split_up),1), 3/4*dy(ss-1):2*dy(ss-1):grid.bd_y, length(up_node_interface),1)+...
            sparse(index_split_down.', ones(length(index_split_down),1), 5/4*dy(ss-1):2*dy(ss-1):grid.bd_y, length(up_node_interface),1);
        
        coor_hor_cell_edges.y = [full(coor_interface_y); coor_hor_cell_edges.y];
        %% CORRECTION, 4/1
        coor_hor_cell_edges.z = [(z_top+dz(ss)/2)*ones(length(up_node_interface),1);coor_hor_cell_edges.z];
        %OLD coor_hor_cell_edges.z = [(z_top+dz(ss-1)/2)*ones(length(up_node_interface),1);coor_hor_cell_edges.z];
        %%
        
        %add surface crossed by vertical network edge
        delta_y_vflux_interface = sparse(index_single.', ones(length(index_single),1), dy(ss-1)*ones(length(index_single),1), length(up_node_interface),1)+...
            sparse(index_split_up.', ones(length(index_split_up),1), 1/2*dy(ss-1)*ones(length(index_split_up),1), length(up_node_interface),1)+...
            sparse(index_split_down.', ones(length(index_split_down),1), 1/2*dy(ss-1)*ones(length(index_split_down),1), length(up_node_interface),1);
        
        delta_y_vflux = [delta_y_vflux_interface; delta_y_vflux];
        
        %add virtual nodes on the interface
        
        index_split_up_virtual = 1:2:grid.n_y(ss-1);
        index_split_down_virtual = 2:2:grid.n_y(ss-1);
        %down node 
        index_node_hudown_interface = reshape(repmat((index_la+1:index_la+grid.n_y(ss)),[2 1]),[ ],1);
        index_node_hddown_interface = sparse(index_split_up_virtual.', ones(length(index_split_up_virtual),1), index_la - grid.n_y(ss-1)+(1:2:grid.n_y(ss-1)), grid.n_y(ss-1),1)+...
            sparse(index_split_down_virtual.', ones(length(index_split_down_virtual),1), index_la- grid.n_y(ss-1)+(2:2:grid.n_y(ss-1)), grid.n_y(ss-1),1);
        index_node_hdown = [index_node_hddown_interface; index_node_hudown_interface; index_node_hdown];
        
        coor_node_h_y_d_up_interface = reshape(repmat((dy(ss):dy(ss):grid.bd_y),[2 1]),[ ],1);
        coor_node_h_y_d_down_interface = sparse(index_split_up_virtual.', ones(length(index_split_up_virtual),1), (dy(ss-1):2*dy(ss-1):grid.bd_y).' , grid.n_y(ss-1),1)+...
            sparse(index_split_down_virtual.', ones(length(index_split_down_virtual),1), (2*dy(ss-1):2*dy(ss-1):grid.bd_y).' , grid.n_y(ss-1),1);
        coor_node_h_y.down = [coor_node_h_y_d_down_interface; coor_node_h_y_d_up_interface; coor_node_h_y.down];
        
         %upnode
        index_uup = reshape(circshift(repmat((index_la+1:index_la+grid.n_y(ss)),[2 1]), 1,2),[ ],1);
        index_node_huup_interface = sparse(index_split_up_virtual.', ones(length(index_split_up_virtual),1), index_uup(index_split_up_virtual) ,  grid.n_y(ss-1),1)+...
            sparse(index_split_down_virtual.', ones(length(index_split_down_virtual),1), index_uup(index_split_down_virtual),grid.n_y(ss-1),1);
        index_node_hdup_interface = sparse(index_split_up_virtual.', ones(length(index_split_up_virtual),1), circshift(index_la - grid.n_y(ss-1)+(2:2:grid.n_y(ss-1)),1,2) ,grid.n_y(ss-1),1)+...
            sparse(index_split_down_virtual.', ones(length(index_split_down_virtual),1), index_la- grid.n_y(ss-1)+(1:2:grid.n_y(ss-1)),grid.n_y(ss-1),1);
        index_node_hup = [index_node_hdup_interface; index_node_huup_interface; index_node_hup];
        
        coor_node_h_y_u_down_interface = coor_node_h_y_d_down_interface - dy(ss-1);
        coor_node_h_y_u_up_interface = coor_node_h_y_d_up_interface -dy(ss);
        coor_node_h_y.up = [coor_node_h_y_u_down_interface; coor_node_h_y_u_up_interface; coor_node_h_y.up];
        
         %one node further up
        index_node_hup_interface = [index_node_hdup_interface.'; index_node_huup_interface.'];
        index_node_hupup_interface = reshape([circshift(index_node_hup_interface(1,:), 1, 2);circshift(index_node_hup_interface(2,:), 2, 2)].', [], 1);
        index_node_hupup = [index_node_hupup_interface; index_node_hupup];
        coor_node_hup_interface =  [coor_node_h_y_u_down_interface.'; coor_node_h_y_u_up_interface.'];
        coor_node_hupup_interface = reshape([coor_node_hup_interface(1,:) - dy(ss-1); coor_node_hup_interface(2,:) - dy(ss)].',[],1);
        coor_node_h_y.upup = [coor_node_hupup_interface; coor_node_h_y.upup];
        
        %one node further down
        index_node_hdown_interface = [index_node_hddown_interface.'; index_node_hudown_interface.'];
        index_node_hdowndown_interface = reshape([circshift(index_node_hdown_interface(1,:), -1, 2); circshift(index_node_hdown_interface(2,:), -2, 2)].', [], 1);
        index_node_hdowndown = [index_node_hdowndown_interface; index_node_hdowndown];
        coor_node_hdown_interface = [coor_node_h_y_d_down_interface.'; coor_node_h_y_d_up_interface.'];
        coor_node_hdowndown_interface = reshape([coor_node_hdown_interface(1,:) +  dy(ss-1); coor_node_hdown_interface(2,:) +  dy(ss)].', [],1);
        coor_node_h_y.downdown = [coor_node_hdowndown_interface; coor_node_h_y.downdown];
        
        %coordinates of virtual nodes
        y_virtual_interface = sort([3/4*dy(ss-1): 2*dy(ss-1):grid.bd_y,5/4*dy(ss-1): 2*dy(ss-1):grid.bd_y]).';
        coor_node_virtual.y = [y_virtual_interface; y_virtual_interface; coor_node_virtual.y];
        %% correction 4/1
        %OLD coor_node_virtual.z = [(z_top+ dz(ss-1))*ones(grid.n_y(ss-1),1);z_top*ones(grid.n_y(ss-1),1); coor_node_virtual.z]; 
        coor_node_virtual.z = [(z_top+ dz(ss-1)/2 + dz(ss)/2)*ones(grid.n_y(ss-1),1);z_top*ones(grid.n_y(ss-1),1); coor_node_virtual.z]; 
    end
    
    %must add all layers together
    fout.coor_nodes.y = [fout.coor_nodes.y; coor_nodes.y];
    fout.coor_nodes.z = [fout.coor_nodes.z; coor_nodes.z];
    fout.Delta_y_cell = [fout.Delta_y_cell; delta_y];
    fout.Delta_z_cell = [fout.Delta_z_cell; Delta_z_cell];
    fout.Delta_z_cell_volume = [fout.Delta_z_cell_volume; Delta_z_cell_volume];
    fout.delta_y_vflux = [fout.delta_y_vflux; delta_y_vflux ];
    fout.Delta_y_edge = [fout.Delta_y_edge; delta_y_edge ];
    
    fout.coor_hor_cell_edges.z = [fout.coor_hor_cell_edges.z; coor_hor_cell_edges.z];
    fout.coor_hor_cell_edges.y = [fout.coor_hor_cell_edges.y; coor_hor_cell_edges.y];
    
    fout.up_node.hor = [fout.up_node.hor; up_node.hor];
    fout.up_node.vert = [fout.up_node.vert; up_node.ver];
    fout.up_node_virtual.vert = [fout.up_node_virtual.vert; up_node_virtual.ver];
    fout.down_node.hor = [fout.down_node.hor; down_node.hor];
    fout.down_node.vert = [fout.down_node.vert; down_node.ver];
    fout.down_node_virtual.vert = [fout.down_node_virtual.vert; down_node_virtual.ver];
    fout.up_node.adv_ver = [fout.up_node.adv_ver; up_node_adv_ver_av];
    fout.down_node.adv_ver = [fout.down_node.adv_ver; down_node_adv_ver_av];
    fout.up_edge.hor = [fout.up_edge.hor; up_edge.hor];
    fout.down_edge.hor = [fout.down_edge.hor; down_edge.hor];
    
    fout.up_node.adv_hor = [fout.up_node.adv_hor; up_node_adv_hor];
    fout.down_node.adv_hor = [fout.down_node.adv_hor; down_node_adv_hor];
    
    fout.n_edges.vert = fout.n_edges.vert + size(coor_hor_cell_edges.y,1);
    
    fout.coor_node_h_y.down = [fout.coor_node_h_y.down; coor_node_h_y.down];
    fout.coor_node_h_y.downdown = [fout.coor_node_h_y.downdown; coor_node_h_y.downdown];
    fout.index_node_h.down = [fout.index_node_h.down; index_node_hdown];
    fout.index_node_h.up = [fout.index_node_h.up; index_node_hup];
    fout.index_node_h.downdown = [fout.index_node_h.downdown; index_node_hdowndown];
    fout.index_node_h.upup = [fout.index_node_h.upup; index_node_hupup];
    fout.coor_node_h_y.upup = [fout.coor_node_h_y.upup; coor_node_h_y.upup];
    fout.coor_node_h_y.up = [fout.coor_node_h_y.up; coor_node_h_y.up];   
    
    fout.coor_node_virtual.y = [fout.coor_node_virtual.y; coor_node_virtual.y];
    fout.coor_node_virtual.z = [fout.coor_node_virtual.z; coor_node_virtual.z];
end

fout.n_edges.vert = length(fout.up_node.vert);

% build connect matrix for vertical fluxes
connect_ver = sparse(fout.n_nodes.tot, fout.n_edges.vert);
for k = 1:fout.n_edges.vert
    index_up = fout.up_node.vert(k) ; 
    index_down = fout.down_node.vert(k); 
    connect_ver = connect_ver + sparse(index_up,k,1,fout.n_nodes.tot, fout.n_edges.vert);
    connect_ver = connect_ver + sparse(index_down,k,-1,fout.n_nodes.tot, fout.n_edges.vert);
end
fout.connect_ver = connect_ver;

fout.n_edges.hor = length(fout.up_node.hor);
% build connect matrix for horizontal fluxes
connect_hor = sparse(fout.n_nodes.tot, fout.n_edges.hor);
for k = 1:fout.n_edges.hor
    index_up = fout.up_node.hor(k) ; 
    index_down = fout.down_node.hor(k); 
    connect_hor = connect_hor + sparse(index_up,k,1,fout.n_nodes.tot, fout.n_edges.hor);
    connect_hor = connect_hor + sparse(index_down,k,-1,fout.n_nodes.tot, fout.n_edges.hor);
end
fout.connect_hor = connect_hor;



end

 