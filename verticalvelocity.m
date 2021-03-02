function [fout,dfout] = verticalvelocity(parameters,dpsi_dy,ddpsidy_dpsi)

% verticalvelocity.m takes as input the psi-dependent part of the vertical
% velocity and its jacobian and reindexes it on the horizontal cell edges
% of the T grid, as required by the temperature solver

psi_nodes = parameters.grid.psi.n_nodes.tot;                                   %number of nodes
T_up_node_ver = parameters.grid.N.up_node.vert;      

if parameters.flag1d == 0  && parameters.grid.N.extra.n_layers>1
    
    index_dpsidy_psitoT = sparse(length(T_up_node_ver),1);
    
    %set up the first layer; the interface always belongs to the layer below
    
    n_dpsidynodes_layer1 = parameters.grid.psi.n_nodes.hor(1).*(parameters.grid.psi.n_nodes.vert(1)-1);
    
    n_dpsidynodes_counter = n_dpsidynodes_layer1;
    n_Tvernetworkedges_counter = n_dpsidynodes_layer1;
    
    index_dpsidy_psitoT(1:n_dpsidynodes_layer1, 1) = 1:n_dpsidynodes_layer1;
    
    index_bulk_splitup = [];
    dpsi_dy_splitup_bulk = [];
    upnode_psi_splitup_bulk = [];
    downnode_psi_splitup_bulk = [];
    upupnode_psi_splitup_bulk = [];
    downdownnode_psi_splitup_bulk = [];
    index_bulk_splitdown = [];
    dpsi_dy_splitdown_bulk = [];
    upnode_psi_splitdown_bulk = [];
    downnode_psi_splitdown_bulk = [];
    upupnode_psi_splitdown_bulk = [];
    downdownnode_psi_splitdown_bulk = [];
    
    for jj= 2:length(parameters.grid.psi.n_nodes.vert)
        if jj<length(parameters.grid.psi.n_nodes.vert), ll =1; else ll=0;end
        n_dpsidynodes_layer = parameters.grid.psi.n_nodes.hor(jj).*(parameters.grid.psi.n_nodes.vert(jj)-ll) +parameters.grid.psi.n_nodes.hor(jj-1) ; %this includes the upper interface
        n_dpsidynodes_counter = n_dpsidynodes_counter + n_dpsidynodes_layer;
        
        n_Tvernetworkedges_layer = parameters.grid.N.n_nodes.hor(jj).*(parameters.grid.N.n_nodes.vert(jj)) + parameters.grid.N.n_nodes.hor(jj -1) ;
        n_Tvernetworkedges_counter = n_Tvernetworkedges_counter + n_Tvernetworkedges_layer;
        
        %deal with the interface first
        index_interface_single = n_Tvernetworkedges_counter - n_Tvernetworkedges_layer +[3:3:3*parameters.grid.N.n_nodes.hor(jj-1)];
        
        index_interface_splitup = n_Tvernetworkedges_counter - n_Tvernetworkedges_layer +[1:3:3*parameters.grid.N.n_nodes.hor(jj-1)];
        upnode_psi_split = n_dpsidynodes_counter - n_dpsidynodes_layer + circshift(1:parameters.grid.psi.n_nodes.hor(jj-1),1,2);
        downnode_psi_split = n_dpsidynodes_counter - n_dpsidynodes_layer + (1:parameters.grid.psi.n_nodes.hor(jj-1));
        upupnode_psi_split = n_dpsidynodes_counter - n_dpsidynodes_layer + circshift(1:parameters.grid.psi.n_nodes.hor(jj-1),2,2);
        downdownnode_psi_split = n_dpsidynodes_counter - n_dpsidynodes_layer + circshift(1:parameters.grid.psi.n_nodes.hor(jj-1),-1,2);
        
        dpsi_dy_splitup = -65/1024 * dpsi_dy(upupnode_psi_split) +715/1024*dpsi_dy(upnode_psi_split) + 429/1024*dpsi_dy(downnode_psi_split) -55/1024 * dpsi_dy(downdownnode_psi_split);
        
        index_interface_splitdown = n_Tvernetworkedges_counter - n_Tvernetworkedges_layer +[2:3:3*parameters.grid.N.n_nodes.hor(jj-1)];
        dpsi_dy_splitdown = -55/1024 * dpsi_dy(upupnode_psi_split) +429/1024*dpsi_dy(upnode_psi_split) + 715/1024*dpsi_dy(downnode_psi_split) - 65/1024 * dpsi_dy(downdownnode_psi_split);
        
        index_dpsidy_psitoT = index_dpsidy_psitoT + sparse(index_interface_single, ones(length(index_interface_single),1), n_dpsidynodes_counter - n_dpsidynodes_layer + (1:parameters.grid.psi.n_nodes.hor(jj-1)), ...
            length(T_up_node_ver),1);
        
        %now deal with the bulk
        index_bulk = (n_Tvernetworkedges_counter - n_Tvernetworkedges_layer +3*parameters.grid.N.n_nodes.hor(jj-1) +1):n_Tvernetworkedges_counter;
        index_psi_bulk = (n_dpsidynodes_counter - n_dpsidynodes_layer + parameters.grid.psi.n_nodes.hor(jj-1)+1):n_dpsidynodes_counter;
        index_dpsidy_psitoT = index_dpsidy_psitoT + sparse(index_bulk, ones(length(index_bulk),1), index_psi_bulk, length(T_up_node_ver),1);
        
        index_bulk_splitup = [index_bulk_splitup;(index_interface_splitup).' ];
        dpsi_dy_splitup_bulk = [dpsi_dy_splitup_bulk; dpsi_dy_splitup];
        upnode_psi_splitup_bulk = [upnode_psi_splitup_bulk; (upnode_psi_split).'];
        downnode_psi_splitup_bulk = [downnode_psi_splitup_bulk; (downnode_psi_split).'];
        upupnode_psi_splitup_bulk = [upupnode_psi_splitup_bulk; (upupnode_psi_split).'];
        downdownnode_psi_splitup_bulk = [downdownnode_psi_splitup_bulk; (downdownnode_psi_split).'];
        
        index_bulk_splitdown = [index_bulk_splitdown;(index_interface_splitdown).' ];
        dpsi_dy_splitdown_bulk = [dpsi_dy_splitdown_bulk; dpsi_dy_splitdown];
        upnode_psi_splitdown_bulk = [upnode_psi_splitdown_bulk; (upnode_psi_split).'];
        downnode_psi_splitdown_bulk = [downnode_psi_splitdown_bulk; (downnode_psi_split).' ];
        upupnode_psi_splitdown_bulk = [upupnode_psi_splitdown_bulk; (upupnode_psi_split).'];
        downdownnode_psi_splitdown_bulk = [downdownnode_psi_splitdown_bulk; (downdownnode_psi_split).'];
        
    end
    [index_bulk_tot,~] = find(index_dpsidy_psitoT~=0);
    connect_indexbulk = sparse(parameters.grid.N.n_edges.vert, parameters.grid.psi.n_edges.hor);
    for k = 1:length(index_bulk_tot)
        index = index_bulk_tot(k) ;
        connect_indexbulk = connect_indexbulk + sparse(index,k,1,parameters.grid.N.n_edges.vert, parameters.grid.psi.n_edges.hor);
    end
    
    dpsi_dy_reindexed = sparse(index_bulk_tot, ones(length(index_bulk_tot),1), dpsi_dy, length(T_up_node_ver),1)+ ...
        sparse(index_bulk_splitup, ones(length(index_bulk_splitup),1), dpsi_dy_splitup_bulk, length(T_up_node_ver),1)+ ...
        sparse(index_bulk_splitdown, ones(length(index_bulk_splitdown),1), dpsi_dy_splitdown_bulk, length(T_up_node_ver),1);
    fout = dpsi_dy_reindexed;
    
    %jacobian of dpsi_dy_reindexed
    %dpsi_dy_splitup = -65/1024 * dpsi_dy(upupnode_psi_split) +715/1024*dpsi_dy(upnode_psi_split) + 429/1024*dpsi_dy(downnode_psi_split) -55/1024 * dpsi_dy(downdownnode_psi_split);
    n_dpsidyTgrid = length(index_bulk_splitdown);
    ddpsidysplitupbulk_ddpsidy_up = 715/1024 * ones(n_dpsidyTgrid,1);
    ddpsidysplitupbulk_ddpsidy_upup = -65/1024 * ones(n_dpsidyTgrid,1);
    ddpsidysplitupbulk_ddpsidy_down = 429/1024 * ones(n_dpsidyTgrid,1);
    ddpsidysplitupbulk_ddpsidy_downdown = -55/1024 * ones(n_dpsidyTgrid,1);
    
    ddpsidysplitupbulk_dpsi = (sparse(index_bulk_splitup,upnode_psi_splitup_bulk,ddpsidysplitupbulk_ddpsidy_up, length(dpsi_dy_reindexed), psi_nodes)+...
        sparse(index_bulk_splitup,downnode_psi_splitup_bulk,ddpsidysplitupbulk_ddpsidy_down, length(dpsi_dy_reindexed), psi_nodes) + ...
        sparse(index_bulk_splitup,upupnode_psi_splitup_bulk,ddpsidysplitupbulk_ddpsidy_upup, length(dpsi_dy_reindexed), psi_nodes)+...
        sparse(index_bulk_splitup,downdownnode_psi_splitup_bulk,ddpsidysplitupbulk_ddpsidy_downdown, length(dpsi_dy_reindexed), psi_nodes))*ddpsidy_dpsi;
    
    %dpsi_dy_splitdown = -55/1024 * dpsi_dy(upupnode_psi_split) +429/1024*dpsi_dy(upnode_psi_split) + 715/1024*dpsi_dy(downnode_psi_split) - 65/1024 * dpsi_dy(downdownnode_psi_split);
    ddpsidysplitdownbulk_ddpsidy_up = 429/1024 * ones(n_dpsidyTgrid,1);
    ddpsidysplitdownbulk_ddpsidy_upup = -55/1024 * ones(n_dpsidyTgrid,1);
    ddpsidysplitdownbulk_ddpsidy_down = 715/1024 * ones(n_dpsidyTgrid,1);
    ddpsidysplitdownbulk_ddpsidy_downdown = -65/1024 * ones(n_dpsidyTgrid,1);
    
    ddpsidysplitdownbulk_dpsi = (sparse(index_bulk_splitdown,upnode_psi_splitdown_bulk,ddpsidysplitdownbulk_ddpsidy_up, length(dpsi_dy_reindexed), psi_nodes)+...
        sparse(index_bulk_splitdown,downnode_psi_splitdown_bulk,ddpsidysplitdownbulk_ddpsidy_down, length(dpsi_dy_reindexed), psi_nodes) + ...
        sparse(index_bulk_splitdown,upupnode_psi_splitdown_bulk,ddpsidysplitdownbulk_ddpsidy_upup, length(dpsi_dy_reindexed), psi_nodes)+...
        sparse(index_bulk_splitdown,downdownnode_psi_splitdown_bulk,ddpsidysplitdownbulk_ddpsidy_downdown, length(dpsi_dy_reindexed), psi_nodes))*ddpsidy_dpsi;
    
    %add the jacobian of the regular nodes
    Ddpsidy_reindexed = connect_indexbulk * ddpsidy_dpsi; 
    
    %sum all three pieces together
    dfout = ddpsidysplitupbulk_dpsi + ddpsidysplitdownbulk_dpsi + Ddpsidy_reindexed;
    
elseif parameters.flag1d == 1 
    fout = zeros(parameters.grid.N.n_edges.vert,1);
    dfout = zeros(parameters.grid.N.n_edges.vert,psi_nodes);
    
elseif  parameters.grid.N.extra.n_layers==1
    
    index_dpsidy_psitoT = sparse(length(T_up_node_ver),1);
    
    %set up the first layer; the interface always belongs to the layer below
    
    n_dpsidynodes_layer1 = parameters.grid.psi.n_nodes.hor(1).*(parameters.grid.psi.n_nodes.vert(1));
    
    index_dpsidy_psitoT(1:n_dpsidynodes_layer1, 1) = 1:n_dpsidynodes_layer1;
    
    [index_bulk_tot,~] = find(index_dpsidy_psitoT~=0);
      connect_indexbulk = sparse(parameters.grid.N.n_edges.vert, parameters.grid.psi.n_edges.hor);
    for k = 1:length(index_bulk_tot)
        index = index_bulk_tot(k) ;
        connect_indexbulk = connect_indexbulk + sparse(index,k,1,parameters.grid.N.n_edges.vert, parameters.grid.psi.n_edges.hor);
    end
   
    
    dpsi_dy_reindexed = sparse(index_bulk_tot, ones(length(index_bulk_tot),1), dpsi_dy, length(T_up_node_ver),1);
    fout = dpsi_dy_reindexed;
    Ddpsidy_reindexed = connect_indexbulk * ddpsidy_dpsi; 
    dfout = Ddpsidy_reindexed;
end

end