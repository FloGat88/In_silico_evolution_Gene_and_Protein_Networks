classdef Network < handle
    properties
        subnetworks Subnetwork
        N_elems int16
        n_sub(1,1) int16
        reactions struct
        decays struct
        
        %% for Ode Solver
        states int8
        initial_elements int16
        initial_densities double
        OdeTransitions int16
        OdeTransitions_rates double
        OdeReactions int16
        OdeReactions_rates double
        OdeDecays int16
        OdeDecays_rates double
        composites int8     % keep this only temporarily in order to check for mass conservation
        
        % for Fitness
        Fitness double
%        Fitness_max double
    end
    
    methods
        function obj = Network()
        %    obj.subnetworks = Subnetwork(false,Element(true,1,rand(),1),-1,-1);
            obj.N_elems = 0;
            obj.n_sub = 0;
            obj.reactions = struct('rate',{},'from1',{},'from2',{},'to',{});
            obj.decays = struct('rate',{},'from',{},'to1',{},'to2',{});
            % create some initial structure
            obj.add_atomic_subnetwork();
            obj.add_derived_element();
            obj.add_derived_element();
            obj.add_atomic_subnetwork();
            obj.add_derived_element();
            obj.add_composite_subnetwork();
            obj.add_derived_element();
            
            obj.Fitness = 0;
 %           obj.Fitness_max = 0;
        end
        
        function N = postincrN(obj)     % post increment operation for N_elems (number of elements in the network)
            N = obj.N_elems;
            obj.N_elems = obj.N_elems + 1;
        end
        
        function ind = choose_subnetwork(obj,attribute)
            if nargin==1
                attribute = '-';
            end
            if strcmp(attribute,'atomic')
                Y = [obj.subnetworks.existent] & [obj.subnetworks.composite]==false;
            elseif strcmp(attribute,'composite')
                Y = [obj.subnetworks.existent] & [obj.subnetworks.composite];
            else
                Y = [obj.subnetworks.existent];
            end
            ind = find_random(Y);
        end  
        
        function X = choose_element(obj,attribute,state)
            if nargin<3
                state = 0;
            end
            if nargin<2
                attribute = '-';
            end
            ind1 = obj.choose_subnetwork(attribute);
            if isempty(ind1)
                X = []; return;
            end
            ind2 = obj.subnetworks(ind1).choose_element(1,state,'non-different');
            if isempty(ind2)
                X = []; return;
            end
            X = [ind1;ind2];
        end
        
        function success = add_transition(obj)
            ind = find_random([obj.subnetworks.existent] & ([obj.subnetworks.n_state1]+[obj.subnetworks.n_state2])>=2);
            if isempty(ind)
                success = 0; return;
            end
            % add transition on subnetwork level and update OdeSystem accordingly
            [success,rate,ode_ind,existent] = obj.subnetworks(ind).add_transition();
            obj.add_OdeTransition(rate,ode_ind,existent);
        end
            
        function success = add_decay(obj,X_from)
            % if no element index is provided choose one randomly
            if nargin==1
                X_from = obj.choose_element('composite');
                if isempty(X_from)
                    success = 0; return;
                end
            end
            % determine state and constituents of decaying element
            state = obj.subnetworks(X_from(1)).elements(X_from(2)).state;
            constituent1 = obj.subnetworks(X_from(1)).constituent1;
            constituent2 = obj.subnetworks(X_from(1)).constituent2;
            % enforce state rules for the state 'to_state' of the decay products 
            if state==1
                to_state = 1;
            else
                to_state = 0;
            end
            % choose two indices within the respective subnetworks obeying the state rules
            ind1 = obj.subnetworks(constituent1).choose_element(1,to_state,'non-different');
            ind2 = obj.subnetworks(constituent2).choose_element(1,to_state,'non-different');
            if isempty(ind1) || isempty(ind2)
                success = 0; return;
            end
            % add decay 
            obj.decays(end+1) = struct('rate',rand(),'from',X_from,'to1',[constituent1;ind1],'to2',[constituent2;ind2]);
            % increment the respective connectivity indices
            obj.subnetworks(constituent1).elements(ind1).increment_connectivity();
            obj.subnetworks(constituent2).elements(ind2).increment_connectivity();  
            % update OdeSystem accordingly
            obj.OdeDecays_rates(end+1) = obj.decays(end).rate;
            obj.OdeDecays(:,end+1) = [obj.net2ode(obj.decays(end).from);obj.net2ode(obj.decays(end).to1);obj.net2ode(obj.decays(end).to2)];
            
            success = 1;
        end
        
        function success = add_composite_subnetwork(obj)
            if obj.n_sub==0
                success = 0; return;
            end
            % choose two random atomic, i.e. non-composite elements that form the components of the new complex. 
            X1 = obj.choose_element('atomic',0);
            X2 = obj.choose_element('atomic',0);
            % determine state of the new complex in accordance with the state-rules
            if obj.subnetworks(X1(1)).elements(X1(2)).state==1 && obj.subnetworks(X2(1)).elements(X2(2)).state==1
                state = 1;
            else
                state = randi(2,1,1);
            end
            % create the new subnetwork either at an emty position of 'subnetworks' or by appending it to the vector of subnetworks 
            ind = find([obj.subnetworks.existent]==false,1);
            if isempty(ind)
                ind = length(obj.subnetworks) + 1;
            end
            obj.subnetworks(ind) = Subnetwork(true,Element(true,state,0,obj.postincrN()),X1(1),X2(1));
            obj.add_state(state);
            % increment the n_sub index
            obj.n_sub = obj.n_sub+1;
            % add a new reaction that connects X1 and X2 to the new element
            obj.reactions(end+1) = struct('rate',rand(),'from1',X1,'from2',X2,'to',[ind;1]);
                % update OdeSystem accordingly
                obj.OdeReactions_rates(end+1) = obj.reactions(end).rate;
                obj.OdeReactions(:,end+1) = [obj.net2ode(X1),obj.net2ode(X2),obj.net2ode([ind;1])];
            % maybe also add a new decay reaction
            obj.add_decay([ind;1]);   
            success = 1;
        end
        
        function add_atomic_subnetwork(obj)
            % create the new subnetwork either at an emty position of 'subnetworks' or by appending it to the vector of subnetworks 
            ind = find([obj.subnetworks.existent]==false,1);
            if isempty(ind)
                ind = length(obj.subnetworks) + 1;
            end
            initial_density = rand()*500;       % choose 600 as maximum initial density... see previous versions
            ode_index = obj.postincrN();
            obj.subnetworks(ind) = Subnetwork(false,Element(true,1,initial_density,ode_index),-1,-1);
            obj.add_state(1);
            obj.add_initial_element(ode_index,initial_density);
            % increment the n_sub index
            obj.n_sub = obj.n_sub + 1;
        end
        
        function success = add_derived_element(obj)
            if obj.n_sub==0
                success=0; return;
            end
            % choose a subnetwork and add a derived element for the chosen subnetwork
            ind1 = obj.choose_subnetwork('-');
            % add derived element and update OdeTransitions accordingly
            [ind2,state,rates,ode_inds,existent] = obj.subnetworks(ind1).add_derived_element(obj.postincrN());
            obj.add_OdeTransition(rates,ode_inds,existent);
            obj.add_state(state);
            % if an element has been added to a composite subnetwork also add a decay of the new element (if possible)
            if obj.subnetworks(ind1).composite==true
                obj.add_decay([ind1;ind2]);
            end
            success = 1;
        end
        
        function delete_reaction(obj,ind)
            obj.reactions(ind) = obj.reactions(end);      
            obj.reactions(end) = [];
        end
        
        function delete_decay(obj,ind)
            ind_elem1 = obj.decays(ind).to1;
            ind_elem2 = obj.decays(ind).to2;
            obj.subnetworks(ind_elem1(1)).elements(ind_elem1(2)).decrement_connectivity();
            obj.subnetworks(ind_elem2(1)).elements(ind_elem2(2)).decrement_connectivity();
            obj.decays(ind) = obj.decays(end);      
            obj.decays(end) = [];
        end
        
        function success = delete_composite_subnetwork(obj,ind)   % delete composite subnetwork and delete all reactions and decays 'downward'
            if nargin == 1
                ind = obj.choose_subnetwork('composite');
                if isempty(ind)
                    success = 0; return;
                end
            end     
            % delete the respective subnetwork by setting existent to false and decrement n_sub 
            obj.subnetworks(ind).existent = false;
            obj.n_sub = obj.n_sub - 1;
            % find the indices of all reactions that go to the composite subnetwork
            ind_reaction = int16(find( colcomp([obj.reactions.to],ind,1) ));
            ind_reaction = fliplr(ind_reaction);    % sort in descending order so that subsequent deletions in the following loop do not intersect with each other
            for i=1:length(ind_reaction)
                obj.delete_reaction(ind_reaction(i));
            end
            % find the indices of all decays that emerge from the composite subnetwork
            ind_decay = int16( find(colcomp([obj.decays.from],ind,1)) );
            ind_decay = fliplr(ind_decay);       % sort in descending order so that subsequent deletions in the following loop do not intersect with each other
            for i=1:length(ind_decay)
                obj.delete_decay(ind_decay(i));
            end
            % if higher order complexes like trimers, quadromers etc. can
            % exist then you check here for builder elements that are composed of this subnetwork and delete the builder subnetworks as well (cf. 'delete_atomic_subnetwork')  
            success = 1;
        end
        
        function success = delete_atomic_subnetwork(obj)   % delete atomic subnetwork and all building subnetworks (which then destroy the connecting reactions and decays)
            if obj.n_sub==0
                success = 0; return;
            end
            % find random existing non-composite (= atomic) subnetwork
            ind = find_random([obj.subnetworks.existent] & [obj.subnetworks.composite]==false);
            % delete the respective subnetwork by setting existent to false and decrement n_sub 
            obj.subnetworks(ind).existent = false;
            obj.n_sub = obj.n_sub - 1;
            % find the indices of all building subnetworks, i.e. those that are composed of the deleted subnetwork and delete them as well
            % all respective connecting reactions and decays are deleted from the building subnetworks by calling 'delete_composite_subnetwork'
            ind_building_subnetworks = find([obj.subnetworks.existent] & ([obj.subnetworks.constituent1]==ind | [obj.subnetworks.constituent2]==ind));
            for i=1:length(ind_building_subnetworks)
                obj.delete_composite_subnetwork(ind_building_subnetworks(i));
            end
            success = 1;
    %{
            % find the indices of all reactions that emerge from the deleted subnetwork
            from1 = [obj.reactions.from1];
            from2 = [obj.reactions.from2];
            ind_reaction = find(from1(1,:)==ind | from2(1,:)==ind);
            % extract the indices of all building subnetworks, i.e. those that are composed of the deleted subnetwork and hence must be deleted as well
            ind_building_subnetworks = obj.reactions(ind_reactions).to;
            ind_building_networks = ind_building_networks(1,:);
            % delete all reactions that emerge from the deleted subnetwork
            for i=1:length(ind_reaction)
                obj.delete_reaction(ind_reaction);
            end
            % find the indices of all decays that connect to the deleted subnetwork
            to = [obj.reactions.from1];
            from2 = [obj.reactions.from2];
            ind_reaction = find(from1(1,:)==ind | from2(1,:)==ind);
            %}
        end
        
        function success = delete_derived_element(obj)
            ind1 = find_random([obj.subnetworks.existent] & ([obj.subnetworks.n_state1]+[obj.subnetworks.n_state2])>1);
            if isempty(ind1)
                success = 0; return;
            end
            ind2 = obj.subnetworks(ind1).delete_derived_element_cascade();
            if obj.subnetworks(ind1).composite==false    % in case the chosen subnetwork is atomic
                % find the indices of all reactions that go from the subnetwork in which the cascade is triggered; 
                % therefore, first create the matrix that contains in its second and third row all the
                % indices of the source elements of all reactions and in its first row contains the indices of the corresponding reactions. 
                if ~isempty(obj.reactions)
                    X_from = [repmat(1:length(obj.reactions),1,2);[obj.reactions.from1],[obj.reactions.from2]];
                    X_from(:,X_from(2,:)~=ind1) = [];     
                    X_from(:,[obj.subnetworks(ind1).elements(X_from(3,:)).existent]) = [];
                    X_to = [obj.reactions(X_from(1,:)).to];
                    if ~isempty(X_to)
                        ind_reactions = sort(X_to(1,:),'descend');      % sort in descending order so that subsequent deletions in the following loop do not intersect with each other 
                        for i=1:length(ind_reactions)
                            obj.delete_composite_subnetwork(ind_reactions(i));
                        end
                    end
                end
                if ~isempty(obj.decays)
                    ind_decays = find( colcomp([obj.decays.to1],[ind1;ind2]) | colcomp([obj.decays.to2],[ind1;ind2]));
                    ind_decays = fliplr(ind_decays);      % sort in descending order so that subsequent deletions in the following loop do not intersect with each other 
                    for i=1:length(ind_decays)
                        obj.delete_decay(ind_decays(i));
                    end
                end
            else    % in case the chosen subnetwork is composite
                % no need to modify reactions in this case because only derived elements can be deleted that are not connected via reactions
                if ~isempty(obj.decays)
                    X_from = [1:length(obj.decays);obj.decays.from];
                    X_from(:,X_from(2,:)~=ind1) = [];     
                    X_from(:,[obj.subnetworks(ind1).elements(X_from(3,:)).existent]) = [];
                    if ~isempty(X_from)
                        ind_decays = sort(X_from(1,:),'descend');
                        for i=1:length(ind_decays)
                            obj.delete_decay(ind_decays(i));
                        end
                    end
                end
            end  
            success = 1;
        end
        
        function success = mutate_init_dens(obj)
            if obj.n_sub==0
                success = 0; return;
            end
            % find some atomic (i.e. non-composite) subnetwork
            ind = obj.choose_subnetwork('atomic');
            % within each subnetwork the base_element is at position 1 in the elements vector; multiply its init_density with a random number between 0 and 2
            obj.subnetworks(ind).elements(1).init_density = obj.subnetworks(ind).elements(1).init_density * 2*rand();
            % update initial_elements and initial_densities accordingly 
            init_el_ind = find([obj.initial_elements]==obj.subnetworks(ind).elements(1).ode_index,1);
            obj.initial_densities(init_el_ind) = obj.subnetworks(ind).elements(1).init_density;
            success = 1;
        end
        
        function success = mutate_transition_rate(obj)
            if obj.n_sub==0
               success = 0; return;
            else
                ind = obj.choose_subnetwork();
                [success,rates] = obj.subnetworks(ind).mutate_transition_rate();
                % update Ode System accordingly; find the modified transition by searching OdeTransitions for the old rate; alternative way to find the corresponding transition:     obj.OdeTransitions_rates(find(obj.OdeTransitions(1,:)==net2ind(transition.from) & obj.OdeTransitions(2,:)==net2ind(transition.to),1)) = transition.rate;
                obj.OdeTransitions_rates(obj.OdeTransitions_rates==rates(1)) = rates(2);
            end      
        end
        
        function success = mutate_reaction_rate(obj)
            if isempty(obj.reactions)
               success = 0; 
            else
                ind = randi(length(obj.reactions),1,1);
                obj.reactions(ind).rate = obj.reactions(ind).rate * 2*rand();
                success = 1;
                % update Ode System accordingly
                obj.OdeReactions_rates(ind) = obj.reactions(ind).rate;
            end      
        end
        
        function success = mutate_decay_rate(obj)
            if isempty(obj.decays)
               success = 0; return;
            else
                ind = randi(length(obj.decays),1,1);
                obj.decays(ind).rate = obj.decays(ind).rate * 2*rand();
                success = 1;
                % update Ode System accordingly
                obj.OdeDecays_rates(ind) = obj.decays(ind).rate;
            end      
        end
        
        % mutate a parameter of the network
        function success = mutate_rate(obj)      
            r = randi(5,1,1);
            if r==1
                success = obj.mutate_init_dens();
            elseif r<=3
                success = obj.mutate_transition_rate();
            elseif r==4
                success = obj.mutate_reaction_rate();
            elseif r==5
                success = obj.mutate_decay_rate(); 
            end
        end
        
        % mutate the structure (topology) of the network by adding nodes (and connections that connect these new nodes)
        function success = mutate_add_element(obj)
            r = randi(7,1,1);
            if r<=4
                success = obj.add_derived_element();    % corresponds to 'add_phosphorylate' in Version 2 of the algorithm
            elseif r<=6
                success = obj.add_composite_subnetwork;    % corresponds to 'add_binding_site' in Version 2 of the algorithm
            elseif r==7
                obj.add_atomic_subnetwork();    % corresponds to 'add_protein' in Version 2 of the algorithm
                success = 1;        % always successful
            end
        end
        
        % mutate the structure (topology) of the network by adding additional connections between existing nodes
        function success = mutate_add_connection(obj)
            r = randi(4,1,1);
            if r<=3
                success = obj.add_transition();
            elseif r==4
                success = obj.add_decay();
            end
        end
        
        % mutate the structure (topology) of the network by removing existing nodes and their connections
        function success = mutate_delete_element(obj)
            r = randi(7,1,1);
            if r<=4
                success = obj.delete_derived_element();    % corresponds to 'delete_phosphorylate' in Version 2 of the algorithm
            elseif r<=6
                success = obj.delete_composite_subnetwork;    % corresponds to 'delete_binding_site' in Version 2 of the algorithm
            elseif r==7
                success = obj.delete_atomic_subnetwork();    % corresponds to 'delete_protein' in Version 2 of the algorithm
            end
        end
        
        function Mutate(obj,times)
            OdeRebuild = false;
            for i=1:times
                success = 0;
                while(success==0)
                    r = randi(14,1,1);
                    if r<=6
                        success = obj.mutate_rate();
                    elseif r<=9
                        success = obj.mutate_add_connection();
                    elseif r<=12
                        success = obj.mutate_add_element();
                    else   % r<=14    % if mutate_delete_element is called activate OdeRebuild to rebuild the OdeSystem below
                        success = obj.mutate_delete_element();
                        if success==true;  OdeRebuild=true;  end
                    end
                end
            end
            if OdeRebuild
                obj.OdeRebuild();
            end 
        end
  
        
        
 %% Processing for ode solver
 
 function determine_N_elems(obj)
     obj.N_elems = sum([obj.subnetworks([obj.subnetworks.existent]).n_state1] + [obj.subnetworks([obj.subnetworks.existent]).n_state2]); 
 end
 function assign_ode_index_make_states(obj)
     index_end=0;
     obj.states = zeros(obj.N_elems,1);
     for subnet=obj.subnetworks([obj.subnetworks.existent])
         index_begin = index_end+1;
         index_end = index_end + subnet.get_numelems();
         obj.states(index_begin:index_end) = subnet.assign_ode_index_process_states(index_begin-1,index_end-1);    % subtract 1 because the ode_indices should be c++ indices that are 0-based
     end
 end
 function make_initial_elements(obj)
    num_initial_el = sum([obj.subnetworks.existent] & ~[obj.subnetworks.composite]);
    obj.initial_elements = int16(zeros(1,num_initial_el));
    obj.initial_densities = zeros(1,num_initial_el);
    i=1;
    for subnet = obj.subnetworks([obj.subnetworks.existent] & ~[obj.subnetworks.composite])
        obj.initial_elements(i) = subnet.elements(1).ode_index;
        obj.initial_densities(i) = subnet.elements(1).init_density;
        i=i+1;
    end
 end
 %%
 function make_composites(obj)
     obj.composites = ones(obj.N_elems,1);
     for subnet = obj.subnetworks([obj.subnetworks.existent] & [obj.subnetworks.composite])
         obj.composites([subnet.elements([subnet.elements.existent]).ode_index]+1) = 2;
     end
 end
 %%
 function ode_vec = net2ode(obj,net_array)
     ode_vec = zeros(1,size(net_array,2));
     for i=1:size(net_array,2)
         ode_vec(i) = obj.subnetworks(net_array(1,i)).elements(net_array(2,i)).ode_index;
     end
 end
 function add_state(obj,state)
     obj.states(end+1) = state;
 end
 function add_initial_element(obj,ode_index,initial_density)
     obj.initial_elements(end+1) = ode_index;
     obj.initial_densities(end+1) = initial_density;
 end
 function add_OdeTransition(obj,rates,ode_inds,existent)
     for i=1:length(existent)
         if existent(i)==true   % if the transition already exists add rate to the existing rate otherwise add a new transition to OdeTransitions  
             ind = find(obj.OdeTransitions(1,:)==ode_inds(1,i) & obj.OdeTransitions(2,:)==ode_inds(2,i),1);
             obj.OdeTransitions_rates(ind) = obj.OdeTransitions_rates(ind) + rates(i);
         else
             obj.OdeTransitions_rates(:,end+1) = rates(i);
             obj.OdeTransitions(:,end+1) = ode_inds(:,i);
         end
     end
 end
 function s = get_numberTransitions(obj)
     s = 0;
     for subnet = obj.subnetworks([obj.subnetworks.existent])
         s = s + length(subnet.transitions);
     end
 end
 function make_OdeTransitions(obj)
     size = obj.get_numberTransitions();
     obj.OdeTransitions = zeros(2,size);
     obj.OdeTransitions_rates = zeros(1,size);
     ind1=int16(0);
     ind2=int16(0);
     for subnet=obj.subnetworks([obj.subnetworks.existent])
         ind1 = ind2+1;
         ind2 = ind2 + length(subnet.transitions());
         [obj.OdeTransitions_rates(ind1:ind2),obj.OdeTransitions(:,ind1:ind2)] = subnet.process_transitions();
     end
 end
 function make_OdeReactions(obj)
     obj.OdeReactions_rates = [obj.reactions.rate];
     obj.OdeReactions = [obj.net2ode([obj.reactions.from1]);obj.net2ode([obj.reactions.from2]);obj.net2ode([obj.reactions.to])];
 end
 function make_OdeDecays(obj)
     obj.OdeDecays_rates = [obj.decays.rate];
     obj.OdeDecays = [obj.net2ode([obj.decays.from]);obj.net2ode([obj.decays.to1]);obj.net2ode([obj.decays.to2])];
 end
 function OdeRebuild(obj)
     obj.determine_N_elems();
     obj.assign_ode_index_make_states();
     obj.make_initial_elements();
     obj.make_OdeTransitions();
     obj.make_OdeReactions();
     obj.make_OdeDecays();
 end
 
 
 %% Fitness calculation
 
 function Fitness = FitnessPolarisation(obj)    % loss function to evolve a polarized state in which one species is inhomogeneously distributed between right and left
    D_cyt = 10;      % cytosolic diffusion constant
    D_mem = 0.1;    % membrane diffusion constant
    L = 10;     % system size
    N_L = int16(20);   % number of discretisation points
    a = L/double(N_L);  % lattice spacing
    T_max = 100;        % in the c++ version: T_max=20
    % T_max = L^2/D_mem;    % maximum simulation time; choose it such that the the system can approximately be transversed once by the slowly diffusing species during the simulation
    D_eff = zeros(1,obj.N_elems);         % effective diffusion vector
    D_eff(obj.states==1) = D_cyt/a^2;
    D_eff(obj.states==2) = D_mem/a^2;
    %% Solution with finite differences
    x0 = zeros(obj.N_elems*N_L,1);      % vector of initial concentrations
    for i=1:length(obj.initial_elements)    % assign initial concentrations and add some noise
        x0(obj.initial_elements(i)+1:obj.N_elems:length(x0)) = obj.initial_densities(i)*(1+0.1*randn(N_L,1));   % note that ode_indices in initial_elements start counting at 0 (according to c++ convention) however for matlab we need counting starting from 1. Therefore add 1 here to initial element
    end
    % Simulation
    dynamics_handle = @(t,x) reaction_diffusion_dynamics(t,x,obj.OdeTransitions_rates,obj.OdeReactions_rates,obj.OdeDecays_rates,obj.OdeTransitions,obj.OdeReactions,obj.OdeDecays,D_eff,[obj.N_elems,N_L]);
    Jacobian_handle = @(t,x) Jacobian(t,x,obj.OdeTransitions_rates,obj.OdeReactions_rates,obj.OdeDecays_rates,double(obj.OdeTransitions+1),double(obj.OdeReactions+1),double(obj.OdeDecays+1),D_eff,double(obj.N_elems),N_L);       % add 1 to all indices to convert from c++ indices (0-based) to matlab indices (1-based); furthermore convert indices to 'double' because for some strange reason the 'sparse' and 'spdiags' commands in 'Jacobian' only accept inputs as doubles
 %   for i=1:10000
 %       dynamics_handle(0,x0);
 %   end    
    [~,X] = ode15s(dynamics_handle, [0 T_max], x0, odeset('RelTol',1e-5,'AbsTol',1e-5,'NonNegative',(1:length(x0)),'Jacobian',Jacobian_handle,'Stats','on' ));  %,'Jacobian',@DavidModel_jacobian));   % ,'JPattern',S 'Jacobian',@Winter_Wonderland_single_noise_jacobian  ,'Jacobian',Jacobian    ,'Stats','on' ,'OutputFcn'
%    X = (x0 + dynamics_handle(0,x0))';
    obj.make_composites();      % remove these lines again... only to check if mass is conserved!
    obj.check_mass_conservation(X(end,:),sum(x0));
    % Fitness evaluation
    probe_vec = [-1*ones(N_L/2,1);ones(N_L/2,1)];
    C = reshape(X(end,:),obj.N_elems,N_L);
    polarization_vec = (C * probe_vec)./vecnorm(C,2,2);
    Fitness = 100 * max(abs(polarization_vec));
    if Fitness<=2       % include a fluctuation boundary. Maybe the boundary must be chosen even larger
        Fitness=0;
    end
    
 end

 function Fitness = FitnessPolarisation_pdepe(obj)     % Simulation with finite elements via pdepe
    if obj.N_elems==0
        Fitness = 0; return;
    end
    D_cyt = 10;      % cytosolic diffusion constant
    D_mem = 0.1;    % membrane diffusion constant
    L = 10;     % system size
    T_max = 100;        % in the c++ version: T_max=20
    % T_max = L^2/D_mem;    % maximum simulation time; choose it such that the the system can approximately be transversed once by the slowly diffusing species during the simulation
    D_mat = zeros(1,obj.N_elems);         % effective diffusion vector
    D_mat(obj.states==1) = D_cyt;
    D_mat(obj.states==2) = D_mem;
    % Simulation
%    return; %!!!!!!
    mesh = linspace(0,L,20);
    N_L = length(mesh);
    t = linspace(0,T_max,10);
    x0 = zeros(N_L,obj.N_elems);      % vector of initial concentrations
    for i=1:length(obj.initial_elements)    % assign initial concentrations and add some noise
        x0(:,obj.initial_elements(i)+1) = obj.initial_densities(i)*(1+0.1*randn(N_L,1));   % note that ode_indices in initial_elements start counting at 0 (according to c++ convention) however for matlab we need counting starting from 1. Therefore add 1 here to initial element
    end
    ic_handle = @(x)interp1(mesh,x0,x)';
    bc_handle = @(xl,ul,xr,ur,t) bc_fun_pdepe(xl,ul,xr,ur,t,obj.N_elems);
    dynamics_handle = @(x,t,u,dudx) reaction_diffusion_dynamics_pdepe(x,t,u,dudx,obj.OdeTransitions_rates,obj.OdeReactions_rates,obj.OdeDecays_rates,obj.OdeTransitions,obj.OdeReactions,obj.OdeDecays,D_mat,obj.N_elems);
    sol = pdepe(0,dynamics_handle,ic_handle,bc_handle,mesh,t);
    X(:,:) = sol(end,:,:); 
    X = reshape(X',1,obj.N_elems*N_L);
    % check mass conservation
    obj.make_composites();      % remove these lines again... only to check if mass is conserved!
    obj.check_mass_conservation(X,sum(sum(x0)))
    % Fitness evaluation
    probe_vec = [-1*ones(N_L/2,1);ones(N_L/2,1)];
    C = reshape(X(end,:),obj.N_elems,N_L);
    norm = vecnorm(C,1,2);
    norm(norm<1) = 0;    % eliminate those components that fluctuate around zero (and usually also have some negative values) and therefore can have a large relative polarization only due to the initial fluctuations  
    norm(norm~=0) = 1./norm(norm~=0);       % inverse all entries that are non-zero.
    polarization_vec = (C * probe_vec).*norm;
    Fitness = 100 * max(abs(polarization_vec));
    if Fitness<=2       % include a fluctuation boundary. Maybe the boundary must be chosen even larger
        Fitness=0;
    end
    obj.Fitness = Fitness;
%    if Fitness > obj.Fitness_max
%        obj.Fitness_max = Fitness;
%    end
 end
 
 
  %% Test classes
        function Test_Decays_Reactions_existent(obj)
            if isempty(obj.decays)
                return;
            end
            % test decays and reactions if all sources and sinks are indeed existent
            nodes = [[obj.decays.from],[obj.decays.to1],[obj.decays.to2],[obj.reactions.from1],[obj.reactions.from2],[obj.reactions.to]];
            for i=1:size(nodes,2)
                if obj.subnetworks(nodes(1,i)).existent == false || obj.subnetworks(nodes(1,i)).elements(nodes(2,i)).existent == false
                    dbstop;
                end
            end
 %           composite = [obj.subnetworks(nodes(1,:)).composite];
 %           if ~all(composite(1:length(obj.decays))) || any(composite(length(obj.decays)+1:end))
 %               dbstop;
 %           end
        end
        
        function relative_mass_conservation = check_mass_conservation(obj,X,initial_mass)
            N_L = length(X)/obj.N_elems;
            mass_have = X * kron(ones(N_L,1),double(obj.composites));
            relative_mass_conservation = (mass_have - initial_mass)/initial_mass;
        end
    end
end

