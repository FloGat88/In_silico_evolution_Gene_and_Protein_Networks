classdef Subnetwork < handle
    properties
        existent(1,1) logical
        composite(1,1) logical
        elements Element
        n_state1(1,1) int16
        n_state2(1,1) int16
        constituent1(1,1) int16
        constituent2(1,1) int16
        transitions struct
    end
    
    methods
        function obj = Subnetwork(composite_,element_,constit1,constit2)
            obj.existent = true;
            obj.composite = composite_;
            obj.elements = element_;
            obj.n_state1 = int8(obj.elements.state==1);
            obj.n_state2 = int8(obj.elements.state==2);
            obj.constituent1 = constit1;
            obj.constituent2 = constit2;
            obj.transitions = struct('rate',{},'from',{},'to',{});
        end
        
        function increment_n_state(obj,state)
            if state==1
                obj.n_state1 = obj.n_state1 + 1;
            else
                obj.n_state2 = obj.n_state2 + 1;
            end
        end
        
        function decrement_n_state(obj,state)
            if state==1
                obj.n_state1 = obj.n_state1 - 1;
            else
                obj.n_state2 = obj.n_state2 - 1;
            end
        end
        
        function N = get_numelems(obj)
            N = obj.n_state1 + obj.n_state2;
        end
        
        function X = choose_element(obj,n,state,attr) 
             n1 = obj.n_state1;     % create copies of the nubmers of posphs of certain state to be able to maipulate these in the following
             n2 = obj.n_state2;
             if state==1    % if only cytosolic phophs desired set n2 to 0
                 n2 = 0;
             elseif state==2    % if only membrane-bound phophs desired
                 n1 = 0;
             end
             if strcmp(attr,'different')
                 if n1+n2<n
                     X = []; return;
                 end
                 indices = randperm(n1+n2,n);
             else
                 indices = randi(n1+n2,1,n);
             end
             if state==1 || state==2
                 X = find([obj.elements.existent] & [obj.elements.state]==state,max(indices));
             else
                 X = find([obj.elements.existent],max(indices));
             end
             X = X(indices);
        end
        
        function X = choose_element_different_from(obj,ind,state)
            X = obj.choose_element(2,state,'different');
            if X(1)~=ind
                X = X(1);
            else
                X = X(2);
            end
        end
        
        function [success,rate,ode_ind,existent] = add_transition(obj,ind_from,ind_to)
            % if no ind_from and ind_to indices are provided choose elements randomly 
            if nargin==1
                X = obj.choose_element(2,0,'different');
                if isempty(X)
                    success = 0; rate=0; ode_ind=-1; return;
                end
                ind_from = X(1);
                ind_to = X(2);
            end
            % check if a transition with the same source and sink already exists. In this case just modify the respective rate, otherwise add a new transition and increase the connectivity of the respective element by 1
            ind = find([obj.transitions.from]==ind_from & [obj.transitions.to]==ind_to,1);     
            if ~isempty(ind)
                obj.transitions(ind).rate = obj.transitions(ind).rate + rand();
                existent = true;
            else 
                obj.transitions(end+1) = struct('rate',rand(),'from',ind_from,'to',ind_to);
                obj.elements(ind_to).increment_connectivity();       % increment connectivity of respective element 
                existent = false;
            end
            % transfer data for OdeSystem
            rate = obj.transitions(end).rate;
            ode_ind = [obj.elements(ind_from).ode_index;obj.elements(ind_to).ode_index];
            success = 1;
        end
        
        function [ind,state,rates,ode_inds,existent] = add_derived_element(obj,ode_ind)
            rates = zeros(1,2);
            ode_inds = int16(zeros(2,2));
            existent = logical(zeros(1,2));
            ind = find([obj.elements.existent]==false,1);
            if isempty(ind)
                ind = length(obj.elements) + 1;
            end
            state = randi(2);
            obj.elements(ind) = Element(false,state,0,ode_ind);
            obj.increment_n_state(state);
            % choose two elements to connect the new one to. transitions are allowed from each state to each other state
            [~,rates(1),ode_inds(:,1),existent(1)] = obj.add_transition(obj.choose_element_different_from(ind,0),ind);       
            [~,rates(2),ode_inds(:,2),existent(2)] = obj.add_transition(ind,obj.choose_element_different_from(ind,0));
        end
        
        function delete_transition(obj,ind)         % delete the transition with index ind and decrement the connectivity of the element the transition goes to
              ind_elem = obj.transitions(ind).to;
              % decrease connectivity of respective element by 1
              obj.elements(ind_elem).decrement_connectivity();      
              % copy last element onto the deleted element and erase the last
              obj.transitions(ind) = obj.transitions(end);      
              obj.transitions(end) = [];
        end
        
        function delete_derived_element(obj,ind)       % delete a derived element and all connecting transitions to and from this element
            % delete element,i.e. set 'existent' to false and decrement the number of corresponding states in 'Subnetwork'
            obj.decrement_n_state(obj.elements(ind).state);
            obj.elements(ind).existent = false;  
            % find and delete all transitions that go from or to this respective element. delete_transition() updates the connectivity index of the elements accordingly
            trans_ind_del = find([obj.transitions.from]==ind | [obj.transitions.to]==ind);     
            trans_ind_del = fliplr(trans_ind_del);   % sort in descending order so that deletion commands in the next foor loop do not intersect
            for i=1:length(trans_ind_del)
                obj.delete_transition(trans_ind_del(i));    
            end
        end
        
        function ind = delete_derived_element_cascade(obj,ind)      % delete a derived element and all other derived elements which, as a result of the deletion, become disconnected from the network; return the index of the element from which the cascade is started (return 0 if no appropriae element can be found)
            if nargin==1
                ind = find_random([obj.elements.existent] & [obj.elements.base_element]==false);
                if isempty(ind)
                    ind = 0; return;
                end
            end
            % delete the given element and hence also all connecting transitions to and from that element which automatically updates the connectivity indices of all other elements accordingly
            obj.delete_derived_element(ind);
            % find all elements (that are existent), which are no base elements and have connectivity zero (as a result of the transition deletions) and hence cannot acquire a nonzero density and delete these elements as well; repeat this step so long until no further such elements can be found anymore  
            while true
                elems_ind_del = find([obj.elements.existent]==true & [obj.elements.connectivity]==0 & [obj.elements.base_element]==false);          
                if isempty(elems_ind_del)    % if elems_ind_del is empty break and finish the function, otherwise delete the found unconnected elements
                    break;
                else
                    for i=1:length(elems_ind_del)
                        obj.delete_derived_element(elems_ind_del(i));
                    end
                end
            end
        end
        
        function [success,rates] = mutate_transition_rate(obj)       % mutate the rate of a randomly chosen transition by multiplying it with a uniformly distributed random factor between 0 and 2
            rates = zeros(1,2);  % [old_rate,new_rate]
            if isempty(obj.transitions)
                success = 0; return;
            else
                ind = randi(length(obj.transitions),1,1);
                rates(1) = obj.transitions(ind).rate;     % old_rate
                obj.transitions(ind).rate = obj.transitions(ind).rate * 2*rand();
                rates(2) = obj.transitions(ind).rate;      % new_rate
                success = 1;
            end
        end
        
        %% Processing for ode solver
 function states = assign_ode_index_process_states(obj,index_begin,index_end)
     t=num2cell(index_begin:index_end);
     [obj.elements([obj.elements.existent]).ode_index] = t{:};   % alternatively you can do this with disperse
     states = [obj.elements([obj.elements.existent]).state];
 end
 function [Rates,Inds] = process_transitions(obj)       % returns the transitions as ode_indexes and the rates of the transitions of the subnetwork in order to make the object OdeTransitions and OdeTransitions_rates on Network level
     Rates = [obj.transitions.rate];
     Inds = [ [obj.elements([obj.transitions.from]).ode_index];[obj.elements([obj.transitions.to]).ode_index] ];
 end
 
    end
end

