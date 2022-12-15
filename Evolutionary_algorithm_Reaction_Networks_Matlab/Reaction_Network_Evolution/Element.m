classdef Element < handle
   properties
      existent(1,1) logical       
      base_element(1,1) logical   % base_element==true -> base element, i.e. either a protein with initial concentration or a composite structure (e.g. dimer) with an ingoing reaction
%      composite(1,1) logical
      state(1,1) int8             % state==1 -> cytosolic, state==2 -> membrane-bound
      connectivity(1,1) int16      % number of ingoing transitions and ingoing decays. if connectivity==0 and the respective element is not a base_element, then the element is completely decoupled from the network and its concentration will always be zero. Therefore it can be deleted. 
      init_density(1,1) single    % initial concentration; only relevant if base_element==true
      ode_index(1,1) int16    
   end
   methods
       function obj = Element(base_element_,state_,init_density_,ode_index_)
           obj.existent = true;
           obj.base_element = base_element_;
           obj.state = state_;
           obj.connectivity = 0;
           obj.init_density = init_density_;
           obj.ode_index = ode_index_;
       end
       
       function increment_connectivity(obj)
            obj.connectivity = obj.connectivity + 1;
        end
        
        function decrement_connectivity(obj)
            obj.connectivity = obj.connectivity - 1;
        end
       
   end
end












