function ind = find_random(vec)             % returns the index of a random non-zero element of the vector vec; ususally vec is a logical vector; to increase efficiency this whole function could also be written as a mex function !!!!! TODO  
         % maybe provide another similar function 'find_n_distinct()' that finds n random distinct indices; with these two functions all calls the 'choose_element' and 'choose_subnetwork' could be replaced making the evoalg more compact and clearer    
             s = sum(logical(vec));
             if s==0
                 ind = []; return;
             end
       %      ind = find_ith(vec,randi(s,1,1))   % for better efficiency write this mex function
             ind = find(vec,randi(s,1,1));
             ind = int16(ind(end));
end

