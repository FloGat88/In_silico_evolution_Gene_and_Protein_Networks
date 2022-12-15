function jac = Jacobian(~,x,OdeTransitions_rates,OdeReactions_rates,OdeDecays_rates,OdeTransitions,OdeReactions,OdeDecays,D_eff,N_el,N_L)
     Jac_Diff = spdiags(D_eff',0,N_el,N_el);
     row = [OdeTransitions(1,:),OdeTransitions(2,:),OdeDecays(1,:),OdeDecays(2,:),OdeDecays(3,:)];   
     col = [OdeTransitions(1,:),OdeTransitions(1,:),OdeDecays(1,:),OdeDecays(1,:),OdeDecays(1,:)]; 
     el = [-OdeTransitions_rates,OdeTransitions_rates,-OdeDecays_rates,OdeDecays_rates,OdeDecays_rates];
     Jac_Ode1 = sparse(row,col,el,N_el,N_el);       % for some strange reason sparse accepts the indices only as doubles 
     Jac_DiagBlocks = cell(N_L,1);
     for i=1:N_L
         x_ = x((i-1)*N_el+1 : i*N_el)';
         x1 = x_(OdeReactions(1,:)) .* OdeReactions_rates;
         x2 = x_(OdeReactions(2,:)) .* OdeReactions_rates;
         row = [OdeReactions(1,:),OdeReactions(1,:),OdeReactions(2,:),OdeReactions(2,:),OdeReactions(3,:),OdeReactions(3,:)];
         col = [OdeReactions(1,:),OdeReactions(2,:),OdeReactions(1,:),OdeReactions(2,:),OdeReactions(1,:),OdeReactions(2,:)];
         el = [-x2,-x1,-x2,-x1,x2,x1];
         if i==1 || i==N_L
            Jac_DiagBlocks{i} = sparse(row,col,el,N_el,N_el) + Jac_Ode1 - Jac_Diff;
         else
             Jac_DiagBlocks{i} = sparse(row,col,el,N_el,N_el) + Jac_Ode1 - 2*Jac_Diff;
         end
     end
 %    jac = blkdiag(Jac_DiagBlocks{:}) + spdiags(repmat(D_eff',N_L-1,2),[-N_el,N_el],N_el*N_L,N_el*N_L);
     jac =  spdiags(repmat(D_eff',N_L,2), [-N_el,N_el], blkdiag(Jac_DiagBlocks{:}));
end
 
 %{
 function jac = Jacobian2(obj,~,x,D_eff,N_L)
     N_el = N_elems;
     Jac_Diff = spdiags(D_eff',0,N_el,N_el);
     row = [OdeTransitions(1,:),OdeTransitions(2,:),OdeDecays(1,:),OdeDecays(2,:),OdeDecays(3,:)];
     col = [OdeTransitions(1,:),OdeTransitions(1,:),OdeDecays(1,:),OdeDecays(1,:),OdeDecays(1,:)]; 
     el = [-OdeTransitions_rates,OdeTransitions_rates,-OdeDecays_rates,OdeDecays_rates,OdeDecays_rates];
     Jac_Ode1 = sparse(row,col,el,N_el,N_el);
     nonzeros = N_L * (nnz(jac_Ode1) + 6*length(OdeReactions_rates) + N_el*3);
     jac = sparse([],[],[],N_el*N_L,N_el*N_L,nonzeros);
     for i=1:N_L
         range = (i-1)*N_el+1 : i*N_el;
         x_ = x(range);
         x1 = x_(OdeReactions(1,:)) .* OdeReactions_rates;
         x2 = x_(OdeReactions(2,:)) .* OdeReactions_rates;
         row = [OdeReactions(1,:),OdeReactions(1,:),OdeReactions(2,:),OdeReactions(2,:),OdeReactions(3,:),OdeReactions(3,:)];
         col = [OdeReactions(1,:),OdeReactions(2,:),OdeReactions(1,:),OdeReactions(2,:),OdeReactions(1,:),OdeReactions(2,:)];
         el = [-x2,-x1,-x2,-x1,x2,x1];
         if i==1 || i==N_L
            jac(range,range) = sparse(row,col,el,N_el,N_el) + Jac_Ode1 - Jac_Diff;
         else
            jac(range,range) = sparse(row,col,el,N_el,N_el) + Jac_Ode1 - 2*Jac_Diff;
         end
     end
     jac = jac + spdiags(repmat(D_eff',N_L,2),[-N_el,N_el]);  
 end
%}