function [pl,ql,pr,qr] = bc_fun_pdepe(~,~,~,~,~,N_elems)
  pl = zeros(N_elems,1);
  ql = ones(N_elems,1);
  pr = zeros(N_elems,1);
  qr = ones(N_elems,1);
end

