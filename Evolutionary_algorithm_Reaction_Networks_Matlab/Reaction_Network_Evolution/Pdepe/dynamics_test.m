function [c,f,s] = dynamics_test(x,t,u,dudx)
    c = ones(8,1);
    f = dudx;
    s = zeros(8,1);
end

