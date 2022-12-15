function X = colcomp(A,b,n_row)         % with n_row specified: compares the columns of (2,n) matrix A with the (2,1) vector b and returns an (1,n) array of logicals indicating whether a column in A equals b  
    if isempty(A)                       % without n_row specified: compares only the row with index n_row of A with the integer b and returns an (1,n) array of logicals
        X=[]; return;
    end
    if nargin==2
        X = (A(1,:)==b(1)) & (A(2,:)==b(2));
    else
        X = (A(n_row,:)==b);
    end
end

