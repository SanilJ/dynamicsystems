%epsilon finder for convergence in 
% 4.1.2
matObj = matfile('adjacency_matrices.mat'); 

A1 = matObj.A1; 
A2 = matObj.A2; 
A3 = matObj.A3; 

e1 = epsilon_finder(A1);
e2 = epsilon_finder(A2);
e3 = epsilon_finder(A3);
disp(e1); 
disp(e2); 
disp(e3);

function e = epsilon_finder(A) 
episilon = 0; 

[L,~,~] = graph_laplacian(A); 
n = length(A); 
In = eye(n); 

xk = ones(n); 
pk = (In - episilon*L) * xk; 
episilon = episilon + .1; 
% first iteration happens outside
k = 1; 
while episilon < 1
    pk_change = (In - episilon*L) * pk; 

    if(pk_change == pk)  
        e = episilon; 
        break; 
    end
    
    pk = pk_change;
    episilon = episilon + .1; 
end

end

function [Laplace,eigenvalue,eigenvector] = graph_laplacian(A)
  
    L = diag(sum(A,2))-A;  

    [V,D] = eig(L);
    Laplace = L;
    eigenvalue = D; 
    eigenvector = V; 
end
