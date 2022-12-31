ic = [1;2;3;4];
N = 4;
epsilon = .4; 

s1 = discrete_switch_simulation(N,.1,epsilon,ic); 
s2 = discrete_switch_simulation(N,.9,epsilon,ic); 


d_p1 = consensus_and_disagreement(4,s1, ic);
d_p2 = consensus_and_disagreement(4,s2, ic); 

matObj = matfile('adjacency_matrices.mat'); 
A2 = matObj.A2;
L = graph_laplacian(A2); 
e = .4;
y = zeros(length(ic));
y(:,1) = ic; % initialize y
for k = 1:19 %20 time steps assumed
    y(:,k+1) = (eye(4) - e*L)  * y(:,k); 
end
d_A2 = consensus_and_disagreement(4,y,ic);

t = 1:20;
figure, 
nexttile
plot(t,d_p1,'k-'), hold on
title("Disagreement Plot for p = .1")       
xlabel('t'), ylabel('x')

nexttile
plot(t,d_p2,'k-'), hold on 
title("Disagreement Plot for p = .9")   
xlabel('t'), ylabel('x')

nexttile
plot(t,d_A2,'k-'), hold on
title("Disagreement Plot for A2")   
xlabel('t'), ylabel('x')

function d = consensus_and_disagreement(N,x,ic)
    D = zeros(20,1);
    sum = 0;
    for j = 1:N
       sum = sum + ic(j,1);
    end
    xbar = (1/N) * sum;
    
    for k = 1:20 %20 time steps assumed
        D(:,k) = norm(x(:,k) - ((xbar)*ones(N,1))); 
    end
    d = (D(1,:)); 
    d = d';
end

function s = discrete_switch_simulation(N,p,e,initial)
[~,L] = random_graph(N,p); 
y1 = zeros(length(initial));
y1(:,1) = initial; % initialize y
for k = 1:19 %20 time steps assumed
    y1(:,k+1) = (eye(4) - e*L)  * y1(:,k); 
end
s = y1; 
end


% random graph script
function [Ar,Lr] = random_graph(N, p)
% Make an NxN matrix of uniformly distributed real numbers between 0 and 1.
% If they are <p, make them 1's. If they are >p, make them 0's.
Ar = rand(N,N)<p;   

% Redefine that matrix as its upper triangle plus its upper triangle
% transposed. This achieves getting a symmetric matrix.
Ar = triu(Ar,1) + triu(Ar,1)';

% Subtract the diagonal to make it all zeros. The line above should already
% do this, but this is just to be sure. 
Ar = Ar-diag(diag(Ar));
    
% Compute the Laplacian
Lr = diag(sum(Ar,2))-Ar;
end

function [Laplace,eigenvalue,eigenvector] = graph_laplacian(A)
  
    L = diag(sum(A,2))-A;  

    [V,D] = eig(L);
    Laplace = L;
    eigenvalue = D; 
    eigenvector = V; 
end