ic = [1;2;3;4];
save('initial_conditions.mat'); 
load('initial_conditions.mat'); 

matObj = matfile('adjacency_matrices.mat'); 

A1 = matObj.A1; 
A2 = matObj.A2; 
A3 = matObj.A3;

e1 = .4; % in
e2 = 1; % out of 
y1_small = zeros(length(ic)); y1_big = zeros(length(ic));
y2_small = zeros(length(ic)); y2_big = zeros(length(ic));
y3_small = zeros(length(ic)); y3_big = zeros(length(ic));

y1_small(:,1) = ic; y1_big(:,1) = ic;   % initialize y
y2_small(:,1) = ic; y2_big(:,1) = ic;   % initialize y
y3_small(:,1) = ic; y3_big(:,1) = ic;   % initialize y

L1 = graph_laplacian(A1);  
L2 = graph_laplacian(A2);   
L3 = graph_laplacian(A3);
for k = 1:19
    y1_small(:,k+1) = (eye(4) - e1*L1)  * y1_small(:,k); 
    y1_big(:,k+1) = (eye(4) - e2*L1)  * y1_big(:,k); 
    
    y2_small(:,k+1) = (eye(4) - e1*L2)  * y2_small(:,k); 
    y2_big(:,k+1) = (eye(4) - e2*L2)  * y2_big(:,k); 

    y3_small(:,k+1) = (eye(4) - e1*L3)  * y3_small(:,k); 
    y3_big(:,k+1) = (eye(4) - e2*L3)  * y3_big(:,k); 

end

t = 1:20;
figure, 

nexttile
plot(t,y1_small,'k-')
title("A1 with Epsilon = .4")       
xlabel('t'), ylabel('x')

nexttile
plot(t,y1_big,'k-')
title("A1 with Epsilon = 1")
xlabel('t'), ylabel('x')

nexttile
plot(t,y2_small,'k-')
title("A2 with Epsilon = .4")       
xlabel('t'), ylabel('x')

nexttile
plot(t,y2_big,'k-')
title("A2 with Epsilon = 1")
xlabel('t'), ylabel('x')

nexttile
plot(t,y3_small,'k-')
title("A3 with Epsilon = .4")       
xlabel('t'), ylabel('x')

nexttile
plot(t,y3_big,'k-')
title("A3 with Epsilon = 1")
xlabel('t'), ylabel('x')

function [Laplace,eigenvalue,eigenvector] = graph_laplacian(A)
  
    L = diag(sum(A,2))-A;  

    [V,D] = eig(L);
    Laplace = L;
    eigenvalue = D; 
    eigenvector = V; 
end


