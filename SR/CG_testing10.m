
set(0,'DefaultFigureWindowStyle','docked')
clear all;
close all;
clc;

%%

NODES = 4;

%DIM_SAMPLES = 2;
%SPARSITY = 1;
%METHOD = "Square";
%POSITIONS = [];
%G = RandomDotProductGraph(NODES, DIM_SAMPLES, SPARSITY, METHOD, POSITIONS);

G = GraphX.cycleGraphX(NODES);
A0 = G.adjacencyMatrix;
A0(1,3) = 1;
A0(3,1) = 1;
G = GraphX(A0);

tiledlayout(1,2);

nexttile;
plot(G.graphObject);
title('Plot of random graph G');

resistanceMatrix = zeros(G.nNodes);

for i=1:G.nNodes
    
    for j=1:G.nNodes
        if i ~= j
            
            resistanceMatrix(i,j) = getBasicResistance(G, [i,j]);
 
        end     
    end      
end

nexttile;
imagesc(resistanceMatrix);
colorbar;
title('Color plot of the standard resistance matrix');

firstEdge = find(full(G.graphLaplacian),2,'first').';

fprintf(strcat('test edge: \n', num2str(firstEdge)));

fprintf('\nStandard conductance formula for i,j:');
resistanceMatrix(firstEdge(1),firstEdge(2))^(-1)

fprintf('Schur complement matrix L_{e^c,e^c}:');
A = SchurComp(full(G.graphLaplacian),firstEdge);
A

fprintf('===================================\n');

DIM_CONNECTION = 2;
CONNECTION_MATRIX = [0, 1; 1,0];

%CONNECTION_MATRIX = RandOrthMat(DIM_CONNECTION);

G = ConnectionGraphX(G,DIM_CONNECTION);
G = G.setEdgeConnection(firstEdge(1), firstEdge(2), CONNECTION_MATRIX);

firstEdge = [2, 3];
fprintf('\nTest edge: {2,3} ');

fprintf('\nStandard conductance formula for i,j:');
resistanceMatrix(firstEdge(1),firstEdge(2))^(-1)

fprintf('\nConn. Schur complement matrix L_{e^c,e^c}:');
A = SchurComp(full(G.connectionLaplacian),[firstEdge(1)*DIM_CONNECTION-1, firstEdge(1)*DIM_CONNECTION, firstEdge(2)*DIM_CONNECTION-1, firstEdge(2)*DIM_CONNECTION]);
A

fprintf('\nConnection Resistance between nodes as (1/2d*trace(C^sigma_ij))^(-1)');

connectionResistanceMatrix = zeros(G.nNodes);

for i=1:G.nNodes
    for j=1:G.nNodes
        if i ~= j
            connectionResistanceMatrix(i,j) = traceConnectionResistance(G, [i,j]);
        end     
    end      
end

connectionResistanceMatrix

fprintf('\nStandard Resistance between nodes as e_ij^T L^+ e_ij');

resistanceMatrix

fprintf('===================================\n');

fprintf('\nConnection Resistance DISTANCE sqrt(R^sigma_ij)');

connectionResistanceMatrix.^(1/2)

[V,D] = eig(G.connectionLaplacian);
v1 = V(firstEdge(1),:);
v2 = V(firstEdge(2),:);

d = (v1 - v2).^2;
d(1) = 0;
%d(2) = 0;
d

lambda = diag(D);
lambda = lambda.^(-1);
lambda(1) = 0;
%lambda(2) = 0;
lambda

lambda.'*(d')


%% Resistance Functions

function r = getBasicResistance(connectionGraphX, testEdge)
   
    L_pinv = pinv(full(connectionGraphX.graphLaplacian));

    id = eye(connectionGraphX.nNodes);
    e = id(:, testEdge(1)) - id(:, testEdge(2));

    r = e' * L_pinv * e;

end

function r = resistanceValuePsuedo(connectionGraphX, testEdge, theta)

        d = 2;

        i = testEdge(1);
        j = testEdge(2);

        rotationTheta = [cosd(theta), -sind(theta) ; sind(theta), cosd(theta)];
        connectionGraphX = connectionGraphX.setEdgeConnection(i, j, rotationTheta);

        L_inv = pinv(connectionGraphX.connectionLaplacian);

        edge_matrix = [rotationTheta; -eye(d)];

        pinvBlock = [   L_inv(d * (i - 1) + 1 : i * d, d * (i - 1) + 1 : i * d), ...
                        L_inv(d * (i - 1) + 1 : i * d, d * (j - 1) + 1 : j * d); ...
                        L_inv(d * (j - 1) + 1 : j * d, d * (i - 1) + 1 : i * d), ...
                        L_inv(d * (j - 1) + 1 : j * d, d * (j - 1) + 1 : j * d) ];
    
        %r = norm( edge_matrix.' * pinvBlock * edge_matrix, 2);

        Psi = connectionGraphX.connectionIncidence * L_inv * connectionGraphX.connectionIncidence';
        r = norm( Psi(1:d, 1:d), 2);
end

function r = resistanceValueSchur(connectionGraphX, testEdge, theta)

        d = 2;

        i = testEdge(1);
        j = testEdge(2);

        rotationTheta = [cosd(theta), -sind(theta) ; sind(theta), cosd(theta)];
        connectionGraphX = connectionGraphX.setEdgeConnection(i, j, rotationTheta);
        L = connectionGraphX.connectionLaplacian;

        a = d * (i - 1) + 1 : i * d;
        b = d * (j - 1) + 1 : j * d;
        Lab = SchurComp(L, [a b]);
        r = 1 ./ norm(Lab(a,b), 2);

end

function r = traceConnectionResistance(connectionGraphX, testEdge)

        DIM_CONNECTION = connectionGraphX.dimConnection;
        
        schurComplementConductance = SchurComp(full(connectionGraphX.connectionLaplacian),[testEdge(1)*DIM_CONNECTION-1, testEdge(1)*DIM_CONNECTION, testEdge(2)*DIM_CONNECTION-1, testEdge(2)*DIM_CONNECTION]);
        
        r = (1 / (2 * DIM_CONNECTION) * trace(schurComplementConductance) )^(-1);
        
end

function S = SchurComp(A,ind)
% Author   : Zhengchao Wan
%            The Ohio State University
%            wan.252@osu.edu
% Last Rev : Fri Mar 5 12:54:08 EST 2021

% Copyright notice: You are free to modify, extend and distribute 
%    this code granted that the author of the original code is 
%    mentioned as the original author of the code.
    [~,w] = size(A);
    nind = setdiff(1:w,ind);
    Saa = A(ind,ind);
    Sab = A(ind,nind);
    Sba = A(nind,ind);
    Sbb = A(nind,nind);
    % The following formula was supposed to be Saa-Sab * pinv(Sbb) * Sba
    % However, we use lsqminnorm(A,b) to efficiently implement pinv(A)*b
    S = Saa-Sab * lsqminnorm(Sbb,Sba);

end

% Set g to identity at various vertices to get different choices for g. Try and calculate effective resistance.
% Fix the randSig and then run code a bunch of times to see how it behaves.
