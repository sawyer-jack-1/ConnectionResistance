
set(0,'DefaultFigureWindowStyle','docked')
clear all;
close all;
clc;

%%

CYCLE_LENGTH = 3;
DIM_CONNECTION = 2;
TEST_EDGE = [1, 2];

%%

G = GraphX.cycleGraphX(CYCLE_LENGTH);
G = ConnectionGraphX(G,DIM_CONNECTION);

randSig2 = ConnectionGraphX.getRandomSOMatrix(DIM_CONNECTION);
randSig = flip(eye(DIM_CONNECTION));
G = G.setEdgeConnection(1,2, randSig);

%% Find the edge/angle synchronization map using manopt.

    % Create the problem structure.
    manifold = stiefelstackedfactory(CYCLE_LENGTH, DIM_CONNECTION, DIM_CONNECTION);
    problem.M = manifold;

    % Define the problem cost function and its gradient.
    S = G.connectionIncidence;
    problem.cost  = @(x) trace( (S*x).' * (S*x) );
    problem.egrad = @(x) 2 * S * x;
    %problem.ehess = @(x, xdot) 2 * (S.' * S);

    % Solve.
    options.verbosity = 2;
    options.maxiter = 1e3;
    [edge_synchronization_map, xcost, info] = trustregions(problem, [], options); 

%% Use edge synchronization map to find the (connection) resistance

resistance_matrix = zeros(CYCLE_LENGTH);
resistance_matrix2 = zeros(CYCLE_LENGTH);
L_inv = pinv(G.connectionLaplacian);

for i=1:CYCLE_LENGTH
    for j=1:CYCLE_LENGTH
        if G.adjacencyMatrix(i,j)==1
            
            d = DIM_CONNECTION;
            n = CYCLE_LENGTH;
            
            displacement_vector = zeros( CYCLE_LENGTH * DIM_CONNECTION, DIM_CONNECTION);
            displacement_vector( d*(i-1) + 1: i * d, :) = edge_synchronization_map(d*(i-1) + 1: i * d,:);
            displacement_vector( d*(j-1) + 1: j * d, :) = -edge_synchronization_map(d*(j-1) + 1: j * d,:);
            
            resistance_matrix(i,j) = norm(displacement_vector.' * L_inv * displacement_vector, 'fro');
            
            displacement_vector( d*(i-1) + 1: i * d, :) = randSig2 * edge_synchronization_map(d*(i-1) + 1: i * d,:);
            displacement_vector( d*(j-1) + 1: j * d, :) = -randSig2 * edge_synchronization_map(d*(j-1) + 1: j * d,:);
            
            resistance_matrix2(i,j) = norm(displacement_vector.' * L_inv * displacement_vector, 'fro');
            
        end
    end
end

resistance_matrix
resistance_matrix2
    
% Set g to identity at various vertices to get different choices for g. Try and calculate effective resistance.
% Fix the randSig and then run code a bunch of times to see how it behaves.
