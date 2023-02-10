
set(0,'DefaultFigureWindowStyle','docked')
clear all;
close all;

%%

CYCLE_LENGTH = 3;
DIM_CONNECTION = 3;
TEST_EDGE = [1, 2];
TEST_EDGE2 = [2, 3];
TEST_EDGE3 = [3, 4];


%%


G = GraphX.cycleGraphX(CYCLE_LENGTH);
G = ConnectionGraphX(G,DIM_CONNECTION);

%randSig1 = ConnectionGraphX.getRandomSOMatrix(DIM_CONNECTION);

theta = 90;
randSig1 = [1, 0, 0; 0, cosd(theta), -sind(theta) ; 0, sind(theta), cosd(theta)];
%randSig2 = [cosd(theta), -sind(theta) ; sind(theta), cosd(theta)];
%randSig3 = [cosd(theta), -sind(theta) ; sind(theta), cosd(theta)];

%randSig = flip(eye(DIM_CONNECTION));
G = G.setEdgeConnection(TEST_EDGE(1),TEST_EDGE(2), randSig1);
%G = G.setEdgeConnection(TEST_EDGE2(1),TEST_EDGE2(2), randSig2);
%G = G.setEdgeConnection(TEST_EDGE3(1),TEST_EDGE3(2), randSig3);


%% Find the edge/angle synchronization map using relaxation.

edge_synchronization_map = G.connectionEigenvectors(:,1:DIM_CONNECTION);

%% Use edge synchronization map to find the (connection) resistance

resistance_matrix = zeros(CYCLE_LENGTH);
standard_resistance = zeros(CYCLE_LENGTH);
L_inv = pinv(G.connectionLaplacian);

for i=1:CYCLE_LENGTH
    for j=1:CYCLE_LENGTH
        if i ~= j
            
            d = DIM_CONNECTION;
            n = CYCLE_LENGTH;
            
            displacement_vector = zeros( CYCLE_LENGTH * DIM_CONNECTION, DIM_CONNECTION);
            displacement_vector( d*(i-1) + 1: i * d, :) = edge_synchronization_map(d*(i-1) + 1: i * d,:);
            displacement_vector( d*(j-1) + 1: j * d, :) = -edge_synchronization_map(d*(j-1) + 1: j * d,:);
            
            resistance_matrix(i,j) = norm(displacement_vector.' * L_inv * displacement_vector, "fro");
            
            id = eye(CYCLE_LENGTH);
            inv = pinv(full(G.graphLaplacian));
            standard_resistance(i,j) = (id(:,i) - id(:,j)).' * inv * (id(:,i) - id(:,j));
            
        end
    end
end

degreeBlock = zeros(CYCLE_LENGTH * DIM_CONNECTION);
for i=1:CYCLE_LENGTH
    for j=1:DIM_CONNECTION
        degreeBlock(d*(i-1) + j, d*(i-1) + j) = G.degreeVector(i);
    end
    
end

resistance_matrix
standard_resistance

norm(G.connectionIncidence * edge_synchronization_map, 'fro')


% Set g to identity at various vertices to get different choices for g. Try and calculate effective resistance.
% Fix the randSig and then run code a bunch of times to see how it behaves.
