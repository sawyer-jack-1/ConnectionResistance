
set(0,'DefaultFigureWindowStyle','docked')
clear all;
close all;

%%

DIM_CONNECTION = 2;
TEST_EDGE = [1, 2];

%%

%randSig1 = ConnectionGraphX.getRandomSOMatrix(DIM_CONNECTION);

%% Use edge synchronization map to find the (connection) resistance

NODES = 10;
DIM_SAMPLES = 2;
SPARSITY = 0.25;
METHOD = "Square";
POSITIONS = [];

G = RandomDotProductGraph(NODES, DIM_SAMPLES, SPARSITY, METHOD, POSITIONS);
G.draw()

for n = nNodes
    ax = nexttile;
    
    CYCLE_LENGTH = n;
    DIM_CONNECTION = 2;
    TEST_EDGE = [1, 2];
    
    G = GraphX.cycleGraphX(CYCLE_LENGTH);
    G = ConnectionGraphX(G,DIM_CONNECTION);
    
    theta = linspace(0,360,1e3);
    resistances = arrayfun(@(x) resistanceValue(G, TEST_EDGE, x), theta);
    
    plot(theta,resistances);
    title(ax, strcat("Cycle length: ",num2str(n))');
    
    waitbar(n/512, w);

end

close(w);
folder = 'G:\My Drive\matlab_projects\connection-graphs';
filename = fullfile(folder, 'cycle_resistances.png');
exportgraphics(t, filename, 'Resolution', 300);

%% Resistance Function

function r = resistanceValue(connectionGraphX, testEdge, theta)

        d = 2;
        dimConnection_ = connectionGraphX.dimConnection;
        n = connectionGraphX.nNodes;

        i = testEdge(1);
        j = testEdge(2);

        randSig1 = [cosd(theta), -sind(theta) ; sind(theta), cosd(theta)];
        connectionGraphX = connectionGraphX.setEdgeConnection(i, j, randSig1);
        L_inv = pinv(connectionGraphX.connectionLaplacian);

        edge_synchronization_map = connectionGraphX.connectionEigenvectors(:,1:d);
        
        displacement_vector = zeros( n * dimConnection_, d);
        displacement_vector( dimConnection_*(i-1) + 1: i * dimConnection_, :) = edge_synchronization_map(dimConnection_*(i-1) + 1: i * dimConnection_,:);
        displacement_vector( dimConnection_*(j-1) + 1: j * dimConnection_, :) = -edge_synchronization_map(dimConnection_*(j-1) + 1: j * dimConnection_,:);
        
        displacement_vector = sqrt(2) * displacement_vector ./ (norm(displacement_vector, 'fro'));

        r = norm(displacement_vector' * L_inv * displacement_vector, "fro");


end

   
% Set g to identity at various vertices to get different choices for g. Try and calculate effective resistance.
% Fix the randSig and then run code a bunch of times to see how it behaves.
