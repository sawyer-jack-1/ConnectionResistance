
set(0,'DefaultFigureWindowStyle','docked')
clear all;
close all;

%%

DIM_CONNECTION = 2;
TEST_EDGE = [1, 2];

%%

%randSig1 = ConnectionGraphX.getRandomSOMatrix(DIM_CONNECTION);

%% Use edge synchronization map to find the (connection) resistance

t = tiledlayout(4, 4);
nNodes = [3, 4, 15, 32];
nThetaSamples = 200;

w = waitbar(0);

for n = nNodes

    ax = nexttile;
    
    CYCLE_LENGTH = n;
    DIM_CONNECTION = 2;
    
    G = GraphX.cycleGraphX(CYCLE_LENGTH);
    G.graphObject.plot()
    title(ax, strcat("Cycle Length: ", num2str(CYCLE_LENGTH)))
end

for n = nNodes

    ax = nexttile;
    
    CYCLE_LENGTH = n;
    DIM_CONNECTION = 2;
    TEST_EDGE = [1, 2];
    
    G = GraphX.cycleGraphX(CYCLE_LENGTH);
    G = ConnectionGraphX(G,DIM_CONNECTION);
    
    theta = linspace(0,360,nThetaSamples);
    resistances = arrayfun(@(x) resistanceValuePsuedo(G, TEST_EDGE, x), theta);
    
    hold on

    plot(theta,resistances);
    plot(theta,repmat(getBasicResistance(G,TEST_EDGE), nThetaSamples));
    title(ax, 'Mark and Fans Definition');
    
    hold off

    waitbar(n/32, w);

end

for n = nNodes

    ax = nexttile;
    
    CYCLE_LENGTH = n;
    DIM_CONNECTION = 2;
    TEST_EDGE = [1, 2];
    
    G = GraphX.cycleGraphX(CYCLE_LENGTH);
    G = ConnectionGraphX(G,DIM_CONNECTION);
    
    theta = linspace(0,360,nThetaSamples);
    resistances = arrayfun(@(x) resistanceValueSchur(G, TEST_EDGE, x), theta);
    
    hold on

    plot(theta,resistances);
    plot(theta,repmat(getBasicResistance(G,TEST_EDGE), nThetaSamples));
    title(ax, 'Schur complement');
    
    hold off

    waitbar(n/32, w);

end


for n = nNodes

    ax = nexttile;
    
    CYCLE_LENGTH = n;
    DIM_CONNECTION = 2;
    TEST_EDGE = [1, 2];
    
    G = GraphX.cycleGraphX(CYCLE_LENGTH);
    G = ConnectionGraphX(G,DIM_CONNECTION);
    
    theta = linspace(0,360,nThetaSamples);
    resistances = arrayfun(@(x) resistanceValue(G, TEST_EDGE, x), theta);
    
    hold on

    plot(theta,resistances);
    plot(theta,repmat(getBasicResistance(G,TEST_EDGE), nThetaSamples));
    title(ax, 'Eigenvector Method');
    
    hold off

    waitbar(n/32, w);

end

close(w);

folder = 'G:\My Drive\matlab_projects\connection-graphs';
filename = fullfile(folder, 'resistance_methods.png');
exportgraphics(t, filename, 'Resolution', 600);

%% Resistance Function

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

function r = resistanceValue(connectionGraphX, testEdge, theta)

        d = 1;
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

function S = SchurComp(A,ind)
% Author   : Zhengchao Wan
%            The Ohio State University
%            wan.252@osu.edu
% Last Rev : Fri Mar 5 12:54:08 EST 2021

% Copyright notice: You are free to modify, extend and distribute 
%    this code granted that the author of the original code is 
%    mentioned as the original author of the code.
    [l,w] = size(A);
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
