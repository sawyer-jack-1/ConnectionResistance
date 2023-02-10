
set(0,'DefaultFigureWindowStyle','docked')
clear all;
close all;
clc;

%%

cycleLength = 8;
dimConnection = 4;

nTestTrials = 1e2;

%%

G = GraphX.cycleGraphX(cycleLength);
G = ConnectionGraphX(G,dimConnection);

%G = G.setEdgeConnection(1,2,flip(eye(dimConnection)));

randSig = ConnectionGraphX.getRandomSOMatrix(dimConnection);
G = G.setEdgeConnection(1,2, randSig);

%det(randSig)
%eig(randSig)

%%

trialResistances = zeros(cycleLength,cycleLength, nTestTrials);

f = waitbar(0, "Please wait...");

for k=1:nTestTrials
    waitbar(k / nTestTrials, f, "Sampling...");
    R = zeros(cycleLength, cycleLength);
    x = randn(dimConnection,1);
    x = x / norm(x);

    for i = 1:cycleLength
        for j = 1:cycleLength
               
            if G.adjacencyMatrix(i,j) == 1
    
                d = dimConnection;
                m = zeros(cycleLength * d, 1);
                m((d*(i-1) + 1):(d * i),:) = x;
                m((d*(j-1) + 1):(d * j),:) = - G.connectionAdjacency((d*(j-1) + 1):(d * j),(d*(i-1) + 1):(d * i)).' * x;

                R(i,j) = (m.' * lsqminnorm(G.connectionLaplacian, m));
                
            end
    
        end
    end

    trialResistances(:,:,k) = R;
end    

waitbar(.5,f,"Calculating stats...");

var(trialResistances,0,3)

R

figure;
histogram(trialResistances(1,2,:), 20, 'Normalization','pdf');
title("Resistances between modified edge \{1,2\}");

figure;
histogram(trialResistances(3,4,:), 20, 'Normalization','pdf');
title("Resistances between other edges");

close(f)



