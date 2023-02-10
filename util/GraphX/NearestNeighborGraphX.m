classdef NearestNeighborGraphX < GraphX

    properties (SetAccess = private)
        latentPositions;       
        nNeighbors;  
        dimLatentPositions;
    end

    
    methods

        %% Constructor method for the RDPG with specified dimension, size, sparsity factor,
        % and method.
        function randomGraph = NearestNeighborGraphX(inputPositions, nNeighbors_)
            
            n = size(inputPositions, 1);
            A = zeros(n,n);

            for i=1:n
                
                posI = inputPositions(i, :);
                distances = zeros(1,n);

                for j=1:n
                    posJ = inputPositions(j, :);
                    distances(j) = norm(posJ - posI);
                end

                distances(i) = 1e10;

                for k=1:nNeighbors_
                    [~, node] = min(distances);
                    A(i,node) = 1;
                    A(node, i) = 1;
                    distances(node) = 1e10;
                end

            end

            randomGraph = randomGraph@GraphX(A);
            randomGraph.dimLatentPositions = size(inputPositions,2);
            randomGraph.latentPositions = inputPositions;
            randomGraph.nNeighbors = nNeighbors_;

        end

        %% Co-opting the matlab illustration feature to visualize a RDPG
        function draw(randomGraph)
            if randomGraph.dimLatentPositions == 1
                xData = randomGraph.latentPositions(:,1);
                plot(randomGraph.graphObject,'XData',xData);
            elseif randomGraph.dimLatentPositions == 2
                xData = randomGraph.latentPositions(:,1);
                yData = randomGraph.latentPositions(:,2);
                plot(randomGraph.graphObject,'XData',xData,'YData',yData);
            elseif randomGraph.dimLatentPositions >= 3
                xData = randomGraph.latentPositions(:,1);
                yData = randomGraph.latentPositions(:,2);
                zData = randomGraph.latentPositions(:,3);
                plot(randomGraph.graphObject,'XData',xData,'YData',yData,'ZData',zData);
            end
        end
    end

end

