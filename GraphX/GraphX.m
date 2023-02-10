classdef GraphX

    properties (SetAccess = private)
        nNodes;                 % Number of nodes in the graph
        nEdges;
        adjacencyMatrix;        % Adjacency matrix of the RDPG
        graphObject;            % Matlab graph object associated to the graph
        graphLaplacian;         % Laplacian matrix
        eigenvalues;            % Eigenvalues of Laplacian 
        eigenvectors;           % Eigenvectors of Laplacian
        degreeVector;
        normalizedLaplacianMatrix;
        normalizedLaplacianEigenvalues;
        normalizedLaplacianEigenvectors;
    end

    methods(Static)

        function cycleGraph = cycleGraphX(nNodes)

            A = zeros(nNodes);

            for i=1:nNodes
                for j=1:nNodes
                    if mod(i-j, nNodes) ==1
                        A(i,j) = 1;
                        A(j,i) = 1;
                    end

                end
            end

            cycleGraph = GraphX(A);
        end

    end

    methods

        %% Constructor based on adj. matrix
        function graphX = GraphX(adjacencyMatrix) 
            graphX.adjacencyMatrix = adjacencyMatrix;
            graphX.nNodes = size(adjacencyMatrix, 1);
            graphX.nEdges = sum(sum(adjacencyMatrix)) / 2;
            graphX.graphObject = graph(adjacencyMatrix); 
            graphX.graphLaplacian = laplacian(graphX.graphObject);
            [graphX.eigenvectors, graphX.eigenvalues] = graphX.getLaplacianEigenvectors(); % Auto-Populate the laplacian spectral info

            graphX.degreeVector = sum(graphX.adjacencyMatrix, 1);

            d = graphX.degreeVector .^(-1/2);
            d(isinf(d)|isnan(d)) = 0; 
            graphX.normalizedLaplacianMatrix = diag(d) * graphX.graphLaplacian * diag(d);
            [graphX.normalizedLaplacianEigenvectors, graphX.normalizedLaplacianEigenvalues] = graphX.getNormalizedLaplacianEigenvectors();

        end

        %% Grab the spectral info for laplacian matrix, sorted by magnitude
        function [sortedEvecs, sortedEvals] = getLaplacianEigenvectors(graphX)
            L = graphX.graphLaplacian;
            L = full(L); % Convert from sparse array
            [evecs, evals] = eig(L);
            [~, ind] = sort(diag(evals),'ComparisonMethod','abs');
            %ind = flip(ind);
            sortedEvals = diag(evals(ind,ind));
            sortedEvecs = evecs(:,ind);
        end

        function [sortedEvecs, sortedEvals] = getNormalizedLaplacianEigenvectors(graphX)
            L = graphX.normalizedLaplacianMatrix;
            L = full(L); % Convert from sparse array
            [evecs, evals] = eig(L);
            [~, ind] = sort(diag(evals),'ComparisonMethod','abs');
            %ind = flip(ind);
            sortedEvals = diag(evals(ind,ind));
            sortedEvecs = evecs(:,ind);
        end

        %% Grab the spectral info for the Adjacency matrix, sorted by magnitude
        function [sortedEvecs, sortedEvals] = getAdjacencyEigenvectors(graphX)
            A = graphX.adjacencyMatrix;
            A = full(A);
            [evecs, evals] = eig(A);
            [~, ind] = sort(diag(evals),'ComparisonMethod','abs');
            ind = flip(ind);
            sortedEvals = diag(evals(ind,ind));
            sortedEvecs = evecs(:,ind);
        end

    end

end