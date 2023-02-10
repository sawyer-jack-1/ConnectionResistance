classdef ConnectionGraphX < GraphX

    properties (SetAccess = private)
        connection;
        dimConnection;
        connectionLaplacian;
        connectionAdjacency;
        connectionIncidence;
        connectionEigenvalues;
        connectionEigenvectors;
    end

    methods(Static)

        function bool = testOrthogonal(matrix)
            
            dim1 = size(matrix,1);
            dim2 = size(matrix,2);

            if dim1 ~= dim2 
                bool = 0;
                return;
            end

            if isequal(matrix.' * matrix, eye(dim1)) 
                bool = 1;
                return;
            end

        end

        function m = getRandomSOMatrix(d)

            success = 0;
            while success < 1
                randSig = RandOrthMat(d);
                if abs((det(randSig) - 1))<1e-3 
                    success = 1; 
                end
            end
            
            m = randSig;

        end

    end

    methods

        function connectionGraph = ConnectionGraphX(connectionGraphX_, dimConnection)
            
            connectionGraph = connectionGraph@GraphX(connectionGraphX_.adjacencyMatrix);
            connectionGraph.dimConnection = dimConnection;
            connectionGraph.connection = zeros(connectionGraphX_.nNodes, connectionGraphX_.nNodes, dimConnection, dimConnection);
                     
            for i = 1:connectionGraphX_.nNodes
                for j = 1:connectionGraphX_.nNodes
                    if connectionGraphX_.adjacencyMatrix(i,j) == 1
                        connectionGraph.connection(i,j,:,:) = eye(connectionGraph.dimConnection);
                    end
                end
            end

            connectionGraph = connectionGraph.setConnectionMatrices(connectionGraph.connection);

        end

        function connectionGraph = setConnectionMatrices(connectionGraph_, connection_)
            
            d = connectionGraph_.dimConnection;
            connectionGraph_.connectionAdjacency = zeros(connectionGraph_.nNodes * connectionGraph_.dimConnection, connectionGraph_.nNodes * connectionGraph_.dimConnection);
            connectionGraph_.connectionIncidence = zeros(connectionGraph_.nEdges * connectionGraph_.dimConnection, connectionGraph_.nNodes * connectionGraph_.dimConnection);
            connectionGraph_.connectionLaplacian = zeros(connectionGraph_.nNodes * connectionGraph_.dimConnection, connectionGraph_.nNodes * connectionGraph_.dimConnection);

            for i = 1:connectionGraph_.nNodes
                for j = 1:connectionGraph_.nNodes
                    if connectionGraph_.adjacencyMatrix(i,j) == 1
                        %fprintf(strcat("fromNode: ", num2str(i), " toNode: ", num2str(j), ".\n"))
                        connectionGraph_.connectionAdjacency((i - 1)*d + 1:( i*d ), (j - 1)*d + 1:( j*d )) = connection_(i,j,:,:);
                        connectionGraph_.connectionLaplacian((i - 1)*d + 1:( i*d ), (j - 1)*d + 1:( j*d )) = -connection_(i,j,:,:);
                    elseif i==j
                        connectionGraph_.connectionLaplacian((i - 1)*d + 1:( i*d ), (j - 1)*d + 1:( j*d )) = connectionGraph_.degreeVector(i) * eye(d);
                    end
                end
            end 

            plainIncidence = full(connectionGraph_.graphObject.incidence).';
            
            for e=1:connectionGraph_.nEdges
                for i=1:connectionGraph_.nNodes
                    if plainIncidence(e,i) == 1
                        connectionGraph_.connectionIncidence((e - 1)*d + 1:( e*d ), (i - 1)*d + 1:( i*d )) = eye(d);
                    elseif plainIncidence(e,i) == -1
                        head = 0;
                        for j=1:connectionGraph_.nNodes
                            if plainIncidence(e,j) == 1
                                head = j;
                            end
                        end
                        connectionGraph_.connectionIncidence((e - 1)*d + 1:( e*d ), (i - 1)*d + 1:( i*d )) = -connection_(head, i, :, :);
                    end
                    
                end
            end
            
            [evecs, evals] = schur(connectionGraph_.connectionLaplacian);
            [~, ind] = sort(diag(evals),'ComparisonMethod','abs');
            connectionGraph_.connectionEigenvalues = diag(evals(ind,ind));
            connectionGraph_.connectionEigenvectors = evecs(:,ind);
            
            connectionGraph = connectionGraph_;
            connectionGraph.connection = connection_;
        end

        function connectionGraph = setEdgeConnection(connectionGraph_, fromNode, toNode, connectionMatrix)
            
            newConnection = connectionGraph_.connection;
            newConnection(fromNode,toNode,:,:) = connectionMatrix;
            newConnection(toNode,fromNode,:,:) = connectionMatrix.';

            connectionGraph = setConnectionMatrices(connectionGraph_, newConnection);
        end

    end

end

