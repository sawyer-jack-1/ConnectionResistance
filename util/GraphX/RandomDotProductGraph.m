classdef RandomDotProductGraph < GraphX

    properties (SetAccess = private)
        dimLatentPositions;     % Dimension of space from which latent positions are sampled
        samplingMethod;         % Distribution used to sample the latent positions; "Square" or "Hypersphere" supported
        sparsityFactor;         % Uniform weight applied to the probability of an edge occuring between two nodes
        latentPositions;        % Matrix of latent positions
        probMatrix;             % Dot product matrix scaled by sparisity factor, each entry is P(i ~ j)
    end

    methods(Static)

        %% Generate the latent position matrix from one of the supported distributions
        % Square: Uniformly sampled from a hypercube with interior in the
        % positive octant with appropriate size
        % Hypersphere: Uniformly sample from the d-dimensional unit hyphersphere
        % within the positive octant
        function X = generateLatentPositions(n,d, method)
            if method=="Square" 
                X = unifrnd(0,1/sqrt(2),n,d); 
            elseif method=="Hypersphere"
                X = abs(sampleSphereSurface(n,d,1)); %Uniformly sample from hypersphere surface
            else 
                fprintf("<ERROR:RandomDotProductGraph> Invalid latent positioning sampling <method>. Cannot calculate latent positions.\n")
                return;
            end
        end

        %% Calculate the ASE U-statistic with given kernel method and parameter sigma
        % Supported kernels: "Gaussian" and "Laplace"
        % Default sigma = 1
        function U = uStatistic(randomGraph1, randomGraph2, kernelMethod, sigma)
            
            d1 = randomGraph1.dimLatentPositions;
            d2 = randomGraph2.dimLatentPositions;
            if d1 ~= d2
                fprintf("<ERROR:RandomDotProductGraph> Incompatible latent dimensions. U Statistic cannot be calculated.\n");
                return;
            end
            
            U = 0;
           
            n1 = randomGraph1.nNodes;
            ase1 = randomGraph1.getAdjacencySpectralEmbedding();
            
            n2 = randomGraph2.nNodes;
            ase2 = randomGraph2.getAdjacencySpectralEmbedding();
            
            function k = kernel(x1,x2)
                if kernelMethod == "Gaussian"
                    k = exp(-1 ./ sigma^2 *(norm(x1-x2,2)^2));
                elseif kernelMethod == "Laplace"
                    k = exp(-1 ./ sigma *(norm(x1-x2,1)));
                else
                    fprintf("<ERROR:RandomDotProductGraph> Unrecognized kernel <method>. U Statistic cannot be calculated.\n");
                    return;
                end
            end

            for i1=1:n1
                for i2=1:n1
                    if i1 ~= i2
                           U = U + 1./(n1*(n1-1)) * kernel(ase1(i1,:),ase1(i2,:));
                    end
                end
            end
            
            for j1=1:n2
                for j2=1:n2
                    if j1 ~= j2
                           U = U + 1./(n2*(n2-1)) * kernel(ase2(j1,:),ase2(j2,:));
                    end
                end
            end

            for i=1:n1
                for j=1:n2
                    U = U - 2./(n1*n2) * kernel(ase1(i,:),ase2(j,:));
                end
            end

        end

        function alignment_matrix = getAlignmentMatrix(inputGraphX1, inputGraphX2, spectralDepth)
            
            initialEigenvectors = inputGraphX1.eigenvectors(:,1:spectralDepth);
            terminalEigenvectors = inputGraphX2.eigenvectors(:,1:spectralDepth);
            
            alignment_matrix = GraphXManOpt.getBestAlignment(initialEigenvectors, terminalEigenvectors, 0);
            
        end

    end
    
    methods

        %% Constructor method for the RDPG with specified dimension, size, sparsity factor,
        % and method.
        function randomGraph = RandomDotProductGraph(n, d, s, method, inputPositions)

            if isempty(inputPositions)
                X = RandomDotProductGraph.generateLatentPositions(n,d,method);
            else
                X = inputPositions;
            end
            probabilityMatrix = s * (X * X.');
            A = zeros(n,n);
            for i = 1:n
                for j = i:n
                    if i~=j  
                        A(i,j) = binornd(1,probabilityMatrix(i,j));
                        A(j,i) = A(i,j);
                    end
                end
            end
            randomGraph = randomGraph@GraphX(A);
            randomGraph.sparsityFactor = s;
            randomGraph.probMatrix = probabilityMatrix;
            randomGraph.samplingMethod = method;
            randomGraph.dimLatentPositions = d;
            randomGraph.latentPositions = X;

        end

        %% Calculate the adjacency spectral embedding from a RDPG
        function ase = getAdjacencySpectralEmbedding(randomGraph)
            [evecs,evals] = getAdjacencyEigenvectors(randomGraph);
            d = randomGraph.dimLatentPositions;
            D = diag(evals(1:d));
            V = evecs(:,1:d);
            ase = V * sqrt(abs(D));
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

