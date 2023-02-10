classdef GraphXFunction

    properties (SetAccess = private)
        nNodes;         % # nodes in parent graph
        parentGraph;    % parent graphX object
        nodeData;       % function values at each node
        fourierData;    % coefficients of function in Laplacian eigenbasis
    end

    methods(Static)

    end

    methods
        
        %% Simple constructor method that auto-populates either fourier data or node data
        function graphFunction = GraphXFunction(domain, data, parentGraph)
            graphFunction.parentGraph = parentGraph;
            graphFunction.nNodes = parentGraph.nNodes;
            if domain == "nodes"
                graphFunction.nodeData = data;
                graphFunction.fourierData = parentGraph.eigenvectors.' * data;
            elseif domain == "fourier" 
                graphFunction.fourierData = data;
                graphFunction.nodeData = parentGraph.eigenvectors * data;
            else
                fprintf("<ERROR:GraphXFunction> Invalid node data specification method. Acceptable options are 'nodes' and 'fourier'.\n ")
            end
        end

    end

end