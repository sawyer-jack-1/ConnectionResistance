classdef GraphXManOpt

    properties (SetAccess = private)
    end

    methods(Static)

        function x = getBestAlignment(U, V, verbose)
            
            dimAlign = size(U,2);
            
            % Create the problem structure.
            manifold = stiefelfactory(dimAlign, dimAlign, 1);
            problem.M = manifold;

            % Define the problem cost function and its gradient.
            problem.cost  = @(x) trace( (U.' - x * V.').' * (U.' - x * V.'));
            problem.egrad = @(x) -2 * (U.' - x * V.') * V;
            problem.ehess = @(x, xdot) V.' * V;
            
            if verbose == 1
                % Numerically check gradient and Hessian consistency.
                figure;
                checkgradient(problem);
                figure;
                checkhessian(problem);
            end

            % Solve. 
            options.verbosity = 1;
            [x, xcost, info] = trustregions(problem, [], options);          %#ok<ASGLU>

            if verbose == 1
                % Display some statistics.
                figure;
                semilogy([info.iter], [info.gradnorm], '.-');
                xlabel('Iteration #');
                ylabel('Gradient norm');
                title('Convergence of the trust-regions algorithm on the sphere');
            end
            
        end

    end

    methods

    end

end