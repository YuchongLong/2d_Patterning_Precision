classdef boundary_conditions
    % functions for different boundary conditions
    
    methods (Static)

        %% function for zero-flux boundary condition

        function C_new = zeroflux_bc(C_old, param, hx, hy, tol)
            
            % C_old: initial concentration
            % param: kinetic parameters {p, d, D}
            %       p: production rate
            %       d: degradation rate
            %       D: diffusion rate
            % hx: grid spacing in x direction
            % hy: grid spacing in y direction
            % tol: error tolerance

            p = param{1};
            d = param{2};
            D = param{3};
            
            % solve the reaction-diffusion equation
            err = inf;
            Cx = NaN(size(C_old));
            Cy = Cx;
            while err > tol
                % zero-flux boundary conditions on all sides
                Cx(1,:) = 2 * C_old(2,:);
                Cx(end,:) = 2 * C_old(end-1,:);
                Cy(:,1) = 2 * C_old(:,2);
                Cy(:,end) = 2 * C_old(:,end-1);

                % 5-point central difference stencil
                Cx(2:end-1,:) = C_old(1:end-2,:) + C_old(3:end,:);
                Cy(:,2:end-1) = C_old(:,1:end-2) + C_old(:,3:end);
                C_new = (hx^2 * hy^2 * p + D .* (Cx * hy^2 + Cy * hx^2)) ./ (hx^2 * hy^2 * d + 2 * D * (hx^2 + hy^2));
                err = max(abs(C_new - C_old), [], 'all');
                C_old = C_new;
            end

        end

        %% function for periodic boundary conditions

        function C_new = periodic_bc(C_old, param, hx, hy, tol)
            
            % C_old: initial concentration
            % param: kinetic parameters {p, d, D}
            %       p: production rate
            %       d: degradation rate
            %       D: diffusion rate
            % hx: grid spacing in x direction
            % hy: grid spacing in y direction
            % tol: error tolerance

            p = param{1};
            d = param{2};
            D = param{3};
            
            % solve the reaction-diffusion equation
            err = inf;
            Cx = NaN(size(C_old));
            Cy = Cx;
            while err > tol
                % zero-flux boundary conditions in x direction
                Cx(1,:) = 2 * C_old(2,:);
                Cx(end,:) = 2 * C_old(end-1,:);
                
                % periodic boundary conditions in y direction
                Cy(:,1) = C_old(:,2) + C_old(:,end);
                Cy(:,end) = C_old(:,end-1) + C_old(:,1);
                
                % 5-point central difference stencil
                Cx(2:end-1,:) = C_old(1:end-2,:) + C_old(3:end,:);
                Cy(:,2:end-1) = C_old(:,1:end-2) + C_old(:,3:end);
                C_new = (hx^2 * hy^2 * p + D .* (Cx * hy^2 + Cy * hx^2)) ./ (hx^2 * hy^2 * d + 2 * D * (hx^2 + hy^2));
                err = max(abs(C_new - C_old), [], 'all');                
                C_old = C_new;
            end

        end


    end


end