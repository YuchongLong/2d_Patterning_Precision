function fig4


% options
simulate = true; % if false, plot results from saved files instead of generating new data
writeresults = true; % when creating new data, should the results be written to output files?
prefix = 'data/'; % filename prefix
fitcosh = true; % fit an exp or cosh to the gradients
checkbounds = true; % make sure fitted lambda and C0 are positive
plotresults = true; % plot resulting data or not
LineWidth = 1;
FontSize = 18;

% parameters
tol = 1e-10; % numerical tolerance for solver and fitting
nruns = 1e3; % number of independent simulation runs
nboot = 1e4; % number of bootstrap samples for error estimation
res = 3; % linear grid resolution of each cell
mu_D = 0.033; % mean morphogen diffusion constant [µm^2/s]
CV = 0.3; % coefficient of variation of the kinetic parameters
ncY = 10; % number of cells in transverse direction to loop over
LP = 50:25:300; % range of patterning domain lengths 
LS = 0.16 * LP; % morphogen source length (fit to data from O. Wartlick et al. 2011)
mu_lambda = 0.11 * LP;  % formula from O. Wartlick et al. 2011, mu_lambda increases with the patterning domain length 
diameter = linspace(4.5,1.5,length(LP));  % cell diameter [µm] (decreases linearly over time)


% log-normal distribution with adjusted mean & CV
logndist = @(mu, CV) makedist('Lognormal', 'mu', log(mu/sqrt(1+CV^2)), 'sigma', sqrt(log(1+CV^2)));
fitopt = statset('TolFun', tol, 'TolX', tol);


if simulate
    sigma_x = NaN(numel(LP));
    sigma_x_SE = sigma_x;
    
    % loop over domain widths
    tic
    for i = 1:length(LP)
        i
        
        % derived variables
        mu_d = mu_D/mu_lambda(i)^2; % mean morphogen degradation rate [1/s]
        mu_p = mu_d; % mean morphogen production rate [substance/(µm^3*s)] 
        ncP = round(LP(i)/diameter(i)); % number of cells along the patterning domain
        ncS = round(LS(i)/diameter(i)); % number of cells along the source domain
        ncX = ncS + ncP; % total number of cells along the patterning axis
        h = diameter(i) / res; % grid spacing

        readouts = 0.4 * LP(i); % readout position at 0.4 of LP 

        % analytical deterministic solution
        C = @(x) mu_p/mu_d * ((x<0) .* (1-cosh(x/mu_lambda(i))) + sinh(LS(i)/mu_lambda(i)) / sinh((LS(i)+LP(i))/mu_lambda(i)) * cosh((LP(i)-x)/mu_lambda(i)));

        % loop over several independent runs
        lambda = NaN(nruns, ncY);
        C0 = lambda;
        x_theta = NaN(nruns, ncY);
        for k = 1:nruns

            % draw random kinetic parameters for each cell
            p = random(logndist(mu_p, CV), ncX, ncY);
            d = random(logndist(mu_d, CV), ncX, ncY);
            D = random(logndist(mu_D, CV), ncX, ncY);
            p(ncS+1:end,:) = 0; % no production in the patterning domain
    
            % inflate them to the grid
            p = repelem(p, res, res);
            d = repelem(d, res, res);
            D = repelem(D, res, res);
    
            % deterministic solution as initial guess
            x = linspace(-LS(i), LP(i)-h, ncX * res)' + h/2;
            C_old = C(x) * ones(1, ncY * res);
            
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
                C_new = (h^2 * p + D .* (Cx + Cy)) ./ (4 * D + h^2 * d);
                err = max(abs(C_new - C_old), [], 'all');
                C_old = C_new;
            end

            % for each cell in y direction
            for y = 1:ncY
                % determine the cellular readout concentrations
                yidx = (y-1)*res+1:y*res;
                C_readout = mean(C_new(:,yidx), 2);
                
                % fit an exponential in log space in the patterning domain
                x = linspace(-LS(i), LP(i)-h, ncX)' + h/2;
                param = polyfit(x, log(C_readout), 1);
                lambda(k,y) = -1/param(1);
                C0(k,y) = exp(param(2));
        
                % fit a hyperbolic cosine in log space in the patterning domain
                if fitcosh
                    logcosh = @(p,x) p(2) + log(cosh((LP(i)-x)/p(1)));
                    mdl = fitnlm(x, log(C_readout), logcosh, [lambda(k,y) log(C0(k,y)) - log(cosh(LP(i)/lambda(k,y)))], 'Options', fitopt);
                    lambda(k,y) = mdl.Coefficients.Estimate(1);
                    C0(k,y) = exp(mdl.Coefficients.Estimate(2)) * cosh(LP(i)/lambda(k,y));
                end

                % discard out-of-bound values
                if checkbounds
                    if lambda(k,y) <= 0 || C0(k,y) <= 0
                        lambda(k,y) = NaN;
                        C0(k,y) = NaN;
                    end
                end
            
                % threshold concentration
                C_theta = C(readouts);

                % determine the readout position
                indices = find(diff(sign(C_readout - C_theta)));
                if isempty(indices)
                    continue
                end
                
                x_theta_all = [];
                for idx = indices
                    x_theta_all = [x_theta_all, interp1(C_readout([idx idx+1]), x([idx idx+1]), C_theta)];
                end
                x_theta(k,y) = mean(x_theta_all);

            end
        end

        % determine the average positional error across all cell-rows, over the independent runs
        % and also their standard errors from bootstrapping
        sigma_x(i) = nanstd(nanmean(x_theta,2));
        sigma_x_SE(i) = nanstd(bootstrp(nboot, @(x) nanstd(nanmean(x,2)), x_theta));
    end
    toc
    
    % write results to output files
    if writeresults
        T = table();
        T.LP = LP';
        T.sigma_x = sigma_x;
        T.sigma_x_SE = sigma_x_SE;
        writetable(T, [prefix 'wing_disc_positional_error.csv']);
    end
else
    % read existing data
    T = readtable([prefix 'wing_disc_positional_error.csv']);
    LP = T.LP';
    sigma_x(:,:) = T.sigma_x;
    sigma_x_SE(:,:) = T.sigma_x_SE;
end

% plot results
if plotresults
    close all
    figure

    errorbar(LP, sigma_x, sigma_x_SE, 'LineWidth', LineWidth)
    xlabel('Patterning length L_p [\mum]')
    ylabel('Positional error \sigma_x [\mum]')
    set(gca, 'LineWidth', LineWidth, 'FontSize', FontSize)
    xlim([0 300])
    grid on

end

end