function fig3bc

% simulate the reaction-diffusion equation with fixed domain size and 
% fixed cell diameter in x direction, then calculate the gradient variability
% and positional error as the diameter in y direction increases.

% options
simulate = true; % if false, plot results from saved files instead of generating new data
writeresults = true; % when creating new data, should the results be written to output files?
prefix = 'data/'; % filename prefix
fitcosh = true; % fit an exp or cosh to the gradients
checkbounds = true; % make sure fitted lambda and C0 are positive
plotresults = true; % plot resulting data or not
zerofluxbc = true; % use zero-flux boundary conditions to solve the reaction-diffusion equation, otherwise periodic boundary conditions
average = false; % calculate average gradient variability and positional error across all gradients
LineWidth = 1;
FontSize = 18;

% parameters
tol = 1e-10; % numerical tolerance for solver and fitting
nruns = 1e3; % number of independent simulation runs
nboot = 1e4; % number of bootstrap samples for error estimation
res = 3; % linear grid resolution of each cell
diameter = 5; % default cell diameter [µm]
diameterX = diameter; % cell diameter in x direction
diameterY = 1:20; % range of cell diameters in y direction
mu_D = 0.033; % mean morphogen diffusion constant [µm^2/s]
mu_lambda = 20; % mean gradient length [µm]
mu_d = mu_D/mu_lambda^2; % mean morphogen degradation rate [1/s]
mu_p = mu_d; % mean morphogen production rate [substance/(µm^3*s)]
CV = 0.3; % coefficient of variation of the kinetic parameters
ncS = 5; % number of cells along the source domain
ncP = 50; % number of cells along the patterning domain
fitrange = 1:10; % ncY index range for curve fitting
readouts = [3 6 9]; % readout positions in units of mu_lambda
colors = [0, 0.4470, 0.7410;...
          0.8500, 0.3250, 0.0980;...
          0.9290, 0.6940, 0.1250]; % plot colors for the different readout positions


% derived variables
ncX = ncS + ncP; % total number of cells along the patterning axis
LS = ncS * diameterX; % source length
LP = ncP * diameterX; % pattern length
hX = diameterX / res; % grid spacing in x direction

% analytical deterministic solution
C = @(x) mu_p/mu_d * ((x<0) .* (1-cosh(x/mu_lambda)) + sinh(LS/mu_lambda) / sinh((LS+LP)/mu_lambda) * cosh((LP-x)/mu_lambda));

% log-normal distribution with adjusted mean & CV
logndist = @(mu, CV) makedist('Lognormal', 'mu', log(mu/sqrt(1+CV^2)), 'sigma', sqrt(log(1+CV^2)));

CVfun = @(x) nanstd(x) / nanmean(x);
SEfun = @(x) nanstd(x) / sqrt(sum(~isnan(x)));

fitopt = statset('TolFun', tol, 'TolX', tol);


if simulate
    CV_lambda = NaN(numel(diameterY), 1);
    CV_lambda_SE = CV_lambda;
    CV_0 = CV_lambda;
    CV_0_SE = CV_lambda;
    sigma_x = NaN(numel(diameterY), numel(readouts));
    sigma_x_SE = sigma_x;
    
    % loop over domain widths
    tic
    for i = 1:length(diameterY)
        diameterY(i)

        LY = 0.2*LP; % fixed width of the domain

        % if the last cell exceeds the domain boundary, determine
        % whether or not the domain width is closer to the target
        % without it. Remove the last cell if the domain width is 
        % closer to the target without it, keep it otherwise.
        if mod(LY/diameterY(i),1) == 0
            ncY = LY/diameterY(i);
        elseif mod(LY/diameterY(i),1) < 0.5
            ncY = floor(LY/diameterY(i));
            LY = ncY * diameterY(i);
        else
            ncY = ceil(LY/diameterY(i));
            LY = ncY * diameterY(i);
        end
        
        hY = diameterY(i) / res; % grid spacing in y direction

        lambda = NaN(nruns, ncY);
        C0 = lambda;
        x_theta = NaN(nruns, ncY, numel(readouts));

        % loop over several independent runs
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
            x = linspace(-LS, LP-hX, ncX * res)' + hX/2;
            C_old = C(x) * ones(1, ncY * res);
            
            param = {p,d,D};
            % solve the reaction-diffusion equation
            if zerofluxbc
                C_new = boundary_conditions.zeroflux_bc(C_old, param, hX, hY, tol);
            else
                C_new = boundary_conditions.periodic_bc(C_old, param, hX, hY, tol);
            end

            % for each cell in y direction
            for y = 1:ncY
                % determine the cellular readout concentrations
                yidx = (y-1)*res+1:y*res;
                C_readout = mean(C_new(:,yidx), 2);
                
                % fit an exponential in log space in the patterning domain
                param = polyfit(x, log(C_readout), 1);
                lambda(k,y) = -1/param(1);
                C0(k,y) = exp(param(2));
        
                % fit a hyperbolic cosine in log space in the patterning domain
                if fitcosh
                    logcosh = @(p,x) p(2) + log(cosh((LP-x)/p(1)));
                    mdl = fitnlm(x, log(C_readout), logcosh, [lambda(k,y) log(C0(k,y)) - log(cosh(LP/lambda(k,y)))], 'Options', fitopt);
                    lambda(k,y) = mdl.Coefficients.Estimate(1);
                    C0(k,y) = exp(mdl.Coefficients.Estimate(2)) * cosh(LP/lambda(k,y));
                end

                % discard out-of-bound values
                if checkbounds
                    if lambda(k,y) <= 0 || C0(k,y) <= 0
                        lambda(k,y) = NaN;
                        C0(k,y) = NaN;
                    end
                end
            
                % loop over readout positions
                for j = 1:numel(readouts)
            
                    % threshold concentration
                    C_theta = C(readouts(j) * mu_lambda);
    
                    % determine the readout position
                    indices = find(diff(sign(C_readout - C_theta)));
                    
                    x_theta_all = [];
                    for idx = indices
                        x_theta_all = [x_theta_all, interp1(C_readout([idx idx+1]), x([idx idx+1]), C_theta)];
                    end
                    x_theta(k,y,j) = mean(x_theta_all);
                end
            end
        end

        % determine the CV of the decay length and the amplitude, and the positional error, over the independent runs
        % and also their standard errors from bootstrapping
        y = ceil(ncY/2); % a single cell-row at the middle of the domain
        CV_lambda(i) = CVfun(lambda(:,y));
        CV_lambda_SE(i) = nanstd(bootstrp(nboot, CVfun, lambda(:,y)));
        CV_0(i) = CVfun(C0(:,y));
        CV_0_SE(i) = nanstd(bootstrp(nboot, CVfun, C0(:,y)));
        for j = 1:numel(readouts)
            sigma_x(i,j) = nanstd(x_theta(:,y,j));            
            sigma_x_SE(i,j) = nanstd(bootstrp(nboot, @nanstd, x_theta(:,y,j)));
        end
        
        % average CV's and sigma_x across all gradients
        if average
            CV_lambda(i) = nanmean(nanstd(lambda,0,1)./nanmean(lambda,1));
            CV_lambda_SE(i) = nanstd(bootstrp(nboot, @(x) nanmean(nanstd(x,0,1)./nanmean(x,1)), lambda));
            CV_0(i) = nanmean(nanstd(C0,0,1)./nanmean(C0,1));
            CV_0_SE(i) = nanstd(bootstrp(nboot, @(x) nanmean(nanstd(x,0,1)./nanmean(x,1)), C0));
            for j = 1:numel(readouts)
                sigma_x(i,j) = nanmean(nanstd(x_theta(:,:,j),0,1));
                sigma_x_SE(i,j) = nanstd(bootstrp(nboot, @(x) nanmean(nanstd(x,0,1)), x_theta(:,:,j)));
            end
        end

    end
    toc
    
    % write results to output files
    if writeresults
        for j = 1:numel(readouts)
            T = table();
            T.diameterY = diameterY';
            T.CV_lambda = CV_lambda;
            T.CV_lambda_SE = CV_lambda_SE;
            T.CV_0 = CV_0;
            T.CV_0_SE = CV_0_SE;
            T.sigma_x = sigma_x(:,j);
            T.sigma_x_SE = sigma_x_SE(:,j);
            if zerofluxbc
                writetable(T, [prefix 'readout_' num2str(readouts(j)) 'lambda_zeroflux.csv']);
            else
                writetable(T, [prefix 'readout_' num2str(readouts(j)) 'lambda_periodic.csv']);
            end
        end
    end
else
    % read existing data
    for j = 1:numel(readouts)
        if zerofluxbc
            T = readtable([prefix 'readout_' num2str(readouts(j)) 'lambda_zeroflux.csv']);
        else
            T = readtable([prefix 'readout_' num2str(readouts(j)) 'lambda_periodic.csv']);
        end
        diameterY = T.diameterY';
        CV_lambda = T.CV_lambda;
        CV_lambda_SE = T.CV_lambda_SE;
        CV_0 = T.CV_0;
        CV_0_SE = T.CV_0_SE;
        sigma_x(:,j) = T.sigma_x;
        sigma_x_SE(:,j) = T.sigma_x_SE;
    end
end

% plot results
if plotresults
    close all
    figure
       
    for j = 1:numel(readouts)
        name = ['x_\theta = ' num2str(readouts(j)) ' \mu_\lambda'];
        
        hold on
        mdl = fitnlm(diameterY(fitrange)./diameterX, sigma_x(fitrange,j) / diameterX, @(p,x) p(1)*sqrt(x), 1);
        plot(diameterY./diameterX, feval(mdl, diameterY./diameterX), '-', 'HandleVisibility', 'off', 'LineWidth', LineWidth, 'Color', colors(j,:))
        errorbar(diameterY./diameterX, sigma_x(:,j)./diameterX, sigma_x_SE(:,j)./diameterX, 'o', 'DisplayName', name, 'LineWidth', LineWidth, 'Color', colors(j,:))
        xlabel('\delta_y / \delta_x')
        ylabel('Relative Positional error  \sigma_x/\delta_x (cells)')
        legend('show')
        set(gca, 'LineWidth', LineWidth, 'FontSize', FontSize, 'XScale', 'log', 'YScale', 'log')
        box on    
        xlim([min(diameterY./diameterX) max(diameterY./diameterX)])
    end
end

end