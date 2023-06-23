function fig4cd

% simulate the reaction-diffusion equation with fixed patterning domain
% length, but increasing cell diameter. 

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
diameter = 1:40; % cell diameters to loop over [µm]
LS = 25; % source length [µm]
LP = 250; % pattern length [µm]
mu_D = 0.033; % mean morphogen diffusion constant [µm^2/s]
mu_lambda = 20; % mean gradient length [µm]
mu_d = mu_D/mu_lambda^2; % mean morphogen degradation rate [1/s]
mu_p = mu_d; % mean morphogen production rate [substance/(µm^3*s)]
CV = 0.3; % coefficient of variation of the kinetic parameters
ncY = [1 4 10 40]; % number of cells in transverse direction to loop over
readouts = 9; % readout positions in units of mu_lambda
colors = [0, 0.4470, 0.7410;...
          0.8500, 0.3250, 0.0980;...
          0.9290, 0.6940, 0.1250;...
          0, 0.5, 0]; % plot colors for the different readout positions

% analytical deterministic solution
C = @(x) mu_p/mu_d * ((x<0) .* (1-cosh(x/mu_lambda)) + sinh(LS/mu_lambda) / sinh((LS+LP)/mu_lambda) * cosh((LP-x)/mu_lambda));


% log-normal distribution with adjusted mean & CV
logndist = @(mu, CV) makedist('Lognormal', 'mu', log(mu/sqrt(1+CV^2)), 'sigma', sqrt(log(1+CV^2)));

CVfun = @(x) nanstd(x) / nanmean(x);
SEfun = @(x) nanstd(x) / sqrt(sum(~isnan(x)));

fitopt = statset('TolFun', tol, 'TolX', tol);


if simulate
    CV_lambda = NaN(numel(diameter), 1);
    CV_lambda_SE = CV_lambda;
    CV_0 = CV_lambda;
    CV_0_SE = CV_lambda;
    sigma_x = NaN(numel(diameter), 1);
    sigma_x_SE = sigma_x;
    
    % loop over domain widths
    
    for i = 1:numel(ncY)

        tic
        for j = 1:numel(diameter)
            diameter(j)
            
            % approximate the number of cells in the morphogen source
            if mod(LS/diameter(j),1) == 0
                ncS = LS/diameter(j);
            elseif mod(LS/diameter(j),1) < 0.5
                ncS = floor(LS/diameter(j));
                LS = ncS * diameter(j);
            else
                ncS = ceil(LS/diameter(j));
                LS = ncS * diameter(j);
            end
            
            % approximate the number of cells in the patterning domain
            if mod(LP/diameter(j),1) == 0
                ncP = LP/diameter(j);
            elseif mod(LP/diameter(j),1) < 0.5
                ncP = floor(LP/diameter(j));
                LP = ncP * diameter(j);
            else
                ncP = ceil(LP/diameter(j));
                LP = ncP * diameter(j);
            end

            % total number of cells along the patterning axis
            ncX = ncP + ncS; 
            
            % grid spacing
            h = diameter(j) / res;
    
            lambda = NaN(nruns, ncY(i));
            C0 = lambda;
            x_theta = NaN(nruns, ncY(i));
    
            % loop over several independent runs
            for k = 1:nruns

                k
    
                % draw random kinetic parameters for each cell
                p = random(logndist(mu_p, CV), ncX, ncY(i));
                d = random(logndist(mu_d, CV), ncX, ncY(i));
                D = random(logndist(mu_D, CV), ncX, ncY(i));
                p(ncS+1:end,:) = 0; % no production in the patterning domain
        
                % inflate them to the grid
                p = repelem(p, res, res);
                d = repelem(d, res, res);
                D = repelem(D, res, res);
        
                % deterministic solution as initial guess
                x = linspace(-LS, LP-h, ncX * res)' + h/2;
                C_old = C(x) * ones(1, ncY(i) * res);
                
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
                for y = 1:ncY(i)
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
                
                    % threshold concentration at readout position
                    C_theta = C(readouts * mu_lambda);
    
                    % determine the readout position
                    indices = find(diff(sign(C_readout - C_theta)));
                    
                    x_theta_all = [];
                    if ~isempty(indices)
                        for idx = indices
                            x_theta_all = [x_theta_all, interp1(C_readout([idx idx+1]), x([idx idx+1]), C_theta)];
                        end
                    end
                    x_theta(k,y) = mean(x_theta_all);
                end
            end
    
            % determine the CV of the decay length and the amplitude, and the positional error, over the independent runs
            % and also their standard errors from bootstrapping
            CV_lambda(j) = nanmean(nanstd(lambda,0,1)./nanmean(lambda,1));
            CV_lambda_SE(j) = nanstd(bootstrp(nboot, @(x) nanmean(nanstd(x,0,1)./nanmean(x,1)), lambda));
            CV_0(j) = nanmean(nanstd(C0,0,1)./nanmean(C0,1));
            CV_0_SE(j) = nanstd(bootstrp(nboot, @(x) nanmean(nanstd(x,0,1)./nanmean(x,1)), C0));
            sigma_x(j) = nanmean(nanstd(x_theta(:,:),0,1));
            sigma_x_SE(j) = nanstd(bootstrp(nboot, @(x) nanmean(nanstd(x,0,1)), x_theta(:,:)));
        end

        % write results to output files
        if writeresults
            T = table();
            T.diameter = diameter';
            T.CV_lambda = CV_lambda;
            T.CV_lambda_SE = CV_lambda_SE;
            T.CV_0 = CV_0;
            T.CV_0_SE = CV_0_SE;
            T.sigma_x = sigma_x;
            T.sigma_x_SE = sigma_x_SE;
            writetable(T, [prefix 'readout_' num2str(readouts) 'lambda_zeroflux_' num2str(ncY(i)) 'ncY.csv']);
        end
        toc
    end

else
    % read existing data
    for i = 1:numel(ncY)
        T = readtable([prefix 'readout_' num2str(readouts) 'lambda_zeroflux_' num2str(ncY(i)) 'ncY.csv']);
        diameter = T.diameter;
        CV_lambda = T.CV_lambda;
        CV_lambda_SE = T.CV_lambda_SE;
        CV_0 = T.CV_0;
        CV_0_SE = T.CV_0_SE;
        sigma_x(:,i) = T.sigma_x;
        sigma_x_SE(:,i) = T.sigma_x_SE;
    end
end

% plot results
if plotresults
    close all
    figure('Position', [0 0 1100 500]);
    
    for i = 1:numel(ncY)
        name = ['N_y = ' num2str(ncY(i)) ' cells'];

        subplot(1,2,1)
        hold on
        plot(diameter,diameter,'--','HandleVisibility','off','LineWidth', LineWidth, 'Color', 'k')
        errorbar(diameter, sigma_x(:,i), sigma_x_SE(:,i), 'o', 'DisplayName', name, 'LineWidth', LineWidth, 'Color', colors(i,:))
        annotation('textbox', [.32 .5 .3 .3], 'String','\sigma_x = \delta_x', 'FitBoxToText', 'on', 'FontSize', FontSize);
        xlabel('\delta_x [\mum]')
        ylabel('Absolute positional error  \sigma_x [\mum]')
        legend('show')
        set(gca, 'LineWidth', LineWidth, 'FontSize', FontSize)
        box on    
        xlim([min(diameter) max(diameter)])
        ylim([0 25])

        subplot(1,2,2)
        hold on
        errorbar(diameter, sigma_x(:,i)./diameter, sigma_x_SE(:,i)./diameter, 'o', 'DisplayName', name, 'LineWidth', LineWidth, 'Color', colors(i,:))
        yline(1,'--','HandleVisibility','off','LineWidth', LineWidth, 'Color', 'k')
        xlabel('\delta_x [\mum]')
        ylabel('Relative positional error  \sigma_x/\delta_x [cells]')
        legend('show')
        set(gca, 'LineWidth', LineWidth, 'FontSize', FontSize)
        box on    
        xlim([0 max(diameter)])
        
        
    end
end

end