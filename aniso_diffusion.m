function aniso_diffusion

% simulate the reaction-diffusion equation with anisotropic diffusion
% constants in x and y directions.

% options
simulate = true; % if false, plot results from saved files instead of generating new data
writeresults = true; % when creating new data, should the results be written to output files?
prefix = 'data/'; % filename prefix
fitcosh = true; % fit an exp or cosh to the gradients
checkbounds = true; % make sure fitted lambda and C0 are positive
plotresults = true; % plot resulting data or not
vary_Dx = false; % if false, then fix Dx, vary Dy (Fig.S4), otherwise fix Dy, vary Dx
LineWidth = 1;
FontSize = 18;

% parameters
tol = 1e-10; % numerical tolerance for solver and fitting
nruns = 1e3; % number of independent simulation runs
nboot = 1e4; % number of bootstrap samples for error estimation
res = 3; % linear grid resolution of each cell
diameter = 5; % cell diameter [µm]
mu_D = 0.033; % mean morphogen diffusion constant [µm^2/s]
mu_lambda = 20; % mean gradient length [µm]
mu_d = mu_D/mu_lambda^2; % mean morphogen degradation rate [1/s]
mu_p = mu_d; % mean morphogen production rate [substance/(µm^3*s)]
CV = 0.3; % coefficient of variation of the kinetic parameters
ncS = 5; % number of cells along the source domain
ncP = 50; % number of cells along the patterning domain
ncY = [1 10 40]; % number of cells in transverse direction to loop over
fitrange = 1:11; % index range for curve fitting
readouts = 6; % readout positions in units of mu_lambda
colors = [0, 0.4470, 0.7410;...
          0.8500, 0.3250, 0.0980;...
          0.9290, 0.6940, 0.1250]; % plot colors for the different readout positions

alpha = logspace(-1,1,11); % degree of orthotropy in morphogen transport, alpha = Dy/Dx

% derived variables
ncX = ncS + ncP; % total number of cells along the patterning axis
LS = ncS * diameter; % source length
LP = ncP * diameter; % pattern length
h = diameter / res; % grid spacing

% log-normal distribution with adjusted mean & CV
logndist = @(mu, CV) makedist('Lognormal', 'mu', log(mu/sqrt(1+CV^2)), 'sigma', sqrt(log(1+CV^2)));

CVfun = @(x) nanstd(x) / nanmean(x);
SEfun = @(x) nanstd(x) / sqrt(sum(~isnan(x)));

fitopt = statset('TolFun', tol, 'TolX', tol);


if simulate
    
    % loop over domain widths
    for i = 1:numel(ncY)
        ncY(i)

        CV_lambda = NaN(numel(alpha), 1);
        CV_lambda_SE = CV_lambda;
        CV_0 = CV_lambda;
        CV_0_SE = CV_lambda;
        sigma_x = NaN(numel(alpha), 1);
        sigma_x_SE = sigma_x;
        
        % loop over alpha
        tic
        for j = 1:numel(alpha)
            j
    
            lambda = NaN(nruns, ncY(i));
            C0 = lambda;
            x_theta = NaN(nruns, ncY(i));
            
            if vary_Dx
                % if Dx varies, mu_D and mu_lambda need to be adjusted
                % accordingly
                mu_D_adjust = mu_D/alpha(j);
                mu_lambda_adjust = sqrt(mu_D_adjust/mu_d);
            else
                % if Dx is fixed while Dy varies, mu_D and mu_lambda are
                % unchanged
                mu_D_adjust = mu_D;
                mu_lambda_adjust = mu_lambda;           
            end
    
            % analytical deterministic solution
            C = @(x) mu_p/mu_d * ((x<0) .* (1-cosh(x/mu_lambda_adjust)) + sinh(LS/mu_lambda_adjust) / sinh((LS+LP)/mu_lambda_adjust) * cosh((LP-x)/mu_lambda_adjust));
    
            % loop over several independent runs
            for k = 1:nruns
    
                % draw random kinetic parameters for each cell
                p = random(logndist(mu_p, CV), ncX, ncY(i));
                d = random(logndist(mu_d, CV), ncX, ncY(i));
                Dx = random(logndist(mu_D_adjust, CV), ncX, ncY(i));
                p(ncS+1:end,:) = 0; % no production in the patterning domain
        
                % inflate them to the grid
                p = repelem(p, res, res);
                d = repelem(d, res, res);
                Dx = repelem(Dx, res, res);
        
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
                    C_new = (h^2 * p + Dx .* (Cx + alpha(j) * Cy)) ./ (2 * (1 + alpha(j)) * Dx + h^2 * d);
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
                
                
                    % threshold concentration
                    C_theta = C(readouts * mu_lambda);
    
                    % determine the readout position
                    indices = find(diff(sign(C_readout - C_theta)));
                    
                    x_theta_all = [];
                    for idx = indices
                        x_theta_all = [x_theta_all, interp1(C_readout([idx idx+1]), x([idx idx+1]), C_theta)];
                    end
                    x_theta(k,y) = mean(x_theta_all);
                end
            end
    
    
            % average CV's and sigma_x
            for i = 1:numel(ncY)
                CV_lambda(j,i) = nanmean(nanstd(lambda,0,1)./nanmean(lambda,1));
                CV_lambda_SE(j,i) = nanstd(bootstrp(nboot, @(x) nanmean(nanstd(x,0,1)./nanmean(x,1)), lambda));
                CV_0(j,i) = nanmean(nanstd(C0,0,1)./nanmean(C0,1));
                CV_0_SE(j,i) = nanstd(bootstrp(nboot, @(x) nanmean(nanstd(x,0,1)./nanmean(x,1)), C0));
                sigma_x(j,i) = nanmean(nanstd(x_theta(:,:),0,1));
                sigma_x_SE(j,i) = nanstd(bootstrp(nboot, @(x) nanmean(nanstd(x,0,1)), x_theta(:,:)));
            end
        end
        toc
        
        % write results to output files
        if writeresults
            for i = 1:numel(ncY)
                T = table();
                T.alpha = alpha';
                T.CV_lambda = CV_lambda(:,i);
                T.CV_lambda_SE = CV_lambda_SE(:,i);
                T.CV_0 = CV_0(:,i);
                T.CV_0_SE = CV_0_SE(:,i);
                T.sigma_x = sigma_x(:,i);
                T.sigma_x_SE = sigma_x_SE(:,i);
                writetable(T, [prefix 'readout_' num2str(readouts) 'lambda_' num2str(ncY(i)) 'ncY.csv']);
            end
        end
    end
else
    % read existing data
    for i = 1:numel(ncY)
        T = readtable([prefix 'readout_' num2str(readouts) 'lambda_' num2str(ncY(i)) 'ncY.csv']);
        alpha = T.alpha';
        CV_lambda(:,i) = T.CV_lambda;
        CV_lambda_SE(:,i) = T.CV_lambda_SE;
        CV_0(:,i) = T.CV_0;
        CV_0_SE(:,i) = T.CV_0_SE;
        sigma_x(:,i) = T.sigma_x;
        sigma_x_SE(:,i) = T.sigma_x_SE;
    end
end


% plot results
if plotresults
    close all
    figure('Position', [0 0 1600 500]);

    % linear plot
    for i = 1:numel(ncY)
        name = ['ncY = ' num2str(ncY(i))];

        subplot(1,3,1)
        hold on
        errorbar(alpha, CV_lambda(:,i), CV_lambda_SE(:,i), 'o', 'DisplayName', name, 'LineWidth', LineWidth, 'Color', colors(i,:))
        xlabel('\alpha = D_y/D_x')
        ylabel('Average decay length variability  CV_\lambda')
        set(gca, 'LineWidth', LineWidth, 'FontSize', FontSize, 'XScale', 'log', 'YScale', 'log')
        box on
        xlim([min(alpha) max(alpha)])
    
        subplot(1,3,2)
        hold on
        errorbar(alpha, CV_0(:,i), CV_0_SE(:,i), 'o', 'DisplayName', name,  'LineWidth', LineWidth, 'Color', colors(i,:))
        xlabel('\alpha = D_y/D_x')
        ylabel('Average amplitude variability  CV_0')
        set(gca, 'LineWidth', LineWidth, 'FontSize', FontSize, 'XScale', 'log', 'YScale', 'log')
        box on
        xlim([min(alpha) max(alpha)])

        subplot(1,3,3)
        hold on
        errorbar(alpha, sigma_x(:,i)./diameter, sigma_x_SE(:,i)./diameter, 'o', 'DisplayName', name, 'LineWidth', LineWidth, 'Color', colors(i,:))
        xlabel('\alpha = D_y/D_x')
        ylabel('Positional error  \sigma_x/\delta_x (cells)')
        legend('show')
        set(gca, 'LineWidth', LineWidth, 'FontSize', FontSize, 'XScale', 'log', 'YScale', 'log')
        box on    
        xlim([min(alpha) max(alpha)])
    end
    
%     % square-root scaling for ncY=40
%     for i = 3
%         figure('Position', [0 0 1600 500]);
%         name = ['ncY = ' num2str(ncY(i))];
% 
%         subplot(1,3,1)
%         hold on
%         mdl = fitnlm(alpha(fitrange), CV_lambda(fitrange,i), @(p,x) p(1)./sqrt(x)+p(2), [1 0]);
%         CV_lambda_inf = mdl.Coefficients.Estimate(2);
%         plot(alpha, feval(mdl, alpha) - CV_lambda_inf, '-', 'HandleVisibility', 'off', 'LineWidth', LineWidth, 'Color', colors(i,:))
%         errorbar(alpha, CV_lambda(:,i) - CV_lambda_inf, CV_lambda_SE(:,i), 'o', 'DisplayName', name, 'LineWidth', LineWidth, 'Color', colors(i,:))
%         xlabel('\alpha = D_y/D_x')
%         ylabel('Average decay length variability  CV_\lambda')
%         set(gca, 'LineWidth', LineWidth, 'FontSize', FontSize, 'XScale', 'log', 'YScale', 'log')
%         box on
%         xlim([min(alpha) max(alpha)])
%     
%         subplot(1,3,2)
%         hold on
%         mdl = fitnlm(alpha(fitrange), CV_0(fitrange,i), @(p,x) p(1)./sqrt(x)+p(2), [1 0]);
%         CV_0_inf = mdl.Coefficients.Estimate(2);
%         plot(alpha, feval(mdl, alpha) - CV_0_inf, '-', 'HandleVisibility', 'off', 'LineWidth', LineWidth, 'Color', colors(i,:))
%         errorbar(alpha, CV_0(:,i) - CV_0_inf, CV_0_SE(:,i), 'o', 'DisplayName', name,  'LineWidth', LineWidth, 'Color', colors(i,:))
%         xlabel('\alpha = D_y/D_x')
%         ylabel('Average amplitude variability  CV_0')
%         set(gca, 'LineWidth', LineWidth, 'FontSize', FontSize, 'XScale', 'log', 'YScale', 'log')
%         box on
%         xlim([min(alpha) max(alpha)])
%         
%         subplot(1,3,3)
%         hold on
%         mdl = fitnlm(alpha(fitrange), sigma_x(fitrange,i)./diameter, @(p,x) p(1)./sqrt(x)+p(2), [1 0]);
%         pe_inf = mdl.Coefficients.Estimate(2);
%         plot(alpha, feval(mdl, alpha) - pe_inf, '-', 'HandleVisibility', 'off', 'LineWidth', LineWidth, 'Color', colors(i,:))
%         errorbar(alpha, sigma_x(:,i)./diameter - pe_inf, sigma_x_SE(:,i)./diameter, 'o', 'DisplayName', name, 'LineWidth', LineWidth, 'Color', colors(i,:))
%         xlabel('\alpha = D_y/D_x')
%         ylabel('Positional error  \sigma_x/\delta_x (cells)')
%         legend('show')
%         set(gca, 'LineWidth', LineWidth, 'FontSize', FontSize, 'XScale', 'log', 'YScale', 'log')
%         box on    
%         xlim([min(alpha) max(alpha)])
%     end

end