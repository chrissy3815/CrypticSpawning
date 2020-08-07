%%% This is a script to make some preliminary plots for Chapt 4. 
% This is with density dependence operating on the combined stock of stage
% 2 and stage 3 individuals.

%% Load necessary files, and set parameter values/ranges

addpath /Users/chrissy/Dropbox/SlopeSeaNotability_and_Matlab/

% parameters that are held constant between models:
sigma_1 = 0.2;
sigma_2 = 0.4;
sigma_3 = 0.8;
gamma = 0.33;
delta = 0.2;

a = sigma_1*gamma/(1-sigma_2*(1-delta));
b = sigma_2*delta/(1-sigma_3);

% Parameters that vary between the two models:

% "Assumed" model:
% phi_3_hat = [1e4, 5e5, 7.5e5, 1e6]; 
phi_3_hat = 1e2:1e2:1e6;
nu_3_hat = 0.05;
alpha_3_hat = 1000;%0.8;
w_hat = 3;

% "True" model:
phi_2_tilde = 0:1e2:5e5;
nu_2_tilde = 0.05;
alpha_2_tilde = 1000;%0.8;
phi_3_tilde = zeros(length(phi_3_hat), length(phi_2_tilde));
nu_3_tilde = 0.05;
alpha_3_tilde = 1000; %0.8;
w_tilde = 3;



%% Numerically find the coincident equilibria (trade-off between phi_2_tilde and phi_3_tilde)

% initialize a figure for these results:
figure

% iterate across the possible values of phi_3_hat:
for i = 1:length(phi_3_hat)
    
    % equilibrium of the assumed model: 
    nhat1_eq_expected = alpha_3_hat/(a*(1+w_hat*b))*(phi_3_hat(i)*nu_3_hat*a*b/(1-sigma_1*(1-gamma))-1);
    nhat_eq_expected = [nhat1_eq_expected, a*nhat1_eq_expected, a*b*nhat1_eq_expected];
    
    % iterate across the possible values of phi_2_tilde:
    for j = 1:length(phi_2_tilde)
        phi_2_j = phi_2_tilde(j);
        params = [sigma_1, sigma_2, sigma_3, gamma, delta, phi_3_hat(i),...
                nu_3_hat, alpha_3_hat, w_hat, nu_2_tilde, alpha_2_tilde, ...
                nu_3_tilde, alpha_3_tilde, w_tilde];

        findSSDj = @(phi_3_tilde) findSSD(params, phi_2_j, phi_3_tilde);

        phi_3_tilde(i,j) = fminbnd(findSSDj,0, 5e6);

    end
    
    plot(phi_2_tilde, phi_3_tilde(i,:), '-')
    hold on
end

title('Reproductive parameters')
xlabel('phi 2 tilde')
ylabel('phi 3 tilde')
%legend({['phi3hat =', num2str(phi_3_hat(1))], ['phi3hat =', num2str(phi_3_hat(2))], ['phi3hat =', num2str(phi_3_hat(3))]});

%% Find MSY in the assumed model:

% Values of fishing mortality to test:
F = 0.1:0.001:1;
Y = zeros(length(phi_3_hat), length(F));
n3_hat_eq = zeros(length(phi_3_hat), length(F));

% set the number of timesteps:
tmax = 500;

% initialize ntilde vectors:
n1_hat = zeros(1,tmax);
n2_hat = zeros(1,tmax);
n3_hat = zeros(1,tmax);

% set the initial condition:
n1_hat(1) = 10;
n2_hat(1) = 10;
n3_hat(1) = 10;

% maximum sustainable yield outputs:
MSY_hat = zeros(length(phi_3_hat), 1);
Bmsy_hat = zeros(length(phi_3_hat), 1);
Fmsy_hat = zeros(length(phi_3_hat), 1);

% iterate across the possible values of phi_3_hat:
for i = 1:length(phi_3_hat)

    for j=1:length(F)

        % run the model for tmax steps:
        for k = 2:tmax
            na = n2_hat(k-1)+w_hat*n3_hat(k-1);
            n1_hat(k) = phi_3_hat(i)*nu_3_hat/(1+na/alpha_3_hat)*n3_hat(k-1) + sigma_1*(1-gamma)*n1_hat(k-1);
            n2_hat(k) = sigma_1*gamma*n1_hat(k-1) + sigma_2*(1-delta)*n2_hat(k-1);
            n3_hat(k) = sigma_2*delta*n2_hat(k-1) + sigma_3*n3_hat(k-1)*(1-F(j));
        end

        n3_hat_eq(i,j) = n3_hat(tmax);
        Y(i,j) = n3_hat(tmax)*F(j);

    end
    
    % find the maximum sustainable yield info:
    MSY_hat(i) = max(Y(i,:));
    Bmsy_hat(i) = n3_hat_eq(i,Y(i,:)==max(Y(i,:)));
    Fmsy_hat(i) = F(Y(i,:)==max(Y(i,:)));
    
end

% make the plot of Yield and Biomass vs. Fishing Mortality Rate:
figure
yyaxis left
plot(F, Y(100,:), '-')
hold on
plot(F, Y(500,:), '--')
plot(F, Y(900,:), ':')
xlabel("Fishing Mortality Rate")
ylabel("Yield")
yyaxis right
plot(F, n3_hat_eq(100,:), '-')
plot(F, n3_hat_eq(500,:), '--')
plot(F, n3_hat_eq(900,:), ':')
ylabel("Stage 3 Biomass")
title('Harvesting in the Assumed model')
legend({['phi3hat =', num2str(phi_3_hat(100))], ['phi3hat =', num2str(phi_3_hat(500))],...
    ['phi3hat =', num2str(phi_3_hat(900))]});

% make a plot of MSY vs. phi_3_hat:
figure
yyaxis left
plot(phi_3_hat, Bmsy_hat, '--')
hold on
plot(phi_3_hat, MSY_hat, '-')
ylabel('Yield or Biomass')
yyaxis right
plot(phi_3_hat, Fmsy_hat, '-')
ylabel('Fmsy')
ylim([0 0.5])
xlabel('phi 3 hat')
title('MSY vs. phi 3 hat')
legend({'Bmsy', 'MSY', 'Fmsy'})

%% Calculate MSY, etc. in the true model


% Values of fishing mortality to test:
F = 0.1:0.02:1;
Y = zeros(length(phi_3_hat), length(phi_2_tilde), length(F));
n3_tilde_eq = zeros(length(phi_3_hat), length(phi_2_tilde), length(F));

% set the number of timesteps:
tmax = 500;

% initialize ntilde vectors:
n1_tilde = zeros(1,tmax);
n2_tilde = zeros(1,tmax);
n3_tilde = zeros(1,tmax);

% set the initial condition:
n1_tilde(1) = 10;
n2_tilde(1) = 10;
n3_tilde(1) = 10;

% maximum sustainable yield outputs:
MSY_tilde = zeros(length(phi_3_hat), length(phi_2_tilde));
Bmsy_tilde = zeros(length(phi_3_hat), length(phi_2_tilde));
Fmsy_tilde = zeros(length(phi_3_hat), length(phi_2_tilde));

for p=1:length(phi_3_hat)
    
    for k=1:length(phi_2_tilde)

        for j=1:length(F)

            % run the model for tmax steps:
            for i = 2:tmax
                na = n2_tilde(i-1)+w_tilde*n3_tilde(i-1);
                n1_tilde(i) = phi_2_tilde(k)*nu_2_tilde/(1+na/alpha_2_tilde)*n2_tilde(i-1) + phi_3_tilde(p,k)*nu_3_tilde/(1+na/alpha_3_tilde)*n3_tilde(i-1) + sigma_1*(1-gamma)*n1_tilde(i-1);
                n2_tilde(i) = sigma_1*gamma*n1_tilde(i-1) + sigma_2*(1-delta)*n2_tilde(i-1);
                n3_tilde(i) = sigma_2*delta*n2_tilde(i-1) + sigma_3*n3_tilde(i-1)*(1-F(j));
            end

            n3_tilde_eq(p,k,j) = n3_tilde(tmax);
            Y(p,k,j) = n3_tilde(tmax)*F(j);

        end
        
        % find the maximum sustainable yield info:
        MSY_tilde(p,k) = max(Y(p,k,:));
        I = find(Y(p,k,:)==max(Y(p,k,:)));
        Bmsy_tilde(p,k) = n3_tilde_eq(p,k,I);
        Fmsy_tilde(p,k) = F(I);
        
    end
    
end

% plot some things?
% make the plot of Yield and Biomass vs. Fishing Mortality Rate:
figure
yyaxis left
plot(F, squeeze(Y(5000,11,:)), '-')
hold on
plot(F, squeeze(Y(5000,21,:)), '--')
plot(F, squeeze(Y(5000, 31,:)), ':')
xlabel("Fishing Mortality Rate")
ylabel("Yield")
yyaxis right
plot(F, squeeze(n3_tilde_eq(5000,11,:)), '-')
plot(F, squeeze(n3_tilde_eq(5000,21,:)), '--')
plot(F, squeeze(n3_tilde_eq(5000,31,:)), ':')
ylabel("Stage 3 Biomass")
title(['Harvesting in the True model, phi3hat =', num2str(phi_3_hat(5000))])
legend({['phi2tilde =', num2str(phi_2_tilde(11))], ['phi2tilde =', num2str(phi_2_tilde(21))],...
    ['phi2tilde =', num2str(phi_2_tilde(31))]});

% make a plot of MSY vs. phi_3_hat:
figure
yyaxis left
plot(phi_2_tilde, Bmsy_tilde(5000,:), '--')
hold on
plot(phi_2_tilde, MSY_tilde(5000,:), '-')
ylabel('Yield or Biomass')
yyaxis right
plot(phi_2_tilde, Fmsy_tilde(5000,:), '-')
ylabel('Fmsy')
xlabel('phi 2 tilde')
title('MSY vs. phi 2 tilde')
legend({'Bmsy', 'MSY', 'Fmsy'})

%% Plots of how yield measures change between the true and assumed models:

% delta Fmsy plot:
figure
hold on
for i=1:length(phi_3_hat)
    toplot = Fmsy_tilde(i,:)-Fmsy_hat(i);
    plot(phi_2_tilde, toplot, '-')
    hold on
end
title('Fmsy tilde - Fmsy hat')
xlabel('phi 2 tilde')
ylabel('delta Fmsy')
legend({['phi3hat =', num2str(phi_3_hat(1))], ['phi3hat =', num2str(phi_3_hat(2))],...
    ['phi3hat =', num2str(phi_3_hat(3))]});

% delta MSY plot:
figure
hold on
for i=1:length(phi_3_hat)
    toplot = MSY_tilde(i,:)-MSY_hat(i);
    plot(phi_2_tilde, toplot, '-')
    hold on
end
title('MSY tilde - MSY hat')
xlabel('phi 2 tilde')
ylabel('delta MSY')
legend({['phi3hat =', num2str(phi_3_hat(1))], ['phi3hat =', num2str(phi_3_hat(2))],...
    ['phi3hat =', num2str(phi_3_hat(3))]});

% delta Bmsy plot:
figure
hold on
for i=1:length(phi_3_hat)
    toplot = Bmsy_tilde(i,:)-Bmsy_hat(i);
    plot(phi_2_tilde, toplot, '-')
    hold on
end
title('Bmsy tilde - Bmsy hat')
xlabel('phi 2 tilde')
ylabel('delta Bmsy')
legend({['phi3hat =', num2str(phi_3_hat(1))], ['phi3hat =', num2str(phi_3_hat(2))],...
    ['phi3hat =', num2str(phi_3_hat(3))]});

% The effect of harvesting the true population at the Fmsy of the assumed model:
figure
hold on
for i=1:length(phi_3_hat)
    I = find(F==Fmsy_hat(i));
    toplot = squeeze(Y(i,:,I))-MSY_hat(i);
    plot(phi_2_tilde, toplot, '-')
    hold on
end
title('Harvesting the true model with Fmsy hat (Yield)')
xlabel('phi 2 tilde')
ylabel('Ytilde(Fmsy hat) - MSY hat')
legend({['phi3hat =', num2str(phi_3_hat(1))], ['phi3hat =', num2str(phi_3_hat(2))],...
    ['phi3hat =', num2str(phi_3_hat(3))]});

% The effect of harvesting the true population at the Fmsy of the assumed model:
figure
hold on
for i=1:length(phi_3_hat)
    I = find(F==Fmsy_hat(i));
    toplot = squeeze(n3_tilde_eq(i,:,I))-Bmsy_hat(i);
    plot(phi_2_tilde, toplot, '-')
    hold on
end
title('Harvesting the true model with Fmsy hat (Biomass)')
xlabel('phi 2 tilde')
ylabel('Btilde(Fmsy hat) - Bmsy hat')
legend({['phi3hat =', num2str(phi_3_hat(1))], ['phi3hat =', num2str(phi_3_hat(2))],...
    ['phi3hat =', num2str(phi_3_hat(3))]});

%% Plot of n3star_hat against phi3_hat

nhat1_eq_nofishing = alpha_3_hat/(a*(1+w_hat*b))*(phi_3_hat*nu_3_hat*a*b/(1-sigma_1*(1-gamma))-1);
n3_hat_eq_nofishing = a*b*nhat1_eq_nofishing;


figure
plot(phi_3_hat, n3_hat_eq_nofishing)
hold on
plot(phi_3_hat, Bmsy_hat)
legend({'n3star, no fishing', 'n3star, with fishing'})
ylim([0 4.5e5])

