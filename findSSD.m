function SSD = findSSD(params, phi_2_tilde, phi_3_tilde)

if length(params)~=14
    error('There should be 14 parameters in the first input argument')
end

sigma_1 = params(1);
sigma_2 = params(2);
sigma_3 = params(3);
gamma = params(4);
delta = params(5);

phi_3_hat = params(6);
nu_3_hat = params(7);
alpha_3_hat = params(8);
w_hat = params(9);

nu_2_tilde = params(10);
alpha_2_tilde = params(11);
nu_3_tilde = params(12);
alpha_3_tilde = params(13);
w_tilde = params(14);


% equilibrium of the assumed model: 
a = sigma_1*gamma/(1-sigma_2*(1-delta));
b = sigma_2*delta/(1-sigma_3);
nhat1_eq_expected = alpha_3_hat/(a*(1+w_hat*b))*(phi_3_hat*nu_3_hat*a*b/(1-sigma_1*(1-gamma))-1);
nhat_eq_expected = [nhat1_eq_expected, a*nhat1_eq_expected, a*b*nhat1_eq_expected];

        % calculate the equilibrium of the true model
        tmax = 200;
    
        % initialize ntilde vectors:
        n1_tilde = zeros(1,tmax);
        n2_tilde = zeros(1,tmax);
        n3_tilde = zeros(1,tmax);
        
        % set the initial condition:
        n1_tilde(1) = 10;
        n2_tilde(1) = 10;
        n3_tilde(1) = 10;
        
        % run the model for tmax steps:
        for i = 2:tmax
            na = n2_tilde(i-1)+w_tilde*n3_tilde(i-1);
            n1_tilde(i) = phi_2_tilde*nu_2_tilde/(1+na/alpha_2_tilde)*n2_tilde(i-1) + phi_3_tilde*nu_3_tilde/(1+na/alpha_3_tilde)*n3_tilde(i-1) + sigma_1*(1-gamma)*n1_tilde(i-1);
            n2_tilde(i) = sigma_1*gamma*n1_tilde(i-1) + sigma_2*(1-delta)*n2_tilde(i-1);
            n3_tilde(i) = sigma_2*delta*n2_tilde(i-1) + sigma_3*n3_tilde(i-1);
        end
        
        SSD = (nhat_eq_expected(1)-n1_tilde(end))^2 + (nhat_eq_expected(2)-n2_tilde(end))^2 + (nhat_eq_expected(3)-n3_tilde(end))^2;