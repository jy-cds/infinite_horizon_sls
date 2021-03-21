% clc; clear; close all;
Nx_list = 10:20:50
actuation = 0.5; % actuation density
alpha = 0.4;
rho = 1.6; % marginally stable
comms = 1; % commm speed between controllers
ta = 1;
% lqr_time = []
% llqr_time = []
% llqr_inf_time = []
d = 4; % d-hop sparsity
T = 8;
    
for i = 1:length(Nx_list)
    Nx = Nx_list(i)

    % Construct system matrices
    [A,B] = generate_dbl_stoch_chain(Nx,alpha,rho,actuation);
    [~,Nu] = size(B); % number of actuators
    C = [speye(Nx); sparse(Nu,Nx)];
    D = [sparse(Nx,Nu); speye(Nu)];

    llqr = @() sf_sls_d_localized_column(A,B,C,D,T,d,comms,ta,'H2')
    llqr_time(i) = timeit(llqr)
    
    lqr = @() sf_sls_d_localized(A,B,C,D,T,d,comms,ta,'H2')
    lqr_time(i) = timeit(lqr)
 
    llqr_inf = @() inf_SLS_cost(Nx,Nu,A,B,d)
    llqr_inf_time(i) = timeit(llqr_inf)
    
    
    
end

%% Plot the one figure

figure(2)
% title('Parallel Time Comparison')
set(gca,'yscale','log')
hold on

scatter(Nx_list, llqr_time,'o','LineWidth',2)
hold on
scatter(Nx_list, llqr_inf_time./Nx_list,'v','LineWidth',2)
legend("FIR SLS Controller-Localized","Proposed Infinite-horizon Controller")
xlabel("Number of Subsystems")
ylabel("Time (s)")
set(gca,'FontSize',14,'fontWeight','bold')
set(findall(gcf,'type','text'),'FontSize',14,'fontWeight','bold')
box on




%% plotting
red =[0.8500, 0.3250, 0.0980]
yellow =[0.9290, 0.6940, 0.1250]

figure(2)
subplot(1,2,1)
title('Parallel Time Comparison')
set(gca,'yscale','log')
hold on

scatter(Nx_list, llqr_time,'ro','LineWidth',2)
hold on
scatter(Nx_list, llqr_inf_time./Nx_list,'bv','LineWidth',2)
legend("FIR SLS Controller-Localized","Proposed Infinite-horizon Controller")
xlabel("Number of Subsystems")
ylabel("Time (s)")
set(gca,'FontSize',14,'fontWeight','bold')
set(findall(gcf,'type','text'),'FontSize',14,'fontWeight','bold')
box on

subplot(1,2,2)
title('Serial Time Comparison')
set(gca,'yscale','log')
hold on
scatter(Nx_list,lqr_time,'k+','LineWidth',2)
hold on
scatter(Nx_list, llqr_time .* Nx_list,'ro','LineWidth',2)
hold on
scatter(Nx_list, llqr_inf_time,'bv','LineWidth',2)
legend("FIR SLS Controller-Centralized","FIR SLS Controller-Localized","Proposed Infinite-horizon Controller")
xlabel("Number of Subsystems")
ylabel("Time (s)")
set(gca,'FontSize',14,'fontWeight','bold')
set(findall(gcf,'type','text'),'FontSize',14,'fontWeight','bold')
box on

