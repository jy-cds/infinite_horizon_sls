% clc; clear; close all;
Nx_list = 10:20:500
actuation = 0.5; % actuation density
alpha = 0.4;
rho = 1.6; % marginally stable
comms = 1; % commm speed between controllers
ta = 1;
FIR_obj = []
obj_inf = []
d = 4; % d-hop sparsity
    T = 8;
    
for i = 1:length(Nx_list)
    Nx = Nx_list(i)

    
    % Construct system matrices
    [A,B] = generate_dbl_stoch_chain(Nx,alpha,rho,actuation);
    [~,Nu] = size(B); % number of actuators
    
    % normal SLS
    C = [speye(Nx); sparse(Nu,Nx)];
    D = [sparse(Nx,Nu); speye(Nu)];

    [R_T,M_T,FIR_obj(i)]  = sf_sls_d_localized(A,B,C,D,T,d,comms,ta,'H2')
    obj_inf(i) = inf_SLS_cost(Nx,Nu,A,B,d)
    
end




%% plotting
figure(2)
hold on
scatter(Nx_list,FIR_obj,'+','LineWidth',2)
hold on
scatter(Nx_list, obj_inf,'o','LineWidth',2)
hold on
scatter(Nx_list(find(isnan(FIR_obj))),ones(1,length(find(isnan(FIR_obj))))*1199, 'v','LineWidth',2)
legend("FIR SLS Controller","Proposed Infinite-horizon Controller")
xlabel("Number of Subsystems")
ylabel("H_2 Cost")
set(gca,'FontSize',14,'fontWeight','bold')
set(findall(gcf,'type','text'),'FontSize',14,'fontWeight','bold')
box on

