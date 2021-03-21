%% Global State Synthesis

clc; clear; close all;

Nx = 10; % states
loc = floor(Nx/2); % specify where disturbance hits

actuation = 0.5; % actuation density
alpha = 0.4;
rho = 1.25; % marginally stable

% Construct system matrices
[A,B] = generate_dbl_stoch_chain(Nx,alpha,rho,actuation);
[~,Nu] = size(B); % number of actuators

% A = [1     1     1     0     0;
%      1     1     0     1     0;
%      0     1     1     0     1;
%      1     0     1     1     0;
%      0     0     1     0     1];
 

% B(1,5) = 1;
% B(2,3) = 1;
% B(4,3) = 1;
% B(4,5) = 1;


% locality
d = 2; % d-hop sparsity
comms = 1; % commm speed between controllers


%% Local K Computation
D = eye(Nx);

nx = {}; % number of internal states for each column
nu = {}; % number of control for each column

sparsity_A = abs(A^(d-1))>0;
sparsity_A_extend = abs(A^(d))>0;
boundary_pattern = abs(sparsity_A_extend - sparsity_A)>0;
sparsity_u = (abs(B')>0)*abs(A^(d))>0;

% sparsity_u = ones(5) %%%%%%%%%%%%%%%%%% DELETE!! %%%%%%%%%%%%%%

e = cell(Nx,1); % modified e_j vector
K = {};
S = {};
internal_list = {};
boundary_list = {};
control_list = {};
phi = cell(Nx,1);
psi = cell(Nx,1);
j_hat = cell(Nx,1);

for i = 1:Nx
    A_nn = [];
    A_bn = [];
    A_bb = [];
    A_nb = [];
    
    B_nn = [];
    B_bb = [];
    D_12 = [];
    
    nx{i} = nnz(sparsity_A(:,i));
    nu{i} = nnz(sparsity_u(:,i));
    
    internal = find(sparsity_A(:,i)); % find which states are internal
    boundary = find(boundary_pattern(:,i)); % find which states are on the boundary
    control = find(sparsity_u(:,i));
    
    e_temp = zeros(nx{i},1);
    e_temp(find(internal==i)) = 1;
    j_hat{i} = find(internal==i);
    e{i} = e_temp;
    
    internal_list{i} = internal;
    boundary_list{i} = boundary;
    control_list{i} = control;
    
    for k = 1:size(internal,1)
        A_nn = [A_nn; A(internal(k),internal)];
        A_nb = [A_nb; A(internal(k),boundary)];
        B_nn = [B_nn;B(internal(k),control)];
    end
    
    for k = 1:size(boundary,1)
        A_bb = [A_bb; A(boundary(k),boundary)];
        A_bn = [A_bn; A(boundary(k),internal)];
        B_bb = [B_bb;B(boundary(k),control)];
    end
    for k = 1:size(internal,1)
        D_12 = [D_12;D(internal(k),control)];
    end
    
    lqr_A = A_nn - B_nn*pinv(full(B_bb))*A_bn;
    lqr_B = B_nn*(eye(size(B_nn,2))-pinv(full(B_bb))*B_bb);
    lqr_Q = (eye(nx{i}) - D_12*pinv(full(B_bb))*A_bn)'*(eye(nx{i}) - D_12*pinv(full(B_bb))*A_bn);
    lqr_P = eye(nu{i});
    %lqr_P = (D_12*(eye(size(B_nn,2))-pinv(full(B_bb))*B_bb))'*(D_12*(eye(size(B_nn,2))-pinv(full(B_bb))*B_bb));
    
    
%     K{i} = dlqr(lqr_A, lqr_B, lqr_Q, lqr_P);
    [K{i},S{i},~] = dlqr(lqr_A, lqr_B, lqr_Q, lqr_P);
    psi{i} = -pinv(full(B_bb))*A_bn -(eye(size(B_nn,2))-pinv(full(B_bb))*B_bb)*K{i};
    phi{i} = lqr_A - lqr_B*K{i};
    
end

%% Controller Implementation
loc = 5;
Tmax = 20;
delta_A = 0*eye(Nx);


% Get the neighbors for each node
neighbor = cell(Nx,1);
for k = 1:Nx
    neighbor{k} = find(sparsity_A(k,:));
end

% global plant variables
w = zeros(Nx,Tmax);
w(loc,1) = 5;
% w(loc+5,1) = 0.5;



x = zeros(Nx,Tmax);
u = zeros(Nu,Tmax);

% local variables
xi_list = {};
w_hat_list = {};
u_local_list ={};

for i = 1:Nx
    xi_list{i} = zeros(nx{i},Tmax);
    w_hat_list{i} = zeros(1,Tmax);
    u_local_list{i} = zeros(nu{i},Tmax);
end

%Time domain simulation

MATx = [];
MATu = [];
for t = 1:Tmax
    t;
    for i = 1:Nx
        i;
        xi = xi_list{i};
        u_local = u_local_list{i};
        w_hat = w_hat_list{i};
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        x_ref = 0;
        neighbors = neighbor{i};
        for j = 1:length(neighbors) % summing all the \xi^j_i
            j; % jth neighbor
            neighbor_xi = xi_list{neighbors(j)};
            x_ref = x_ref + neighbor_xi(find(internal_list{neighbors(j)}==i),t);
        end
        w_hat(t) = x(i,t) - x_ref;
        w_hat_list{i}=w_hat;
        
        
        w_hat = w_hat_list{i};
        u_local(:,t) = psi{i}*(xi(:,t) + e{i}*w_hat(t));
        xi(:,t+1) = phi{i}*(xi(:,t)+e{i}*w_hat(t));
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        xi_list{i} = xi ;
        u_local_list{i} =  u_local ;
        
        
        % global control action computation
        temp = zeros(Nu,1);
        temp(control_list{i}) = u_local(:,t);
        temp
        u(:,t) = u(:,t) + temp;
        
    end
    
    x(:,t+1) = (A+delta_A)*x(:,t) + w(:,t) +  B*u(:,t) ;
    
    MATx = [MATx,log10(abs(x(:,t)))];
    MATu = [MATu,log10(abs(B*u(:,t)))];
    
    
    
    
    
end


%% Plotting state and control trajectory
for i = [5]
    figure(i)
    
suptitle('x_5(t) and u_5(t)')
    subplot(2,1,1);
    plot(1:Tmax,x(i,1:end-1),'bo')
    %plot(Tstart:Tmax,w_d(loc,Tstart:end),'mo');
   % legend(['x' num2str(i)])
    xlabel('t')
    ylabel(['x' num2str(i) '(t)'])
    
    subplot(2,1,2);
    plot(1:Tmax,u(i,1:end),'bo')
    hold on
   % legend('control effort')
    ylabel(['u {' num2str(i) '}(t)'])
    xlabel('t')
    currentFigure = gcf;
    title(currentFigure.Children(end), 'Closed Loop');
    hold on
end



%% Heatmap for Infinite Horizon
figure
suptitle('Infinite Horizon')
subplot(1,2,1);
handle = imagesc(((MATx)));
colormap jet
caxis([-3 1])
title('log10(|x|)')
xlabel('Time')
ylabel('Space')
set(gca,'FontSize',16,'fontWeight','bold')
set(findall(gcf,'type','text'),'FontSize',16,'fontWeight','bold')

subplot(1,2,2);
handle = imagesc(((MATu)));
colorbar
colormap jet
caxis([-3 1])
title('log10(|u|)')
xlabel('Time')
set(gca,'FontSize',16,'fontWeight','bold')
set(findall(gcf,'type','text'),'FontSize',16,'fontWeight','bold')


%% Normal deadbeat SLS

C = [speye(Nx); sparse(Nu,Nx)];
D = [sparse(Nx,Nu); speye(Nu)];
ta = 1;

T = 8;
d = 4;
comms = 1;
[R_T,M_T,obj_T]  = sf_sls_d_localized(A,B,C,D,T,d,comms,ta,'H2');

[x_T,u_T] = make_heat_map(A+delta_A,B,T,Nx,Nu,R_T,M_T,loc,Tmax+T,'',0)

figure(5)
hold on
subplot(2,1,1);
hold on
plot(1:Tmax,x_T(5,T+1:end),'r+')
%plot(Tstart:Tmax,w_d(loc,Tstart:end),'mo');
legend('x_{infinite}','x_{finite}')
xlabel('t')
ylabel(['x' num2str(i) '(t)'])
set(gca,'FontSize',13,'fontWeight','bold')
set(findall(gcf,'type','text'),'FontSize',13,'fontWeight','bold')


subplot(2,1,2);
plot(1:Tmax,u_T(5,T+1:end),'r+')
hold on
legend('u_{infinite}','u_{finite}')
ylabel(['u{' num2str(i) '}(t)'])
xlabel('t')
currentFigure = gcf;
hold on
set(gca,'FontSize',13,'fontWeight','bold')
set(findall(gcf,'type','text'),'FontSize',13,'fontWeight','bold')



