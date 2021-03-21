function obj_inf = inf_SLS_cost(Nx,Nu,A,B,d)

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

% % Find sparsity patterns of each columns of A in order to put back elements
% neighbor = cell(Nx,1);
% for k = 1:Nx
%     neighbor{k} = find(sparsity_A(k,:));
% end



obj_inf = 0;
% parfor i = 1:Nx
%     obj_inf = S{i}(find(internal_list{i}==i),find(internal_list{i}==i)) + obj_inf;
% end
end