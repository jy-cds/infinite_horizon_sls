function [R_concat, M_concat, clnorm] = sf_sls_d_localized(A,B,C,D,T,d,comms,ta,obj);

% System level synthesis with d-localizable constraints
% Built around cvx
% Returns system response (R,M) as defined in Thm 1 of https://arxiv.org/pdf/1610.04815.pdf
Nx = size(A,1);
Nu = size(B,2);

[Rsupport,Msupport,count,total_count] = make_d_localized_constraints(A,B,T,d,comms,ta);

clnorm = 0;

R_concat = {};
M_concat = {};
identity = eye(Nx);

for j = 1:Nx
    R = {};
    M = {};
    
    cvx_begin
    %cvx_solver sedumi
    variable X(count(j))
    expression Rs(Nx,1,T)
    expression Ms(Nu,1,T)

    %populate decision variables
    %locality constraints automatically enforced by limiting support of R and M
    spot = 0;
    for t = 1:T
        R{t} = Rs(:,:,t);
        supp = find(Rsupport{t}(:,j));
        num = sum(sum(Rsupport{t}(:,j)));
        R{t}(supp) = X(spot+1:spot+num);
        spot = spot + num;

        M{t} = Ms(:,:,t);
        supp = find(Msupport{t}(:,j));
        num = sum(sum(Msupport{t}(:,j)));
        M{t}(supp) = X(spot+1:spot+num);
        spot = spot + num;
    end
    objective = 0;

    %set up objective function
    switch obj
        case 'H2'
            objective = compute_H2(R,M,C,D,T);
        case 'Hinf'
            % todo
        case 'L1'
            % todo
        case 'L1T'
            % todo
        otherwise
            % todo: throw a warning that only fidning a feasible solution
    end
    cvx_precision low
    minimize(objective)
    subject to
        %achievability constraints
        R{1} == identity(:,j);

        for t= 1:T-1
            R{t+1} == A*R{t} + B*M{t};
        end
        R{T} == zeros(Nx,1);
    cvx_end
    
    clnorm = clnorm + objective;
    R_concat{j} = R;
    M_concat{j} = M;
    
end 
