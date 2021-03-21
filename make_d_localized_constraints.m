function [Rsupport,Msupport,count, total_count] = make_d_localized_constraints(A,B,T,d,comms,ta);

Comm_speed =comms; 

Comms_Adj = abs(A)>0;
LocalityR = Comms_Adj^(d-1)>0;

if isempty(B)
    B = zeros(size(A,1),1);
end

count = zeros(1,size(A,1));
total_count = 0

for t = 1:T
    Rsupport{t} = LocalityR>0;
    Msupport{t} = (abs(B)'*Rsupport{t})>0;
    
    for j = 1:size(A,1)
        count(j) = count(j) + sum(Rsupport{t}(:,j)) + sum(Msupport{t}(:,j))
    end 
    total_count = total_count + sum(sum(Rsupport{t}))+sum(sum(Msupport{t}));
end