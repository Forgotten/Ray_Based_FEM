function [ux,uy] = num_derivative(u,h,order)
%% numerical central difference schemes to compute the derivative of u
% u is a N x 1 vector
N = length(u);
n = round(sqrt(N));
u = reshape(u,n,n);

ux = u; uy = u;

if order == 2
    ux(:,2:n-1) = (u(:,3:n) - u(:,1:n-2))/(2*h);
    ux(:,n) = 2*ux(:,n-1) - ux(:,n-2);
    ux(:,1) = 2*ux(:,2) - ux(:,3);    
    
    uy(2:n-1,:) = (u(3:n,:) - u(1:n-2,:))/(2*h);
    uy(n,:) = 2*uy(n-1,:) - uy(n-2,:);
    uy(1,:) = 2*uy(2,:) - uy(3,:); 
end

if order == 4   % n >= 8
    ux(:,3:n-2) = ( -u(:,5:n) + 8*u(:,4:n-1) - 8*u(:,2:n-3) + u(:,1:n-4) )/(12*h);
    ux(:,n-1) = 4*ux(:,n-3) - ux(:,n-2) - ux(:,n-4) - ux(:,n-5);
    ux(:,n) = 4*ux(:,n-2) - ux(:,n-1) - ux(:,n-3) - ux(:,n-4);
    ux(:,2) = 4*ux(:,4) - ux(:,3) - ux(:,5) - ux(:,6);
    ux(:,1) = 4*ux(:,3) - ux(:,2) - ux(:,4) - ux(:,5);
    
    uy(3:n-2,:) = ( -u(5:n,:) + 8*u(4:n-1,:) - 8*u(2:n-3,:) + u(1:n-4,:) )/(12*h);
    uy(n-1,:) = 4*uy(n-3,:) - uy(n-2,:) - uy(n-4,:) - uy(n-5,:);
    uy(n,:) = 4*uy(n-2,:) - uy(n-1,:) - uy(n-3,:) - uy(n-4,:);
    uy(2,:) = 4*uy(4,:) - uy(3,:) - uy(5,:) - uy(6,:);
    uy(1,:) = 4*uy(3,:) - uy(2,:) - uy(4,:) - uy(5,:);
end

ux = ux(:);  uy = uy(:);