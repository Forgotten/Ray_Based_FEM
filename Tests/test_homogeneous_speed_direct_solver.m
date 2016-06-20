%% One point source inside domain (homogeneous medium):  iterative idea
% script that computes the low and high frequency systems using the 
% method of polarized traces

xs = -0.3;   ys = -0.3;             % point source location
speed = @(x) ones(size(x,1),1);    % medium speed

plt = 1;                           % plot solution or not
fquadorder = 6;                    % numerical quadrature order
Nray = 1;
Rest = 1;
pde = [];
pct = 1/5;
data = 'num';
opt = 1;

high_omega = 100*pi;
low_omega = sqrt(high_omega);
NPW = 10;

h = 1/(20*round(high_omega*NPW/(2*pi*20)));
% high frequency h
hTilde = 1/(20*max(round(low_omega*NPW/(2*pi*20)),1));
% low frequency h (tildeh)

wavelength = 2*pi/high_omega;  
wpml = wavelength;               % width of PML
sigmaMax = 25/wpml;                % Maximun absorbtion
r = 8*wpml;

lg_a = 1;
md_a = 0.65;

sm_a = 1/2;

fprintf(['\n' '-'*ones(1,80) '\n']);
fprintf('omega/(2*pi) = %d,   1/h = %d   1/ch = %d,  NPW = %d \n',high_omega/(2*pi), 1/h, 1/hTilde, NPW);

sigma = 1/100;
k = high_omega;
source = @(p) -4*sqrt(k)*1i*1/(2*pi*sigma^2)*exp( -( (p(:,1)-xs).^2 + (p(:,2)-ys).^2  )/(2*sigma^2) );


%% Step 1: Solve the Hemholtz equation with the same source but with a relative low frequency sqrt(\omega) by Standard FEM, mesh size \omega*h = constant
fprintf(['\n' '-'*ones(1,80) '\n']);
fprintf('Step1: S-FEM, low frequency\n');

a = lg_a;
[nodeLow,elemLow] = squaremesh([-a,a,-a,a],h);
omega = low_omega;

M0 = model;
M0.init(nodeLow,elemLow,low_omega,wpml,h/NPW, speed, fquadorder)

% factorizing the matrix within the model class
M0.LU_factorization()


fquadorder = 4;  
% solving the PDE using the global model 
f = assemble_RHS(nodeLow,elemLow,source,fquadorder);
u_std = M0.solve(f);
    
% plotting the global solution 
if plt 
    M0.showresult(real(u_std));
end

% getting only the interior degrees of freedom
fInt = f(M0.freeNode);

%[u_std] = Standard_FEM_PML_PointSource(lnode,lelem,omega,wpml,sigmaMax,xs,ys,speed,fquadorder,plt);


%% Step 2: Use NMLA to find ray directions d_c with low frequency sqrt(\omega)
fprintf(['\n' '-'*ones(1,80) '\n']);
fprintf('\nStep2: NMLA, low frequency\n');

[ux,uy] = num_derivative(u_std,h,2);

a = md_a;
[mnode,melem] = squaremesh([-a,a,-a,a],h);
[cnode,celem] = squaremesh([-a,a,-a,a],hTilde);
cN = size(cnode,1);
cnumray = zeros(cN,Nray);
cray = ex_ray(cnode,xs,ys);
cr1 = zeros(cN,1);
cd1 = cr1;

fprintf('NMLA time: \n');
tic;
for i = 1:cN
    x0 = cnode(i,1);  y0 = cnode(i,2);
    d0 = sqrt((x0-xs)^2 + (y0-ys)^2);
    if d0 <= r
        cnumray(i,:) = cray(i,:);
    else
        Rest = d0;
        if d0 <= 2*r
            Rest = 2*d0;
        end
        c0 = speed(cnode(i,:));
        %[cnumray(i,:),~,cr1(i)] = NMLA_2D_2nd(x0,y0,c0,omega,Rest,nodeLow,elemLow,u_std,ux,uy,pde,pct,Nray,data,opt,plt);
        [cnumray(i,:)] = NMLA_2D_2nd(x0,y0,c0,omega,Rest,nodeLow,elemLow,u_std,ux,uy,pde,pct,Nray,data,opt,plt);
        cd1(i) = cr1(i) - d0;
%         if r1 > d0
%             cr1(i) = 1;
%         end
    end
end
toc;

clear lnode lelem;

cdiffang = angle_error(cnumray,cray);
norm(cdiffang,2)/norm(cray)
norm(cdiffang,inf)

cnumray = exp(1i*cnumray);
numray1 = interpolation(cnode,celem,mnode,cnumray);

ray = ex_ray(mnode,xs,ys);
ray = exp(1i*ray);
md = sqrt((mnode(:,1)-xs).^2 + (mnode(:,2)-ys).^2);
ray = ray.*(1 - (md<eps));

numray1 = numray1.*(md>r) + ray.*(md<=r);
diffray1 = numray1 - ray;
% figure(1);
% FJ_showresult(mnode,melem,real(diffray1));
figure(1);
FJ_showresult(cnode,celem,cd1);

%% Step 3: Solve the original Helmholtz equation by Ray-based FEM with ray directions d_c
fprintf(['\n' '-'*ones(1,80) '\n']);
fprintf('\nStep3: Ray-FEM, high frequency \n');

omega = high_omega;

MHigh = model;
MHigh.initRay(mnode,melem,high_omega,wpml,h/NPW, speed, numray1, fquadorder);

% % performing LU factorization (only if we want ot compare against the true
% % solution 
MHigh.LU_factorization()

% solving the PDE using the global model 
fHigh = assemble_RHS_with_ray(mnode,melem, high_omega,source, speed,numray1,fquadorder);
u1 = MHigh.solve(fHigh);


N = size(mnode,1);       % number of grid points
Nray = size(numray1,2);     % number of rays crossing at each grid node
Ndof = N*Nray; 

if plt
    
    grad = numray1(:);
    grad = [real(grad),imag(grad)];
    repnode = repmat(mnode,Nray,1);
    temp = grad(:,1).*repnode(:,1) + grad(:,2).*repnode(:,2);

    c = speed(mnode);    % medium speed
    k = omega./c;       % wavenumber
    kk = repmat(k,1,Nray);

    u = u1.*exp(1i*kk(:).*temp);
    u = reshape(u,N,Nray);
    u = sum(u,2);

% showing the real part of the solution
figure(2); 
MHigh.showresult(real(u));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% until here is should be OK
% 
% 
% %% Step 4: NMLA to find original ray directions d_o with wavenumber k
% fprintf(['\n' '-'*ones(1,80) '\n']);
% fprintf('\nStep4: NMLA, high frequency \n');
% 
% a = sm_a;
% [node,elem] = squaremesh([-a,a,-a,a],h);
% [cnode,celem] = squaremesh([-a,a,-a,a],hTilde);
% cN = size(cnode,1);
% cnumray = zeros(cN,Nray);
% cray = ex_ray(cnode,xs,ys,0);
% cr2 = zeros(cN,1);
% cd2 = cr2;
% 
% [ux,uy] = num_derivative(u1,h,2);
% 
% fprintf('NMLA time: \n');
% tic;
% for i = 1:cN
%     x0 = cnode(i,1);  y0 = cnode(i,2);
%     d0 = sqrt((x0-xs)^2 + (y0-ys)^2);
%     if d0 <= r
%         cnumray(i,:) = cray(i,:);
%     else
%         Rest = d0;
%         if d0 <= 2*r
%             Rest = 2*d0;
%         end
%         c0 = speed(cnode(i,:));
%         [cnumray(i,:)] = NMLA_2D_2nd(x0,y0,c0,omega,Rest,mnode,melem,u1,ux,uy,pde,pct,Nray,data,opt,plt);
%         cd2(i) = cr2(i) - d0;
% %         if cr > d0
% %             cr2(i) = 1;
% %         end
%     end
% end
% toc;
% 
% % clear mnode melem;
% 
% cdiffang = angle_error(cnumray,cray);
% norm(cdiffang,2)/norm(cray)
% norm(cdiffang,inf)
% 
% cnumray = exp(1i*cnumray);
% numray2 = interpolation(cnode,celem,node,cnumray);
% 
% ray = ex_ray(node,xs,ys,0);
% ray = exp(1i*ray);
% d = sqrt((node(:,1)-xs).^2 + (node(:,2)-ys).^2);
% ray = ray.*(1 - (d < eps) );
% 
% numray2 = numray2.*(d>r) + ray.*(d<=r);
% diffray2 = numray2 - ray;
% % figure(2);
% % FJ_showresult(node,elem,real(diffray2));
% 
% numray_dir = [real(numray2), imag(numray2)];
% numray_dir = atan2(numray_dir(:,2), numray_dir(:,1));
% numray_dir = numray_dir + 2*pi*(numray_dir<0);
% 
% diffang = angle_error(numray_dir,ex_ray(node,xs,ys,0));
% diffang = diffang.*(d>r);
% % figure(3);
% % FJ_showresult(node,elem,real(diffang));
% % title('NMLA angle error');
% figure(3);
% FJ_showresult(cnode,celem,cr2);
% 
% %% Step 5: Solve the original Helmholtz equation by Ray-based FEM with ray directions d_o
% fprintf(['\n' '-'*ones(1,80) '\n']);
% fprintf('\nStep5: Ray-FEM, high frequency \n');
% 
% omega = high_omega;
% [u2,~,v2] = Ray_FEM_PML_1_PointSource(node,elem,omega,wpml,sigmaMax,xs,ys,speed,numray2,fquadorder,plt);
% 
% 
% %% 
% % [X,Y] = meshgrid(-a:h:a,-a:h:a);
% % [m,n] = size(X);
% % uh = reshape(u2,m,n);
% % 
% % figure(4)
% % contour(X,Y,real(uh));
% % title('Level set of Ray-FEM solution u_{ray}');
% 
% 
% %% map to polar
% figure(5);
% % r1 = 10.3*wpml;
% % r2 = 12.5*wpml;
% r1 = 0.834;
% r2 = 0.834 + 2.1*wavelength;
% theta1 = pi/4 - pi/16;
% theta2 = pi/4 + pi/16;
% subplot(2,1,1);
% mapto_polar(node,elem,omega,speed,v2,numray2,xs,ys,r1,r2,theta1,theta2);
% 
% subplot(2,1,2);
% mapto_polar(node,elem,omega,speed,v2,numray2,xs,ys,r1,r2,theta1,theta2);
% axis equal
% 
% 
% 
% %% 
% % figure(6);
% % % r1 = 0.81;
% % % r2 = 0.81 + 2.1*wavelength;
% % % omega = 80*pi;
% % [X,Y] = meshgrid(-a:1/1000:a,-a:1/1000:a);
% % [m,n] = size(X);
% % xynode = [X(:),Y(:)];
% % uu = ray_solution(node,elem,omega,speed,v2,numray2,xynode);
% % uu = reshape(uu,m,n);
% % mesh(X,Y,real(uu));
% % az = 0;
% % el = 90;
% % view(az, el);
% % axis equal; axis tight;
% % hold on;
% % 
% % theta = theta1:1/10000:theta2;
% % rr = r1: 1/5000:r2;
% % x1 = r1*cos(theta) + xs;   y1 = r1*sin(theta) + ys;
% % x2 = r2*cos(theta) + xs;   y2 = r2*sin(theta) + ys;
% % x3 = rr*cos(theta1) + xs;  y3 = rr*sin(theta1) + ys;
% % x4 = rr*cos(theta2) + xs;  y4 = rr*sin(theta2) + ys;
% % p = plot(x1,y1,'r-');
% % p.LineWidth = 3;  hold on;
% % p = plot(x2,y2,'r-');
% % p.LineWidth = 3;  hold on;
% % p = plot(x3,y3,'r-');
% % p.LineWidth = 3;  hold on;
% % p = plot(x4,y4,'r-');
% % p.LineWidth = 3;  hold on;
% % xlabel('x');
% % ylabel('y');
