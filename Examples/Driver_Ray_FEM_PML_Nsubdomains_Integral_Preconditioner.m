%% Driver for the ray based finite element with absorbing boundary conditions 
% in this case we use the correct ray information 
% implemented via PML using N different subdomains
% the system is solved using a variation of the method of polarized traces
% We test the preconditioner (without encapsulation)

%% This is necessesary to run the script remotely using nohup
% addpath('../Helmholtz_data/')
% addpath('../Methods/')
% addpath('../NMLA')
% addpath('../Classes')
% addpath(genpath('../../ifem-ifem-8af722848372'))


%% Set up 

xc = 1/10; yc = 1/10;

% in this case we use the exact ray information given by the eikonal
% equatiion solved in a constant media
speed  = @(x) ones(size(x(:,1)));

%% Set up
plt = 0;                   % show solution or not
fquadorder = 10;            % numerical quadrature order

% size of the physical domain
a = 1/2;
% physical domain between -(a, a) to (a, a)
NPW = 8;

mult = 2;

% number of wavelenghts inside the domain
numberWavelengths = mult*16;

omega = 2*pi*numberWavelengths;
% discretization step (the constant needs to be modified)
h = 1/2^(round(log2((NPW*omega)/(2*pi))));

% physical width of the PML
wpml = 2*pi/omega + mult*h/NPW;
npml = round(ceil(wpml/h));
% maximum absorption of the PML
sigmaMax = 25/wpml;

% Number of subdomains
nSub = mult*10 ;

% computational domain [-a,a]^2 (we add the pml in this case)
a = a + wpml ;         

[node,elem] = squaremesh([-a,a,-a,a],h);

fprintf('Number of degrees of freedom is %d \n', size(node,1));

%% building the source 
k = omega;
sigma = 1/k;

source = @(p) -4*sqrt(k)*1i*1/(2*pi*sigma^2)*exp( -( (p(:,1)-xc).^2 + (p(:,2)-yc).^2  )/(2*sigma^2) );



%% Exact ray information
xx = node(:,1)-xc; yy = node(:,2)-yc;
ray = atan2(yy,xx);
ray = exp(1i*ray);
ray = ray.*(1 - (abs(xx)<eps).*(abs(yy)<eps));


N = size(node,1);       % number of grid points
Nray = size(ray,2);     % number of rays crossing at each grid node
Ndof = N*Nray; 

%% global problem to be solved

% initialazing the global model 
M0 = model;
M0.initRay(node,elem,omega,wpml,h/NPW, speed, ray, fquadorder);

% % performing LU factorization (only if we want ot compare against the true
% % solution 
%M0.LU_factorization()

% solving the PDE using the global model 
f = assemble_RHS_with_ray(node,elem,omega,source, speed,ray,fquadorder);
%v = M0.solve(f);

%% reconstructing the solution 
% This needs to be encapsulated for the ray solve 
% In this case the easiest way is to define a bigger class that
% encapsulates the layered stripping


if plt 
    grad = ray(:);
    grad = [real(grad),imag(grad)];
    repnode = repmat(node,Nray,1);
    temp = grad(:,1).*repnode(:,1) + grad(:,2).*repnode(:,2);

    c = speed(node);    % medium speed
    k = omega./c;       % wavenumber
    kk = repmat(k,1,Nray);

    u = v.*exp(1i*kk(:).*temp);
    u = reshape(u,N,Nray);
    u = sum(u,2);



    M0.showresult(real(u));
end

% getting only the interior degrees of freedom
fInt = f(M0.freeNode);

%% Local problems ( We need to find a way to perform the partitioning)
% This needs to be properly encapsulated
% geometrical information 
x = unique(node(:,1));
y = unique(node(:,2));

% gotta be careful with the PML (they have to be outside the 
% interior domain
xInt = x(npml:end-npml);
nInt = length(xInt);

% number of points per subdomains
nSubInt = round(nInt/nSub);

% the limits for the domaind decomposition
nIndLim = round(linspace(1,nInt-1,nSub+1));

% indices at which each subdomain starts (and finishes)
indn = npml + nIndLim(2:end);
ind1 = npml + nIndLim(1:end-1)+1;
ind1(1) = ind1(1) -1;

% we need to add this to the ray construction

% last and first coordinate
xn = x(indn);
x1 = x(ind1);

nodeArray = {};
elemArray = {};
for ii = 1:nSub
    [nodeArray{ii},elemArray{ii}] = squaremesh([x(ind1(ii)-npml),...
                                                x(indn(ii)+npml),...
                                                -a,a],h);
end

% we define the physical location of the interface boundaries
separatorN = xn(1:end-1);
separator1 = x1(2:end);

nx = size(x,1);
rayArray = {};
rayMatrix = (reshape(ray,nx,nx)).';

for ii = 1:nSub
     aux = (rayMatrix(ind1(ii)-npml:indn(ii)+npml,:)).';
     rayArray{ii} = aux(:);
end


% building the array of objects
MArray = {};
SubArray = {};
for ii = 1:nSub
    % we will need to make sure that the wpml has the correct width
    MArray{ii} = model;
    % we need to be carefull when defining the PML
    if ii == 1
        wpmlvec = [wpml, wpml-2*h, wpml, wpml];
    elseif ii == nSub
        wpmlvec = [wpml-2*h, wpml, wpml, wpml];
    else
        wpmlvec = [wpml-2*h, wpml-2*h, wpml, wpml];
    end     
    MArray{ii}.initRay(nodeArray{ii},elemArray{ii},omega,...
                    wpmlvec,h/NPW,speed,rayArray{ii},fquadorder);

end

%% we need to define all the local (and global) local Indices
% We start by defining the indices associated to each boundary (locally)
indxn = {};
indxnp = {};
indx0  = {};
indx1  = {};

for ii = 1:nSub
    if ii ~= 1
        indx0{ii} = find(MArray{ii}.node(MArray{ii}.freeNode,1) == separatorN(ii-1) );
        indx1{ii} = find(MArray{ii}.node(MArray{ii}.freeNode,1) == separator1(ii-1) );
    end
    
    if ii ~= nSub
        indxn{ii} = find(MArray{ii}.node(MArray{ii}.freeNode,1) == separatorN(ii) );
        indxnp{ii}= find(MArray{ii}.node(MArray{ii}.freeNode,1) == separator1(ii) );
    end
end
indxn{nSub} = [];
indxnp{nSub} = [];

% building the the local-global indices

indIntGlobal = {};
for ii = 1:nSub
    if ii == 1
        indIntGlobal{ii} = find(M0.node(M0.freeNode,1) <= separatorN(ii) );
    elseif ii == nSub
        indIntGlobal{ii} = find(M0.node(M0.freeNode,1) >= separator1(ii-1) );
    else
        indIntGlobal{ii} = find((M0.node(M0.freeNode,1) <= separatorN(ii) ).*  ...
                                (M0.node(M0.freeNode,1) >= separator1(ii-1)));
    end
end
  
 
% building the local-local indices

indIntLocal = {};
for ii = 1:nSub
    if ii == 1
        indIntLocal{ii} = find(MArray{ii}.node(MArray{ii}.freeNode,1) <= separatorN(ii) );
    elseif ii == nSub
        indIntLocal{ii} = find(MArray{ii}.node(MArray{ii}.freeNode,1) >= separator1(ii-1) );
    else
        indIntLocal{ii} = find((MArray{ii}.node(MArray{ii}.freeNode,1) <= separatorN(ii) ).*  ...
                                (MArray{ii}.node(MArray{ii}.freeNode,1) >= separator1(ii-1)));
    end
end

%% building the subdomains from the models and the local indices
for ii = 1:nSub
    if ii == 1
        position = 'N';
    elseif ii == nSub
        position = 'S';
    else
        position = 'C'; %center
    end
    
    SubArray{ii} = subdomain;    
    SubArray{ii}.init_Model(MArray{ii},indIntLocal{ii}, indIntGlobal{ii}, indxn{ii}, indxnp{ii}, indx0{ii}, indx1{ii}, position)
end

%% Factorizing the local matrices
tic();
for ii = 1:nSub
   SubArray{ii}.LU_factorization();
end
toc();

%% bulding the global array 
DDM = layer_stripping;

DDM.init_model(SubArray,M0, nSub);

%% Solving local problems

uArray =  DDM.solve_local(fInt);

fGamma = DDM.formRHS(uArray);

%% Testing the application of the matrix M
% 
% DDM.MBase.LU_factorization();
% uGlobal = DDM.MBase.solveInt(fInt);
% uGamma = DDM.extractGlobalTraces(uGlobal);
%  
% res =  DDM.applyM(uGamma) + fGamma;
% 
% norm(res)

%% testing the application of the  matrix MM
fGammaPol = DDM.formPolarizedRHS(uArray);

MMu = DDM.applyMpolarized(fGammaPol );

DDM.initPermutationMatrix();
MM = @(u) DDM.applyMpolarized(u);

%% testing the Dinv 
LL = @(u) DDM.applyGSinv(u);

tic();
uGmres = gmres(MM, -fGammaPol,[], 1e-6, 30, LL );
toc()
%% testing that the solution coincides

uBdy = uGmres(1:end/2) + uGmres(end/2+1:end);

res = DDM.applyMpolarized(uGmres) + fGammaPol;
norm(res)


%% reconstructing the solution within the subdomain
UsolArray = DDM.solveLocal(fInt, uBdy);
U = DDM.reconstruct(UsolArray);
 
UGlobal = zeros(size(f,1),1);
UGlobal(M0.freeNode) = U;


grad = ray(:);
grad = [real(grad),imag(grad)];
repnode = repmat(node,Nray,1);
temp = grad(:,1).*repnode(:,1) + grad(:,2).*repnode(:,2);

c = speed(node);    % medium speed
k = omega./c;       % wavenumber
kk = repmat(k,1,Nray);

u = UGlobal.*exp(1i*kk(:).*temp);
u = reshape(u,N,Nray);
u = sum(u,2);



% showing the real part of the solution
figure(1); 
M0.showresult(real(u));
