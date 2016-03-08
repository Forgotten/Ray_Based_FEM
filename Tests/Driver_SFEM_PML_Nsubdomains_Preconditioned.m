%% Driver for the standard finite element with absorbing boundary conditions 
% implemented via PML using N different subdomains
% We test the preconditioner (without encapsulation)
% addpath('../Helmholtz_data/')
% addpath('../Methods/')
% addpath('../NMLA')
% addpath('../Classes')
% addpath(genpath('../../ifem-ifem-8af722848372'))


%% Load source data
% we want to use a non-constant wave-speed
pde = Helmholtz_data9;

global omega;
global a;
global xc;
global yc;
global h;

xc = 1/8;   yc = 1/10;

%% Set up
plt = 0;                   % show solution or not
fquadorder = 4;            % numerical quadrature order

% size of the physical domain
a = 1/2;
% physical domain between -(a, a) to (a, a)
NPW = 10;

% number of wavelenghts inside the domain
numberWavelengths = 3*46;

omega = 2*pi*numberWavelengths;
% discretization step (the constant needs to be modified)
h = 1/2^(round(log2((NPW*omega)/(2*pi))));

% physical width of the PML
wpml = 2*pi/omega + 4*h;
npml = round(ceil(wpml/h));
% maximum absoption of the PML
sigmaMax = 25/wpml;

% Number of subdomains
nSub =45;

a = a + wpml ;         % computational domian [-a,a]^2


%% global problem to be solved
[node,elem] = squaremesh([-a,a,-a,a],h);


% initialazing the  global model 
M0 = model;
M0.init(node,elem,omega,wpml,sigmaMax, pde,fquadorder);

% performing LU factorization
%M0.LU_factorization()

% solving the PDE using the global model 
f = assemble_RHS(node,elem,pde.f,fquadorder);
%v = M0.solve(f);
    
% plotting the global solution 
if plt 
    M0.showresult(real(v));
end

% getting only the interior degrees of freedom
fInt = f(M0.freeNode);

%% Local problems ( We need to find a way to perform the partitioning)
%  This needs to tbe properly encapsulated
% geometrical information 
x = unique(node(:,1));
y = unique(node(:,2));

% gotta be careful with the PML (they have to be outside the 
% interior domain
xInt = x(npml:end-npml);
nInt = length(xInt);

nSubInt = round(nInt/nSub);

nIndLim = round(linspace(1,nInt-1,nSub+1));

% indices at which each subdomain starts (and finishes)
indn = npml + nIndLim(2:end);
ind1 = npml + nIndLim(1:end-1)+1;
ind1(1) = ind1(1) -1;

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

% building the array of objects
MArray = {};

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
    MArray{ii}.init(nodeArray{ii},elemArray{ii},omega,...
                    wpmlvec,sigmaMax,pde,fquadorder);

end

% we need to define all the local (and global) local Indices

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

%% Factorizing the local matrices
tic();
for ii = 1:nSub
   MArray{ii}.LU_factorization();
end
toc();
%% Testing the Preconditioner
% we need to wrap and encapsulate as much as possible
%wrapper = partitioning(MArray, indx0, indx1, indxn ,indxnp, indIntGlobal, indIntLocal);

% defining the wrapped preconditioner function                        
%GSPrecond = @(x) GSpreconditioner(MArray, indx0, indx1, indxn ,indxnp, indIntGlobal, indIntLocal, x );

% optimized version of the preconditioner (one less for loop) and some
% modifications on the Array
GSPrecondOpt = @(x) GSpreconditionerOpt(MArray, indx0, indx1, indxn ,indxnp, indIntGlobal, indIntLocal, x );

% solving the linear sysmtem with GMRES
%tic();
%uGmres = gmres(M0.H, fInt,[], 1e-6, 10, GSPrecond );
%toc()

tic();
uGmres = gmres(M0.H, fInt,[], 1e-6, 10, GSPrecondOpt );
toc()


%% Padding the preconditioned solution 

UPrecond = zeros(size(f,1),1);
UPrecond(M0.freeNode) = uGmres;

% showing the real part of the solution
figure(1); 
M0.showresult(real(UPrecond));