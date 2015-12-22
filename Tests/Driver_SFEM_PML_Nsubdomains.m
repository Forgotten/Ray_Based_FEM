%% Driver for the standard finite element with absorbing boundary conditions 
% implemented via PML using N different subdomains
% The main aim of this small script is to provided some details on how the
% domains are glued together.
%clc; clear;


%% Load source data
pde = Helmholtz_data8;

global omega;
global a;
global xc;
global yc;
global h;

xc = 1/8;   yc = 1/10;

%% Set up
plt = 1;                   % show solution or not
fquadorder = 4;            % numerical quadrature order

% size of the physical domain
a = 1/2;
% physical domain between -(a, a) to (a, a)
NPW = 5;

omega = 2*pi*3;
% discretization step (the constant needs to be modified)
h = 1/2^(round(log2((NPW*omega)/(2*pi))));

% physical width of the PML
wpml = 2*pi/omega;
npml = round(ceil(wpml/h));
% maximum absoption of the PML
sigmaMax = 25/wpml;

% Number of subdomains
nSub = 4;


a = a + wpml ;         % computational domian [-a,a]^2


%% global problem to be solved
[node,elem] = squaremesh([-a,a,-a,a],h);


% initialazing the model 
M0 = model;
M0.init(node,elem,omega,wpml,sigmaMax, pde,fquadorder);

% solving the PDE
f = assemble_RHS(node,elem,pde.f,fquadorder);
v = M0.solve(f);
    
% plotting the solution 
if plt 
    M0.showresult(real(v));
end

%% Local problems ( We need to find a way to perform the partitioning)
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

% we suppose at each interface we have 
% separator1
% separatorN

MArray = {};

for ii = 1:nSub
    % we will need to make sure that the wpml has the correct width
    MArray{ii} = model;
    MArray{ii}.init(nodeArray{ii},elemArray{ii},omega,...
                    wpml-h,sigmaMax,pde,fquadorder);

end

% we need to define all the local (and global) local Indices


% defining local model (i.e. local subdomains this would need to be
% properly encapsulated 
M1 = model;
M1.init(node1,elem1,omega,wpml,sigmaMax, pde,fquadorder);

% indices of the boundary points
ind1xn  = find(M1.node(M1.freeNode,1) == separator0 );
ind1xnp = find(M1.node(M1.freeNode,1) == separator1 );

% indices of the local interior points
ind1Int = find(M1.node(M1.freeNode,1) <= separator0 );


M2 = model;
M2.init(node2,elem2,omega,wpml,sigmaMax, pde,fquadorder);

ind2x0  = find(M2.node(M2.freeNode,1) == separator0 );
ind2x1  = find(M2.node(M2.freeNode,1) == separator1 );
ind2Int = find(M2.node(M2.freeNode,1) >= separator1 );

H1 = M1.H;
H2 = M2.H;

%% testing the reconstruction
% interior souce
fInt = f(M0.freeNode);
% interior solution
vInt = v(M0.freeNode);

v0 = vInt(find(M0.node(M0.freeNode,1) == separator0 ))
v1 = vInt(find(M0.node(M0.freeNode,1) == separator1 ))

indInt1 = find(M0.node(M0.freeNode,1) <= separator0 );
indInt2 = find(M0.node(M0.freeNode,1) >= separator1);

%(we need to obtain the different set of indices to make the correct trnasition)

f1 = zeros(length(M1.freeNode),1);
f2 = zeros(length(M2.freeNode),1);

f1(ind1Int) = fInt(indInt1);
f2(ind2Int) = fInt(indInt2);

% adding the GRF for the first su
f1(ind1xn) =  f1(ind1xn)  - M1.H(ind1xn ,ind1xnp)*v1;
f1(ind1xnp) = f1(ind1xnp) + M1.H(ind1xnp,ind1xn )*v0;

u1 = M1.solveInt(f1);

uu1 = zeros(size(M1.node,1),1);
uu1(M1.freeNode) = u1;


figure(2); clf(); 
M1.showresult(real(uu1));

f2(ind2x0) = f2(ind2x0) + M2.H(ind2x0,ind2x1)*v1;
f2(ind2x1) = f2(ind2x1) - M2.H(ind2x1,ind2x0 )*v0;

u2 = M2.solveInt(f2);

uu2 = zeros(size(M2.node,1),1);
uu2(M1.freeNode) = u2;

figure(3); clf(); 
M2.showresult(real(uu2));

%% adding the local solutions together

U = zeros(size(M0.freeNode,1),1);

U(indInt1) = u1(ind1Int);
U(indInt2) = u2(ind2Int);

UU = zeros(size(M0.node,1),1);
UU(M0.freeNode) = U;
figure(4); clf();
M0.showresult(real(UU))

fprintf('Error in the reconstruction is %d \n', norm(v - UU)/norm(v))