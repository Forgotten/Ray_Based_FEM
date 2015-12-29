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
NPW = 10;

omega = 2*pi*6;
% discretization step (the constant needs to be modified)
h = 1/2^(round(log2((NPW*omega)/(2*pi))));

% physical width of the PML
wpml = 2*pi/omega;
npml = round(ceil(wpml/h));
% maximum absoption of the PML
sigmaMax = 25/wpml;

% Number of subdomains
nSub = 6;


a = a + wpml ;         % computational domian [-a,a]^2


%% global problem to be solved
[node,elem] = squaremesh([-a,a,-a,a],h);


% initialazing the model 
M0 = model;
M0.init(node,elem,omega,wpml,sigmaMax, pde,fquadorder);

% performing LU factorization
M0.LU_factorization()

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
ind1(1) = ind1(1) -1

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

%% Preliminary tests
% testing the continuity of the discrete operator within each subproblem
for ii = 2:nSub

    fprintf('misfit is given by %d \n',...
            norm(max(abs(MArray{ii}.H(indx0{ii},indx1{ii}) ...
            - MArray{ii-1}.H(indxn{ii-1},indxnp{ii-1}))))) ;
end

for ii = 1:nSub-1

    fprintf('misfit is given by %d \n',...
            norm(max(abs(MArray{ii+1}.H(indx0{ii+1},indx1{ii+1}) ...
            - MArray{ii}.H(indxn{ii},indxnp{ii}))))) ;
end


%% building the local-global indices

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
  
  
%% source partitioning

fInt = f(M0.freeNode);
fIntLocal = {};

for ii = 1:nSub
    fIntLocal{ii} = fInt(indIntGlobal{ii});
end

% testing the the parititioning is done properly

fTest = [];
for ii = 1:nSub
    fTest = [fTest; fIntLocal{ii}];
end

fprintf('misfit between the paritionned source %d \n', norm(fTest - fInt))

% building the local solutions 

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

%% preconditioner 


fInt = f(M0.freeNode);
fIntLocal = {};
fExtLocal = {};

for ii = 1:nSub
    fIntLocal{ii} = zeros(size(MArray{ii}.node(MArray{ii}.freeNode,1),1),1);
    fIntLocal{ii}(indIntLocal{ii}) = fInt(indIntGlobal{ii});
end

%% Factorizing the local matrices

for ii = 1:nSub
   MArray{ii}.LU_factorization();
end

%Defining the local traces 

% size of the local traces
n = size(indxn{3},1);

u_0  = zeros(n*nSub,1);
u_1  = zeros(n*nSub,1);
u_n  = zeros(n*nSub,1);
u_np = zeros(n*nSub,1);

index = 1:n;


%% Solving for each subdomain 
uArray = {};
for ii = 1:nSub
    uArray{ii} = MArray{ii}.solveInt(fIntLocal{ii});
    localSizes(ii) = size(indIntGlobal{ii},1);
end

localLim = [0 cumsum(localSizes)];

%% local solve + downwards sweep
for ii = 1:nSub
   
        % making a local copy of the local rhs
    rhsLocaltemp = fIntLocal{ii};

    if ii ~=1
        rhsLocaltemp(indx0{ii}) =  rhsLocaltemp(indx0{ii}) +...
                   MArray{ii}.H(indx0{ii},indx1{ii})*u_np((ii-2)*n + index);
        rhsLocaltemp(indx1{ii}) = rhsLocaltemp(indx1{ii}) -...
                    MArray{ii}.H(indx1{ii},indx0{ii})*u_n((ii-2)*n + index);
    end

        % solving the rhs
        vDown = MArray{ii}.solveInt(rhsLocaltemp);

        % extracting the traces
        if ii ~= nSub
            u_n((ii-1)*n  + index) = vDown(indxn{ii});
            u_np((ii-1)*n + index) = vDown(indxnp{ii});
        end
end

%% upwards sweep + reflections + reconstruction

uPrecond = zeros(size(fInt,1),1);

     for ii = nSub:-1:1
      
        % making a copy of the parititioned source
         rhsLocaltemp = fIntLocal{ii};

        % adding the source at the boundaries
        if ii~= 1
            % we need to be carefull at the edges
            rhsLocaltemp(indx1{ii})  = rhsLocaltemp(indx1{ii}) - ...
               MArray{ii}.H(indx1{ii},indx0{ii})*u_n((ii-2)*n + index);
            rhsLocaltemp(indx0{ii})  =  rhsLocaltemp(indx0{ii}) ...
                + MArray{ii}.H(indx0{ii},indx1{ii})*u_np((ii-2)*n + index);
        end
        if ii~= nSub
            rhsLocaltemp(indxnp{ii}) = rhsLocaltemp(indxnp{ii}) + ...
                MArray{ii}.H(indxnp{ii},indxn{ii})*u_0((ii)*n + index);
            rhsLocaltemp(indxn{ii})  = rhsLocaltemp(indxn{ii}) - ...
                MArray{ii}.H(indxn{ii},indxnp{ii})*u_1((ii)*n + index);
        end
 
        % solving the local problem
        uUp = MArray{ii}.solveInt(rhsLocaltemp);

        if ii > 1
            u_0((ii-1)*n + index) = uUp(indx0{ii});
            u_1((ii-1)*n + index) = uUp(indx1{ii}) - u_np((ii-2)*n + index);
        end

        % reconstructing the problem on the fly
        uPrecond(localLim(ii)+1:localLim(ii+1)) = uUp(indIntLocal{ii});
    end

%% Padding the preconditioned solution 

UPrecond = zeros(size(v,1),1);
UPrecond(M0.freeNode) = uPrecond;
norm(v - UPrecond)
