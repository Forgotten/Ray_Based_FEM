%% Driver for the standard finite element with absorbing boundary conditions 
% implemented via PML 
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
plt = 0;                   % show solution or not
fquadorder = 4;            % numerical quadrature order

% size of the physical domain
a = 1/2;
% physical domain between -(a, a) to (a, a)
NPW = 10;

omega = 2*pi*5;
% discretization step (the constant needs to be modified)
h = 0.1*1/2^(round(log2((NPW*omega)/(2*pi))));

% physical width of the PML
wpml = 2*pi/omega;
% maximum absoption of the PML
sigmaMax = 25/wpml;
    
fprintf('\ncase %d: \nomega/(2*pi) = %d,      1/h = %d \n\n', ii, omega/(2*pi), 1/h);

a = a + wpml;         % computational domian [-a,a]^2
[node,elem] = squaremesh([-a,a,-a,a],h);

% initialazing the model 
M0 = model;
M0.init(node,elem,omega,wpml,sigmaMax, pde,fquadorder);

% solving the PDE
v = M0.solve(pde.f);
    
% plotting the solution 
if plt 
    M0.showresult(real(v));
end