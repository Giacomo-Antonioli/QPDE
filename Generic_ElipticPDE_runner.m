%% Main Script to Test GenericQuantumDiffusion
clear; clc; close all;

%% 1. Simulation Configuration
dim = 2;            % Set to 2 or 3 to test different dimensions
n  = 6;
N=2^n;% Grid points (Power of 2, e.g., 16, 32, 64)
dt  = 1e-3;         % Time step
steps = 50;         % Number of time steps

%% 2. Define Diffusion Tensor (Matrix A)
% Must be d x d. Using Identity for simple diffusion, or a coupled matrix.
if dim == 2
    % Example 2D: Anisotropic diffusion
    A = eye(2);
elseif dim == 3
    % Example 3D: Simple isotropic diffusion
    A = eye(3); 
end

%% 3. Define Source Term (f) and Initial Condition (u_init)
% These must accept 'd' arguments (x,y) or (x,y,z)

if dim == 2
% Source Term: f(x,y) = cos(2*pi*x) * sin(-4*pi*y)
    f_handle = @(x,y) cos(2*pi*x) .* sin(-4*pi*y);
    
 
elseif dim == 3
    % --- 3D Functions ---
    % Source: Oscillating source in center
    f_handle = @(x,y,z) 5 * sin(2*pi*x) .* sin(2*pi*y) .* sin(2*pi*z);
    
 
end

%% 4. Run the Generalized Solver
fprintf('Starting %dD Simulation with N=%d...\n', dim, N);

GenericElliptic_QPDE(f_handle, A, N, dim)

fprintf('Done.\n');