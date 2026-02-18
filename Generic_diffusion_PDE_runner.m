%% Main Script to Test GenericQuantumDiffusion
clear; clc; close all;

%% 1. Simulation Configuration
dim = 3;            % Set to 2 or 3 to test different dimensions
n  = 3;
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
    
    % Initial Condition: u0(x,y) = cos(2*pi*x)sin(8*pi*y) + 2sin(6*pi*y) + 3sin(10*pi*x)cos^2(12*pi*y)
    u_handle = @(x,y) cos(2*pi*x) .* sin(8*pi*y) + ...
                      2 * sin(6*pi*y) + ...
                      3 * sin(10*pi*x) .* (cos(12*pi*y).^2);

elseif dim == 3
    % --- 3D Functions ---
    % Source: Oscillating source in center
    f_handle = @(x,y,z) 5 * sin(2*pi*x) .* sin(2*pi*y) .* sin(2*pi*z);
    
    % Initial: Gaussian ball in the center
    u_handle = @(x,y,z) exp( -((x-0.5).^2 + (y-0.5).^2 + (z-0.5).^2) / 0.05 );
end

%% 4. Run the Generalized Solver
fprintf('Starting %dD Simulation with N=%d...\n', dim, N);

% Make sure 'GenericDiffusion_QPDE.m' and 'QPDE_Generator.m' are in your path
GenericDiffusion_QPDE(f_handle, u_handle, A, N, dim, dt, steps);

fprintf('Done.\n');