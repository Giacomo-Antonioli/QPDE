function GenericDiffusion_QPDE(f_handle, u_init_handle, A, N, d, dt, steps,f_analytical)
% GENERICDIFFUSION_QPDE Generic Quantum Diffusion Solver (2D/3D).
% Uses implicit time-stepping for classical verification.
 flagUtrue=true;    

 disp(nargin)
    if nargin < 8

           flagUtrue=false;
    end
    %% 1. Setup & Grid Generation
    L = 1; 
    dx = L / N;
    x_vec = (0:N-1) * dx;
 
    % Generate Grids
    ndgrid_inputs = repmat({x_vec'}, 1, d);
    grids = cell(1, d);
    [grids{:}] = ndgrid(ndgrid_inputs{:});
    
    % Evaluate Source (f) and Initial Condition (u)
    f_vals = f_handle(grids{:});
    u_init = u_init_handle(grids{:});
    
    % N vector for solver (e.g. [32, 32])
    N_vecs = repmat(N, 1, d);

    %% 2. Operator Construction
    n = log2(N);
    fprintf('Building Quantum Operator for d=%d, N=%d...\n', d, N);
    
    % A) Quantum Operator (Propagator)
    % Ensure QPDE_Generator.m is in your path
    Q_Op = QPDE_Generator_Diffusion(A, n, dt);
    
    % B) Compute Normalization Alpha
    % We compute alpha by running the solver once on a dummy impulse or 
    % by estimating the max eigenvalue of the inverse operator.
    % For diffusion (I - dt*L)^-1, the max value is 1 (at k=0).
    alpha = 1.0; 

    %% 3. Time Evolution Loop
    u_class = u_init;
    u_quant = u_init;
    % Initialize Energy History
    E_class_history = zeros(1, steps);
    E_quant_history = zeros(1, steps);
    fprintf('Running %d steps...\n', steps);
    
    for t = 1:steps
        % --- Classical Step ---
        % Solving: (I - dt * Div(A * Grad)) u_new = u_old - dt * f
        rhs = u_class - dt * f_vals;
        u_class = solver_Diffusion_generic(rhs, grids, A, N_vecs, dx, dt);
        
        % --- Quantum Step ---
        % 1. Apply Source: v = u - dt * f
        v = u_quant - dt * f_vals; 
        
        % 2. Apply Operator: u_new = Q_Op * v
        u_vec_next = (Q_Op * v(:)) * alpha;
        
        % 3. Reshape
        u_quant = reshape(real(u_vec_next), N_vecs);

        % --- Compute Energy for Both ---
        E_class_history(t) = energy(u_class, A, N, dx,d);
        E_quant_history(t) = energy(u_quant, A, N, dx,d);
    end
    disp(flagUtrue)
    if flagUtrue
        u_true=f_analytical(grids{:});
        visualize_simulation_results(u_class, u_quant, d, N,u_true);
    else
        visualize_simulation_results(u_class, u_quant, d, N);
    end
    %% 4. Visualization & Metrics
    
    
    
    visualize_energy(E_class_history, E_quant_history,dt)

   timestamp = char(datetime('now','Format','yyyyMMdd_HHmmSS'));
folderName = sprintf('Results/Results_d%d_%s', d, timestamp);

if ~exist(folderName, 'dir')
    mkdir(folderName);
end

h5FileName = fullfile(folderName, 'simulation_data.h5');

% --- Save Solutions ---
if d == 3
    u_class_save = double(reshape(real(u_class), N, N, N));
    u_quant_save = double(reshape(real(u_quant), N, N, N));
else
    u_class_save = double(reshape(real(u_class), N, N));
    u_quant_save = double(reshape(real(u_quant), N, N));
end

h5create(h5FileName, '/u_classical', size(u_class_save));
h5write(h5FileName,  '/u_classical', u_class_save);

h5create(h5FileName, '/u_quantum', size(u_quant_save));
h5write(h5FileName,  '/u_quantum', u_quant_save);

% --- Save Energy Histories ---
h5create(h5FileName, '/E_class_history', size(E_class_history));
h5write(h5FileName,  '/E_class_history', double(E_class_history));

h5create(h5FileName, '/E_quant_history', size(E_quant_history));
h5write(h5FileName,  '/E_quant_history', double(E_quant_history));

% --- Save Ground Truth if available ---
if flagUtrue
    if d == 3
        gt_save = double(reshape(real(u_true), N, N, N));
    else
        gt_save = double(reshape(real(u_true), N, N));
    end
    h5create(h5FileName, '/u_true', size(gt_save));
    h5write(h5FileName,  '/u_true', gt_save);
    fprintf('Saved u_classical, u_quantum, u_true, and energy histories to %s\n', folderName);
else
    fprintf('Saved u_classical, u_quantum, and energy histories to %s\n', folderName);
end

% --- Relative Error Calculation ---
u_class_flat = u_class_save(:);
u_quant_flat = u_quant_save(:);

if flagUtrue
    gt_flat = gt_save(:);
    rel_err_classical = norm(u_class_flat - gt_flat) / norm(gt_flat);
    rel_err_quantum   = norm(u_quant_flat - gt_flat) / norm(gt_flat);

    fprintf('\n--- Relative Errors (vs Ground Truth) ---\n');
    fprintf('Classical Relative Error: %.4e\n', rel_err_classical);
    fprintf('Quantum   Relative Error: %.4e\n', rel_err_quantum);
else
    rel_err_quantum = norm(u_quant_flat - u_class_flat) / norm(u_class_flat);

    fprintf('\n--- Relative Error (vs Classical Solution) ---\n');
    fprintf('Quantum Relative Error: %.4e\n', rel_err_quantum);
end

% --- Energy History Summary ---
fprintf('\n--- Energy History Summary ---\n');
fprintf('Classical: initial = %.4e,  final = %.4e\n', E_class_history(1), E_class_history(end));
fprintf('Quantum:   initial = %.4e,  final = %.4e\n', E_quant_history(1), E_quant_history(end));
end