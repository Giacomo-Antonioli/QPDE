%% Quantum Generic Diffusion Solver with Detailed Plots
clear; clc; close all;

%% 1. Setup & Parameters
n = 5;              % Qubits per dimension (N = 32)
N = 2^n;
dt = 1e-3;
steps = 40;
dim = 2;            % Dimension of the problem (2D)

% Physical Grid [0, 1)
L = 1;
dx = L / N;
x = (0:N-1) * dx; 
[X, Y] = ndgrid(x, x);

% --- GENERIC SQUARE MATRIX A ---
A = [3.0, 1.0; 
     1.0, 2.0];

% Source Term (f) and Initial Condition (u)
f = cos(2*pi*X) .* sin(-2*pi*Y);
u_init = cos(2*pi*X) .* sin(8*pi*Y) + 2*sin(6*pi*Y) + 3*sin(10*X).*cos(12*Y).^2;

%% 2. Operator Construction (Pre-computation)
fprintf('Building Quantum Operator for Generic A...\n');

% Wave Vectors
k_ind = [0:N/2-1, -N/2:-1]';
k_base = 2i * pi * k_ind;       
[K_grids{1:dim}] = ndgrid(k_base);

% Construct Elliptic Operator
Elliptic_diag = zeros(size(K_grids{1}));
for i = 1:dim
    for j = 1:dim
        Elliptic_diag = Elliptic_diag + A(i,j) * K_grids{i} .* K_grids{j};
    end
end

% Create Filter
Denominator = 1 - dt * Elliptic_diag;
Filter = 1 ./ Denominator;

% Quantum Encoding
alpha = max(abs(Filter(:)));           
DiagGate = MakeUnitary(diag(Filter(:) / alpha));

% Circuit: QFT -> Diagonal -> iQFT
FG = GroupFourier(dim, n);               
GF = FG.ctranspose();                  

QC = qclab.QCircuit(dim*n + 1);
QC.push_back(FG);
QC.push_back(qclab.qgates.MatrixGate(0:dim*n, DiagGate, "D"));
QC.push_back(GF);

% Extract Matrix
FullMat = QC.matrix;
Q_Op = FullMat(1:N^dim, 1:N^dim);

%% 3. Time Evolution Loop
u_class = u_init;
u_quant = u_init;

energy_hist_class = zeros(1, steps+1);
energy_hist_quant = zeros(1, steps+1);

calc_energy = @(u_in) 0.5 * real(sum(sum( conj(fft2(u_in)) .* (-Elliptic_diag) .* fft2(u_in) ))) / N^(2*dim);

energy_hist_class(1) = calc_energy(u_class);
energy_hist_quant(1) = calc_energy(u_quant);

fprintf('Running %d steps (N=%d)...\n', steps, N);

for t = 1:steps
    % Classical
    u_h = fftn(u_class);
    f_h = fftn(f);
    u_class = real(ifftn((u_h - dt * f_h) .* Filter));

    % Quantum
    v = u_quant - dt*f;
    u_vec_next = (Q_Op * v(:)) * alpha;
    u_quant = reshape(real(u_vec_next), repmat(N, 1, dim));
    
    energy_hist_class(t+1) = calc_energy(u_class);
    energy_hist_quant(t+1) = calc_energy(u_quant);
end

%% 4. Visualization
diff_norm = norm(u_class(:) - u_quant(:), 'fro');
abs_err_map = abs(u_class - u_quant);
rel_err_map = abs_err_map ./ (abs(u_class) + 1e-10); % Avoid division by zero
energy_err = abs(energy_hist_class - energy_hist_quant);

fprintf('Final Frobenius Error: %.3e\n', diff_norm);

figure('Position', [50, 50, 1600, 800], 'Name', 'Generic Diffusion Results');

% --- PLOT 1: Energy Evolution ---
subplot(2, 4, 1);
plot(0:steps, energy_hist_class, 'b-', 'LineWidth', 2); hold on;
plot(0:steps, energy_hist_quant, 'r--', 'LineWidth', 2);
title('Energy Evolution');
xlabel('Time Step'); ylabel('Energy');
legend('Classical', 'Quantum'); grid on;

% --- PLOT 2: Absolute Energy Error ---
subplot(2, 4, 2);
semilogy(0:steps, energy_err, 'k-o', 'LineWidth', 1.5, 'MarkerSize', 4);
title('Abs. Energy Error |E_{cl} - E_{qu}|');
xlabel('Time Step'); ylabel('Error'); grid on;

% --- PLOT 3: Classical Solution ---
subplot(2, 4, 3);
imagesc(u_class'); 
axis square; colorbar;
title('Classical Solution');
xlabel('x'); ylabel('y'); set(gca, 'YDir', 'normal');

% --- PLOT 4: Quantum Solution ---
subplot(2, 4, 4);
imagesc(u_quant'); 
axis square; colorbar;
title('Quantum Solution');
xlabel('x'); ylabel('y'); set(gca, 'YDir', 'normal');

% --- PLOT 5: Absolute Error Map ---
subplot(2, 4, 5);
imagesc(abs_err_map'); 
axis square; colorbar;
title(sprintf('Absolute Error\nMax: %.2e', max(abs_err_map(:))));
xlabel('x'); ylabel('y'); set(gca, 'YDir', 'normal');

% --- PLOT 6: Relative Error Map ---
subplot(2, 4, 6);
imagesc(rel_err_map'); 
axis square; colorbar;
title('Relative Error');
xlabel('x'); ylabel('y'); set(gca, 'YDir', 'normal');

% --- PLOT 7: Cross-Section (Middle of Y) ---
subplot(2, 4, 7);
mid_idx = floor(N/2);
plot(x, u_class(:, mid_idx), 'b-', 'LineWidth', 2); hold on;
plot(x, u_quant(:, mid_idx), 'r--', 'LineWidth', 2);
title(sprintf('Cross-section at y=%.2f', x(mid_idx)));
legend('Classical', 'Quantum'); grid on;
xlabel('x'); ylabel('u(x, y_{mid})');

% --- PLOT 8: Error Metrics (Text Box) ---
subplot(2, 4, 8);
axis off;
text(0.1, 0.8, sprintf('Run Metrics (Steps=%d)', steps), 'FontWeight', 'bold', 'FontSize', 12);
text(0.1, 0.6, sprintf('Frobenius Norm: %.3e', diff_norm));
text(0.1, 0.5, sprintf('Max Abs Error: %.3e', max(abs_err_map(:))));
text(0.1, 0.4, sprintf('Mean Abs Error: %.3e', mean(abs_err_map(:))));
text(0.1, 0.3, sprintf('Final Energy (Cl): %.4f', energy_hist_class(end)));
text(0.1, 0.2, sprintf('Final Energy (Qu): %.4f', energy_hist_quant(end)));
text(0.1, 0.1, sprintf('Encoding Alpha: %.4f', alpha));
box on;