function totalMat = QPDE_Generator_Diffusion(A, n, dt)
% QPDE_GENERATOR_DIFFUSION Generates the quantum diffusion operator.
%
%   U = (I - dt * Div(A * Grad))^-1
%
% Inputs:
%   A  : Diffusion coefficient matrix (d x d)
%   n  : Number of qubits per dimension (N = 2^n)
%   dt : Time step size

    d = size(A, 1);
    N = 2^n;

    % --- 1. Gates (QFT and iQFT) ---
    FG = GroupFourier(d, n);   % Forward QFT
    GF = FG.ctranspose();      % Inverse QFT

    % --- 2. Build Spectral Derivatives ---
    % Eigenvalues of the derivative operator d/dx in spectral space:
    % D_k = 2*pi*i*k / L  (Assuming L=1 for normalized operator)
    % We create a vector of these eigenvalues.
    k_vec = [0:N/2-1, -N/2:-1]'; 
    D_vals = 2i * pi * k_vec; 
    
    % Identity vector for dimensions not being differentiated
    I_vals = ones(N, 1);

    % --- 3. Construct Laplacian Eigenvalues via Kronecker Products ---
    % We accumulate the sum: L_hat = sum_{i,j} A_{ij} * (k_i * k_j)
    % We use vectors instead of full matrices to save memory.
    OP_vals = zeros(N^d, 1);

    for i = 1:d
        for j = 1:d
            if A(i,j) == 0
                continue;
            end
            
            % Initialize accumulator for the current term
            % We will perform kron(current, K) to match MATLAB's 
            % column-major ordering (Dimension 1 is fastest/innermost).
            term_vals = 1; 
            
            % Loop d down to 1 ensures dim 1 is innermost in the Kron product
            for k = d:-1:1
                if (k == i) && (k == j)
                    % Second derivative: d^2 / dx_k^2
                    current_op = D_vals.^2;
                elseif (k == i) || (k == j)
                    % First derivative: d / dx_k
                    current_op = D_vals;
                else
                    % Identity
                    current_op = I_vals;
                end
                
                % Accumulate: New dimension goes to the "left" (slower index)
                term_vals = kron(current_op, term_vals);
            end
            
            % Add weighted term to total operator eigenvalues
            OP_vals = OP_vals + A(i,j) * term_vals;
        end
    end
    
    % --- 4. Build Diffusion Propagator ---
    % Implicit Euler Step: (I - dt * Laplacian)^-1
    % In spectral space: 1 ./ (1 - dt * OP_vals)
    Inv_OP_vals = 1 ./ (1 - dt * OP_vals);
    
    % Normalize for Unitary Encoding (Must be <= 1)
    alpha = max(abs(Inv_OP_vals));
    Normalized_Vals = Inv_OP_vals / alpha;
    
    % Create sparse diagonal matrix
    DiagMat = sparse(1:N^d, 1:N^d, Normalized_Vals, N^d, N^d);
    
    % Encode into Unitary
    % (Using MakeUnitary helper assumed to be available)
    DiagEncoding = MakeUnitary(DiagMat);

    % --- 5. Circuit Assembly ---
    % Sequence: QFT (Spatial->Spectral) -> Diagonal Op -> iQFT (Spectral->Spatial)
    totalCircuit = qclab.QCircuit(d*n + 1);
    totalCircuit.push_back(FG); 
    totalCircuit.push_back(qclab.qgates.MatrixGate(0:d*n, DiagEncoding, "Diagonal"));
    totalCircuit.push_back(GF);

    % Extract Matrix
    totalMat = totalCircuit.matrix;
    totalMat = totalMat(1:2^(d*n), 1:2^(d*n));
end