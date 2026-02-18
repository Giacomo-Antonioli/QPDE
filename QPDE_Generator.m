function totalMat = QPDE_Generator(A,n)
d=size(A,1);
N=2^n;
FG=GroupFourier(d,n);
GF=FG.ctranspose();

% Build derivative operator in spectral domain
D = diag(spectral_eigenvalues(N));
% Identity matrix
I = eye(N);

OP=sparse(N^d,N^d);

  for i = 1:d
        for j = 1:d

            % initialize Kronecker product
            K = 1;

            for k = 1:d
                if (k == i) && (k == j)
                    % second derivative in direction k
                    K = kron(K, D * D);
                elseif (k == i) || (k == j)
                    % first derivative in direction k
                    K = kron(K, D);
                else
                    K = kron(K, I);
                end
            end

            % accumulate contribution
            OP = OP+ A(i,j) * K;

        end
    end
OP(1,1)=1;

invOP=diag(1./diag(OP));


DiagEncoding=MakeUnitary(invOP);

totalCircuit=qclab.QCircuit(d*n+1);
totalCircuit.push_back(GF);
totalCircuit.push_back(qclab.qgates.MatrixGate(0:d*n,DiagEncoding,"Diagonal"))
totalCircuit.push_back(FG);

totalMat=totalCircuit.matrix;
totalMat=totalMat(1:2^(d*n),1:2^(d*n));

end