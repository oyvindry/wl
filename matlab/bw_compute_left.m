function [W,Wtilde, A_pre_inv, Atilde_pre_inv, N0, C, Ctilde]=bw_compute_left(g0, g1, N, K, g0tilde, g1tilde, Ntilde, Ktilde, symbolic)
    normalize = 1;
    isortho = (isequal(g0, g0tilde) & isequal(g1, g1tilde));  
    Nprime = max(N,Ntilde);
    [L,R]=findsupports(g0);
    [Ltilde,Rtilde]=findsupports(g0tilde);
    N0=max( [ceil((2*N-K-1+Rtilde)/2) ceil((2*Ntilde-Ktilde-1+R)/2) R Rtilde]); % Equation (5.1)
    
    S = [1:(K-N+R) (K-N+R+2):2:(K-N-R+2*N0)]; S = S'; % Defined in Lemma 5.2
    Stilde = [1:(Ktilde-Ntilde+Rtilde) (Ktilde-Ntilde+Rtilde+2):2:(Ktilde-Ntilde-Rtilde+2*N0)]; Stilde = Stilde';
    
    T      = N0 + floor(-(R+Ltilde+1)/2); % Equation (5.3)
    Ttilde = N0 + floor(-(Rtilde+L+1)/2);
    T = max(T, Nprime); % TODO: add to the paper.
    Ttilde = max(Ttilde, Nprime);

    % Find initial boundary functions
    G_I_X = Gsegment(g0, L:R, (-R+1):(K-1), 2*((-R+1):(K-1)), symbolic); % Defined in Lemma 4.1
    G_I_Z = Gsegment(g0, L:R, K:(2*K+R-2),  2*((-R+1):(K-1)), symbolic);
    [C,C_c]=findc(R, K, N, symbolic);
    Xe = inv(C'*C)*C'* G_I_X *C; % Equation (4.2). Replace inv(C'*C)*C' with transpose if matrix is orthogonal.
    Ze = G_I_Z*C;
    
    G_I_X = Gsegment(g0tilde, Ltilde:Rtilde, (-Rtilde+1):(Ktilde-1),     2*((-Rtilde+1):(Ktilde-1)), symbolic);
    G_I_Z = Gsegment(g0tilde, Ltilde:Rtilde, Ktilde:(2*Ktilde+Rtilde-2), 2*((-Rtilde+1):(Ktilde-1)), symbolic);
    [Ctilde,C_ctilde]=findc(Rtilde, Ktilde, Ntilde, symbolic);
    Xtildee = inv(Ctilde'*Ctilde)*Ctilde'*G_I_X *Ctilde;
    Ztildee = G_I_Z*Ctilde;
    
    % Make phi and phitilde staggered
    [Qmatr,Rmatr] = qr_exact((flipud(C((R+K-N):end,:)))', symbolic);
    P1 = fliplr(Qmatr);
    [Qmatr,Rmatr] = qr_exact((flipud(Ctilde((Rtilde+Ktilde-Ntilde):end,:)))', symbolic);
    P2 = fliplr(Qmatr);
    [Xe,Ze,C,Xtildee,Ztildee,Ctilde]=update_vars(Xe,Ze,C,Xtildee,Ztildee,Ctilde,P1,P2);
    
    % Bi-orthogonalize phi and phitilde
    
    if Ntilde == N
        [Ze,Ztildee]=expand_cols_smallest(Ze,Ztildee);
        rhs = Ze'*Ztildee;
    elseif Ntilde > N  % Address when Ntilde and N are different
        [Ze1,Ze2]=expand_cols_smallest( Ze( (Ntilde-N+1):end, :), Ztildee);
        [zefirst,Ctildelast]=expand_cols_smallest(Ze, Ctilde( (end-(Ntilde-N)+1):end, :));
        rhs2 = zefirst'*Ctildelast*Xtildee; 
        rhs = Ze1'*Ze2 + rhs2;
    else
        [Ze1,Ze2]=expand_cols_smallest( Ze, Ztildee( (N-Ntilde+1):end, :));
        [ztildeefirst,Clast]=expand_cols_smallest(Ztildee, C( (end-(N-Ntilde)+1):end, :)); % Check
        rhs2 = Xe'*Clast'*ztildeefirst,
        rhs = Ze1'*Ze2 + rhs2;
    end
    
    ls = eye(N*Ntilde) - kron(Xe',Xtildee');
    rs = reshape(rhs', N*Ntilde, 1);
    Y  = reshape(ls\rs, Ntilde, N); % Solve Equation (4.12)
    Y=Y'; % Gramm matrix found
    
    if isortho
        P1 = inv(chol(Y));
        P2 = P1;
    else
        [Lmatr,Umatr]= lu_nopivot(Y, symbolic);
        if Ntilde>N
            Umatr( (end-(Ntilde-N)+1):end, : ) = Ctilde( (end-(Ntilde-N)+1):end, :);
        elseif N>Ntilde
            Lmatr( :, (end-(N-Ntilde)+1):end ) = C( (end-(N-Ntilde)+1):end, : )';
        end
        
        P1 = (inv(Lmatr))';
        P2 = inv(Umatr);
    end    
    
    [Xe,Ze,C,Xtildee,Ztildee,Ctilde]=update_vars(Xe,Ze,C,Xtildee,Ztildee,Ctilde,P1,P2);
    
    % Mother wavelets
    
    % Construct Xo
    newcols = Gsegment(g0, (L:R) + K - N, 0:(2*T + K - N + R), 2*(N:T), symbolic);
    [Rmatr, newcols] = expand_cols_smallest([Xe; Ze], newcols);
    G = [Rmatr newcols];
    
    newcols = Gsegment(g0tilde, (Ltilde:Rtilde) + K - N, 0:(2*T + K - N + Rtilde), 2*(Ntilde:T), symbolic);
    [Rmatr, newcols] = expand_cols_smallest([Xtildee; Ztildee], newcols);
    Gtilde = [Rmatr newcols];
    
    lastmatr = G*(Gtilde(S,:))';
    idmatr = eye(S(N0+K-N));
    idmatr = idmatr(:,S);
    [idmatr,lastmatr]=expand_cols_smallest(idmatr,lastmatr);
    Xo = idmatr - lastmatr;
    
    % Make psi staggered
    [Rmatr,jb] = rref_exact((flipud(Xo))', symbolic);
    jb = sort(size(Xo,1) + 1 - jb);
    [Lmatr,Umatr,Pmatr] = lu((flipud(Xo(jb,:)))');
    invP = flipud(Lmatr'*Pmatr);
    P = inv(invP);
    Xo = Xo*P;
    
    % Construct Xtildeo
    newcols = Gsegment(g0tilde, (Ltilde:Rtilde) + K - N, 0:(2*Ttilde + K - N + Rtilde), 2*(Ntilde:Ttilde), symbolic);
    [Rmatr, newcols] = expand_cols_smallest([Xtildee; Ztildee], newcols);
    Gtilde = [Rmatr newcols];
    
    newcols = Gsegment(g0, (L:R) + K - N, 0:(2*Ttilde + K - N + R), 2*(N:Ttilde), symbolic);
    [Rmatr, newcols] = expand_cols_smallest([Xe; Ze], newcols);
    G = [Rmatr newcols];
    
    lastmatr = Gtilde*(G(Stilde,:))';
    idmatr = eye(Stilde(N0+K-N));
    idmatr = idmatr(:,Stilde);
    [idmatr,lastmatr]=expand_cols_smallest(idmatr,lastmatr);
    Xtildeo = idmatr - lastmatr;
    
    % Make psitilde staggered
    [Rmatr,jb] = rref_exact(double((flipud(Xtildeo))'), symbolic);
    jb = sort(size(Xtildeo,1) + 1 - jb);
    [Lmatr,Umatr,Pmatr] = lu((flipud(Xtildeo(jb,:)))');
    invP = flipud(Lmatr'*Pmatr);
    P = inv(invP);
    Xtildeo = Xtildeo*P;
    
    % Bi-orthogonalize psi and psitilde
    [Xo,Xtildeo]=expand_cols_smallest(Xo,Xtildeo);
    Y = Xo'*Xtildeo;
    
    if isortho
        P1 = inv(chol(Y));
        P2 = P1;
    else
        [Lmatr,Umatr]=lu_nopivot(Y, symbolic);
        P1 = (inv(Lmatr))';
        P2 = inv(Umatr);
    end
    Xo = Xo*P1;
    Xtildeo = Xtildeo*P2;
    
    if ~isortho & normalize % Normalize. First phi-functions, then psi-functions. Work numerically from now on.
        Xe = double(Xe); Ze = double(Ze); C = double(C);
        Xtildee = double(Xtildee); Ztildee = double(Ztildee); Ctilde = double(Ctilde);
        [phi1_phi1, norm_phi_internal_squared]=find_phi_grammian(L, R, K, Xe, Ze, N, C, g0);
        P1phi = eye(N); P2phi = eye(Ntilde);
        for k=1:min(N,Ntilde)
            P1phi(k,k) = sqrt(norm_phi_internal_squared/phi1_phi1(k,k));
            P2phi(k,k) = sqrt(phi1_phi1(k,k)/norm_phi_internal_squared);
        end
        phi1_phi1 = P1phi*phi1_phi1*P1phi;
        [Xe,Ze,C,Xtildee,Ztildee,Ctilde]=update_vars(Xe,Ze,C,Xtildee,Ztildee,Ctilde,P1phi,P2phi);
    
        Xo = double(Xo); 
        Xo(1:min(N,Ntilde),:) =           double(P2phi(1:min(N,Ntilde),1:min(N,Ntilde))*Xo(1:min(N,Ntilde),:));
        Xtildeo(1:min(N,Ntilde),:) = double(P1phi(1:min(N,Ntilde),1:min(N,Ntilde))*Xtildeo(1:min(N,Ntilde),:)); 
        [psi1_psi1, norm_psi_internal_squared]=find_psi_grammian(L, R, K, Xo(1:N,:), Xo((N+1):end,:), phi1_phi1, C, g0, g1);
        P1psi = diag(sqrt( norm_psi_internal_squared./diag(psi1_psi1) ));
        P2psi = diag(sqrt( diag(psi1_psi1)/norm_psi_internal_squared  ));
        Xo = Xo*P1psi;
        Xtildeo = Xtildeo*P2psi;
    end
    
    A_pre_inv = C( (end-N+1):end, (end-N+1):end );
    Atilde_pre_inv = Ctilde( (end-Ntilde+1):end, (end-Ntilde+1):end );
    
    [Xe,Xtildee]=expand_cols_smallest([Xe; Ze],[Xtildee; Ztildee]);
    
    % Assemble W and Wtilde
    [Xe, Xo] = expand_cols_smallest(Xe, Xo); X = [Xe Xo];
    W=     assembleW(Xe,      Xo,      g0,      L:R,           g1,      (-Rtilde):(-Ltilde), K - N, N,      N0, symbolic);
    
    [Xtildee, Xtildeo] = expand_cols_smallest(Xtildee, Xtildeo); Xtilde = [Xtildee Xtildeo];
    Wtilde=assembleW(Xtildee, Xtildeo, g0tilde, Ltilde:Rtilde, g1tilde, (-R):(-L),           K - N, Ntilde, N0, symbolic);
    
    [W, Wtilde]=expand_cols_smallest(W, Wtilde);
    % only dyadic fractions?
    % Wtilde % only dyadic fractions?
    % res = W'*Wtilde; % Should be identity matrix
end    

function [Xe,Ze,C,Xtildee,Ztildee,Ctilde]=update_vars(Xe,Ze,C,Xtildee,Ztildee,Ctilde,P1,P2)
    Xe = inv(P1)*Xe*P1;
    Ze = Ze*P1;
    C  = C*P1;
    
    Xtildee = inv(P2)*Xtildee*P2;
    Ztildee = Ztildee*P2;
    Ctilde  = Ctilde*P2;
end

function [L,R]=findsupports(g0)
     % Find L, R, from g0
    if mod(length(g0),2) == 0 
        k = length(g0)/2;
        L = -k+1; R = k;
    else
        k = (length(g0)-1)/2;
        L = -k; R = k;
    end
end

function W=assembleW(Xe, Xo, g0, suppg0, g1, suppg1, KminN, N, N0, symbolic)
    numcols = KminN + 2*max(N,N0);
    if symbolic
        W = sym(zeros(size(Xe,1),numcols));
    else
        W = zeros(size(Xe,1),numcols);
    end
    W(:,1:(KminN)) = Xo(:,1:(KminN));                       % K-N psi-functions at the beginning.
    W(:,KminN     + 2*(1:N0)) = Xo(:,(KminN+1):end);  % the remaining psi-functions
    W(:,KminN - 1 + 2*(1:N))  = Xe;                 % all phi functions
    
    if N > N0  % Add internal psi-functions in W
        % Add coordinates of \bpsi_{0,N_0},...,\bpsi_{0,N'-1}^b in phi_1^b. These come at columns (K-N+2N0+1):(K-N+2(N'-1)+1)
        insertpsi = Gsegment(g1, suppg1, 0:(KminN + 2*(N-1) + 1 + suppg1(end)), KminN + 2*(N0:(N-1)) + 1, symbolic);
        [W, insertpsi] = expand_cols_smallest(W, insertpsi);
        W( :, KminN + 2*(N0:(N-1)) + 2) = insertpsi;
    end
    if N0 > N % Add internal phi-functions in W
        % Add coordinates of \bphi_{0,N'},...,\bphi_{0,N0-1}^b in phi_1^b. These come at columns (K-N+2N'):(K-N+2(N0-1))
        insertphi = Gsegment(g0, suppg0, 0:(KminN + 2*(N0-1) + suppg0(end)), KminN + 2*(N:(N0-1)), symbolic);
        [W, insertphi] = expand_cols_smallest(W, insertphi);
        W( :, KminN + 2*(N:(N0-1)) + 1) = insertphi;
    end
end

function [Anew,Bnew]=expand_cols_smallest(A,B)
    if size(A,1) > size(B,1)
        Anew = A;
        Bnew = [ B; zeros(size(A,1)-size(B,1),size(B,2))];
    else
        Anew = [ A; zeros(size(B,1)-size(A,1),size(A,2))];
        Bnew = B;
    end
end

function val=Gsegment(g0,supp,rowrange,colrange,symbolic)
    if symbolic 
        val = sym(zeros(length(rowrange),length(colrange)));
    else
        val = zeros(length(rowrange),length(colrange));
    end
    k=1;
    for col_ind = colrange
        actualinds =  supp + col_ind;
        [intersec,i1,i2] = intersect(rowrange,actualinds);
        val(i1,k) = g0( actualinds(i2) - actualinds(1) + 1 );
        k = k+1;
    end
end

function [C,C_c]=findc(R,K,N, symbolic)
    C_c = sym(zeros(N));
    C_c(1,1) = sym(1);

    L1 = 2/sym((R+K-1)); L0 = sym((R-K)/(R+K-1));
    C_c(1,2) = L0; C_c(2,2) = L1;

    for n = 1:(N-2) % fill out column n+2. n is as in paper
        betaval = (sym(n^2)/sym(4*n^2-1))*(1-sym(n^2)/sym((R+K-1)^2)) ;
        C_c(1,n+2) = L0*C_c(1,n+1); 
        C_c(2:N,n+2) = L1*C_c(1:(N-1),n+1) + L0*C_c(2:N,n+1);
        C_c(:,n+2) = C_c(:,n+2) - betaval*C_c(:,n);
    end

    C_0 = sym(zeros(R+K-1,N));
    for k=0:(N-1)
        C_0(:,k+1) = sym((-R+1):(K-1)).^k;
    end

    C = C_0*C_c;
    if ~symbolic
        C = double(C);
        C_c = double(C_c);
    end
    % C'*C check for orthogonality
end

function x=find_x(L, R, g0)
    J = (L-R+1):(R-L-1);
    I = (L + J(1)):(R+ J(end));
    seg1 = Gsegment(g0, L:R, I, 2*J, 0);
    seg2 = Gsegment(g0, L:R, I,   J, 0);
    [V,D] = eig((2/sum(g0)^2)*seg1'*seg2); % Modified.
    [M, I] = min(abs(diag(D)-1));          % Find the eigenvalue closest to 1.
    x = V(:, I(1));
    x = x/sum(x); % Scale so that \int\phi(t)dt=1. 
end

function [phi1_phi1, norm_phi_internal_squared]=find_phi_grammian(L, R, K, Xe, Ze, N, C, g0)
    x0 = find_x(L, R, g0); x1 = (2/sum(g0)^2)*x0;
    phi1_phi2 = C'*Gsegment(x1, (L-R+1):(R-L-1), (-R+1):(K-1), K:(K + size(Ze, 1) - 1), 0);
    phi2_phi2 =    Gsegment(x1, (L-R+1):(R-L-1), 1:size(Ze, 1) , 1:size(Ze, 1), 0);
    norm_phi_internal_squared = x0(R-L);
    
    combined = Xe'*phi1_phi2*Ze;
    rhs = Ze'*phi2_phi2*Ze + combined + combined';
    ls = eye(N^2) - (2/sum(g0)^2)*kron(Xe',Xe'); % Modified
    rs = reshape(rhs, N^2, 1);
    phi1_phi1 = reshape(ls\rs, N, N);
end

function [psi1_psi1, norm_psi_internal_squared]=find_psi_grammian(L, R, K, Xo, Zo, phi1_phi1, C, g0, g1)
    x0 = find_x(L, R, g0); x1 = (2/sum(g0)^2)*x0;
    phi1_phi2 = C'*Gsegment(x1, (L-R+1):(R-L-1), (-R+1):(K-1), K:(K + size(Zo, 1) - 1), 0);
    phi2_phi2 =    Gsegment(x1, (L-R+1):(R-L-1), 1:size(Zo, 1) , 1:size(Zo, 1), 0);
    
    seg=Gsegment(x1, (L-R+1):(R-L-1), 1:length(g1), 1:length(g1), 0);
    norm_psi_internal_squared = g1*seg*g1';
    
    combined = Xo'*phi1_phi2*Zo;
    psi1_psi1 = (2/sum(g0)^2)*Xo'*phi1_phi1*Xo + Zo'*phi2_phi2*Zo + combined + combined';
end

function [L,U]=lu_nopivot(A, symbolic)
    [m,n] = size(A); s= min(m,n);
    if symbolic
        L = sym(zeros(m,s)); U = sym(zeros(s,n));
    else
        L = zeros(m,s); U = zeros(s,n);
    end
    for k=1:s
        L(:,k) = A(:,k) /A(k,k);
        U(k,:) = A(k,:);
        A = A - L(:,k)*U(k,:);
    end
    if m>n
        L = [L [zeros(n,m-n); eye(m-n)]];
    end
    if m<n
        U = [U; [zeros(n-m,m) eye(n-m)]];
    end
end

function [Q,R]=qr_exact(A, symbolic)
    if symbolic
        n=size(A,1);
        Q= sym(zeros(n)); Q(:,1)=A(:,1);
        R=sym(eye(n));
        for j=2:n
            Q(:,j) = A(:,j) - Q(:,1:(j-1))*R(1:(j-1),j);
            R(j,(j+1):n) = Q(:,j)'*A(:,(j+1):n);
            R(j,(j+1):n) = R(j,(j+1):n)/(Q(:,j)'*Q(:,j));
        end
    else
        [Q,R]=qr(A);
    end
end

function [R,jb] = rref_exact(A, symbolic)
    if symbolic 
        R = rref(A);
        [m,n]=size(A);
        jb=sym([]);
        for k=1:m
            nonzeroes=find(R(k,:));
            if ~isempty(nonzeroes)
                jb = [jb nonzeroes(1)];
            else
                break;
            end
        end
    else 
        [R,jb] = rref(A);
    end
end