function [W,Wtilde,A_pre_inv,Atilde_pre_inv, N0]=bw_compute_left(g0, g1, N, K, g0tilde, g1tilde, Ntilde, Ktilde)
    symbolic = 0;
    normalize = 1;
    isortho = (isequal(g0, g0tilde) & isequal(g1, g1tilde));  
    A_pre_inv = eye(N);
    Atilde_pre_inv = eye(Ntilde);
    Nprime = max(N,Ntilde);
    [L,R]=findsupports(g0);
    [Ltilde,Rtilde]=findsupports(g0tilde);
    N0=max( [ceil((Nprime+N-K-1+R)/2) ceil((Nprime+N-K-1+Rtilde)/2) R Rtilde]); % Equation (5.1)
    
    S = [1:(K-N+R) (K-N+R+2):2:(K-N-R+2*N0)]; S = S'; % Defined in Lemma 5.2
    Stilde = [1:(Ktilde-Ntilde+Rtilde) (Ktilde-Ntilde+Rtilde+2):2:(Ktilde-Ntilde-Rtilde+2*N0)]; Stilde = Stilde';
    
    T      = N0 + floor(-(R+Ltilde+1)/2); % Equation (5.3)
    Ttilde = N0 + floor(-(Rtilde+L+1)/2);
    T = max(T, Nprime); % TODO: add to the paper.
    Ttilde = max(Ttilde, Nprime);

    % Find initial boundary functions
    G_I_X = Gsegment(g0,L:R,(-R+1):(K-1),(-2*R+2):2:(2*K-2), symbolic); % Defined in Lemma 4.1
    G_I_Z = Gsegment(g0,L:R,K:(2*K+R-2),(-2*R+2):2:(2*K-2), symbolic);
    [C,C_c]=findc(R,K,N);
    if ~symbolic
        C = double(C);
    end
    Xe = inv(C'*C)*C'* G_I_X *C; % Equation (4.2). Replace inv(C'*C)*C' with transpose if matrix is orthogonal.
    Ze = G_I_Z*C;
    
    G_I_X = Gsegment(g0tilde,Ltilde:Rtilde,(-Rtilde+1):(Ktilde-1),(-2*Rtilde+2):2:(2*Ktilde-2), symbolic);
    G_I_Z = Gsegment(g0tilde,Ltilde:Rtilde,Ktilde:(2*Ktilde+Rtilde-2),(-2*Rtilde+2):2:(2*Ktilde-2), symbolic);
    [Ctilde,C_ctilde]=findc(Rtilde,Ktilde,Ntilde);
    if ~symbolic
        Ctilde = double(Ctilde);
    end
    Xtildee = inv(Ctilde'*Ctilde)*Ctilde'*G_I_X *Ctilde;
    Ztildee = G_I_Z*Ctilde;
    
    % Make phi staggered
    [Qmatr,Rmatr] = qr((flipud(C((R+K-N):end,:)))');
    P = fliplr(Qmatr);
    invP = flipud(Qmatr');
    
    
    %[Lmatr,Umatr] = lu((flipud(C((R+K-N):end,:)))');% old
    %invP = flipud(Lmatr');
    %P = inv(invP);
    Xe = inv(P)*Xe*P;
    Ze = Ze*P;
    if N == Ntilde
        A_pre_inv = C((K+R-N):end,:)*P;
    end
    C=C*P;
    
    % Make phitilde staggered
    [Qmatr,Rmatr] = qr((flipud(Ctilde((Rtilde+Ktilde-Ntilde):end,:)))');
    P = fliplr(Qmatr);
    invP = flipud(Qmatr');
    
    %[Lmatr,Umatr] = lu((flipud(Ctilde((Rtilde+Ktilde-Ntilde):end,:)))'); % old code
    %invP = flipud(Lmatr');
    %P = inv(invP);
    Xtildee = inv(P)*Xtildee*P;
    Ztildee = Ztildee*P;
    if N == Ntilde
        Atilde_pre_inv = Ctilde((Ktilde+Rtilde-Ntilde):end,:)*P;
    end
    Ctilde=Ctilde*P;
    
    % Address when Ntilde and N are different
    if Ntilde > N    
        G_I_Xp = Gsegment(g0,L:R,(2*K+L):(Ktilde-1),(2*K):2:(2*Ktilde-2), symbolic); % Defined in Lemma 4.2
        Z_I_Xp = Gsegment(g0,L:R,Ktilde:(2*Ktilde+R-2),(2*K):2:(2*Ktilde-2), symbolic);
        
        Xe = [Xe; Ze(1:(Ntilde-N),:)];
        Xe = [Xe [zeros(K+L+N,Ntilde-N); G_I_Xp]]; % Equation (4.7)
        [Ze,Z_I_Xp] = expand_cols_smallest(Ze((Ntilde-N+1):end,:), Z_I_Xp);
        Ze = [Ze Z_I_Xp]; % Equation (4.8)
    elseif N > Ntilde
        G_I_Xp = Gsegment(g0tilde,Ltilde:Rtilde,(2*Ktilde+Ltilde):(K-1),(2*Ktilde):2:(2*K-2), symbolic);
        Z_I_Xp = Gsegment(g0tilde,Ltilde:Rtilde,K:(2*K+Rtilde-2),(2*Ktilde):2:(2*K-2), symbolic);
        
        Xtildee = [Xtildee; Ztildee(1:(N-Ntilde),:)];
        Xtildee = [Xtildee [zeros(Ktilde+Ltilde+Ntilde,N-Ntilde); G_I_Xp]];
        [Ztildee,Z_I_Xp] = expand_cols_smallest(Ztildee((N-Ntilde+1):end,:), Z_I_Xp);
        Ztildee = [Ztildee Z_I_Xp];
    end
    
    % Bi-orthogonalize phi and phitilde
    ls = eye(Nprime^2) - kron(Xe',Xtildee');
    
    [Ze,Ztildee]=expand_cols_smallest(Ze,Ztildee);
    rs = reshape(Ztildee'*Ze, [Nprime^2, 1]);
    Y  = reshape(ls\rs, [Nprime, Nprime]); % Solve Equation (4.12)
    Y=Y';
    
    if isortho
        P1 = inv(chol(Y));
        P2 = P1;
    else
        [Lmatr,Umatr]= lu_nopivot(Y);
        P1 = (inv(Lmatr))';
        P2 = inv(Umatr);
    end    
    Xe = inv(P1)*Xe*P1;
    Ze = Ze*P1;
    Xtildee = inv(P2)*Xtildee*P2;
    Ztildee = Ztildee*P2;
    C=C*P1; 
    Ctilde=Ctilde*P2;
    % Y=P1'*Y*P2
    
    
    if N == Ntilde
        A_pre_inv = A_pre_inv*P1; A_pre_inv = double(A_pre_inv);
        Atilde_pre_inv = Atilde_pre_inv*P2; Atilde_pre_inv = double(Atilde_pre_inv);
    end
    
    if ~isortho & normalize
        [phi1_phi1, phi1_phi2, phi2_phi2]=find_phi_grammian(L, R, K, Xe, Ze, Nprime, C, g0);
        P1 = diag(1./sqrt(diag(phi1_phi1)));
        P2 = diag(sqrt(diag(phi1_phi1)));
        Xe = inv(P1)*Xe*P1;
        Ze = Ze*P1;
        phi1_phi1 = P1*phi1_phi1*P1; phi1_phi2 = P1*phi1_phi2; % Fixed this
        Xtildee = inv(P2)*Xtildee*P2;
        Ztildee = Ztildee*P2;
        C=C*P1; 
        Ctilde=Ctilde*P2;
        % Y=P1'*Y*P2
    
        if N == Ntilde
            A_pre_inv = A_pre_inv*P1; A_pre_inv = double(A_pre_inv);
            Atilde_pre_inv = Atilde_pre_inv*P2; Atilde_pre_inv = double(Atilde_pre_inv);
        end
    end
    
    % psi-funksjoner
    
    % Construct Xo
    newcols = Gsegment(g0, L:R, 0:(2*T + K - N + R), K - N + 2*(Nprime:T), symbolic);
    [Rmatr, newcols] = expand_cols_smallest([Xe;Ze], newcols);
    G = [Rmatr newcols];
    
    newcols = Gsegment(g0tilde, Ltilde:Rtilde, 0:(2*T + K - N + Rtilde),K - N + 2*(Nprime:T), symbolic);
    [Rmatr, newcols] = expand_cols_smallest([Xtildee;Ztildee], newcols);
    Gtilde = [Rmatr newcols];
    
    lastmatr = G*(Gtilde(S,:))';
    idmatr = eye(S(N0+K-N));
    idmatr = idmatr(:,S);
    [idmatr,lastmatr]=expand_cols_smallest(idmatr,lastmatr);
    Xo = idmatr - lastmatr;
    
    % Make psi staggered
    [Rmatr,jb] = rref(double((flipud(Xo))'));
    jb = sort(size(Xo,1) + 1 - jb);
    [Lmatr,Umatr,Pmatr] = lu((flipud(Xo(jb,:)))');
    invP = flipud(Lmatr'*Pmatr);
    P = inv(invP);
    Xo = Xo*P;
    
    % Construct Xtildeo
    newcols = Gsegment(g0tilde, Ltilde:Rtilde, 0:(2*Ttilde + K - N + Rtilde), K - N + 2*(Nprime:Ttilde), symbolic);
    [Rmatr, newcols] = expand_cols_smallest([Xtildee;Ztildee], newcols);
    Gtilde = [Rmatr newcols];
    
    newcols = Gsegment(g0, L:R, 0:(2*Ttilde + K - N + R), K - N + 2*(Nprime:Ttilde), symbolic);
    [Rmatr, newcols] = expand_cols_smallest([Xe;Ze], newcols);
    G = [Rmatr newcols];
    
    lastmatr = Gtilde*(G(Stilde,:))';
    idmatr = eye(Stilde(N0+K-N));
    idmatr = idmatr(:,Stilde);
    [idmatr,lastmatr]=expand_cols_smallest(idmatr,lastmatr);
    Xtildeo = idmatr - lastmatr;
    
    % Make psitilde staggered
    [Rmatr,jb] = rref(double((flipud(Xtildeo))'));
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
        [Lmatr,Umatr]=lu_nopivot(Y);
        P1 = (inv(Lmatr))';
        P2 = inv(Umatr);
    end
    Xo = Xo*P1;
    Xtildeo = Xtildeo*P2;
        
    if ~isortho & normalize % Normalize. Same principle as for the phi functions
        Ypsi=find_psi_grammian(Xo(1:Nprime,:), Xo((Nprime+1):end,:), phi1_phi1, phi1_phi2, phi2_phi2);
        P1 = diag(1./sqrt(diag(Ypsi)));
        P2 = diag(sqrt(diag(Ypsi)));
        Xo = Xo*P1;
        Xtildeo = Xtildeo*P2;
        % P1'*Ypsi*P1
    end
    
    % Assemble W and Wtilde
    Xe = [Xe; Ze]; Xtildee = [Xtildee; Ztildee];
    W=     assembleW(Xe,      Xo,      g0, L:R,           g1, (-Rtilde):(-Ltilde), K - N, Nprime, N0, symbolic);
    Wtilde=assembleW(Xtildee, Xtildeo, g0tilde, Ltilde:Rtilde, g1tilde, (-R):(-L),           K - N, Nprime, N0, symbolic);
    [W, Wtilde]=expand_cols_smallest(W, Wtilde);
    % only dyadic fractions?
    % Wtilde % only dyadic fractions?
    %res = W'*Wtilde; % Should be identity matrix
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

function W=assembleW(Xe, Xo, g0, suppg0, g1, suppg1, KminN, Nprime, N0, symbolic)
    [Xe, Xo] = expand_cols_smallest(Xe, Xo);
    numcols = KminN + 2*max(Nprime,N0);
    if symbolic
        W = sym(zeros(size(Xe,1),numcols));
    else
        W = zeros(size(Xe,1),numcols);
    end
    W(:,1:(KminN)) = Xo(:,1:(KminN));                       % K-N psi-functions at the beginning.
    W(:,KminN     + 2*(1:N0))     = Xo(:,(KminN+1):end);  % the remaining psi-functions
    W(:,KminN - 1 + 2*(1:Nprime)) = Xe;                 % all phi functions
    
    if Nprime > N0  % Add internal psi-functions in W
        % Add coordinates of \bpsi_{0,N_0},...,\bpsi_{0,N'-1}^b in phi_1^b. These come at columns (K-N+2N0+1):(K-N+2(N'-1)+1)
        insertpsi = Gsegment(g1, suppg1, 0:(KminN + 2*(Nprime-1) + 1 + suppg1(end)), KminN + 2*(N0:(Nprime-1)) + 1, symbolic);
        [W, insertpsi] = expand_cols_smallest(W, insertpsi);
        W( :, KminN + 2*(N0:(Nprime-1)) + 2) = insertpsi;
    end
    if N0 > Nprime % Add internal phi-functions in W
        % Add coordinates of \bphi_{0,N'},...,\bphi_{0,N0-1}^b in phi_1^b. These come at columns (K-N+2N'):(K-N+2(N0-1))
        insertphi = Gsegment(g0, suppg0, 0:(KminN + 2*N0-2 + suppg0(end)), KminN + 2*(Nprime:(N0-1)), symbolic);
        [W, insertphi] = expand_cols_smallest(W, insertphi);
        W( :, KminN + 2*(Nprime:(N0-1)) + 1) = insertphi;
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

function [C,C_c]=findc(R,K,N)
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
    C = double(C);
    C_c = double(C_c);
    % C'*C check for orthogonality
end

function g=find_g(L, R, g0)
    colinds = 2*((L-R+1):(R-L-1));
    rowinds = (L + colinds(1)):(R+ colinds(end));
    seg1 = Gsegment(g0, L:R, rowinds, colinds, 0);
    colinds = (L-R+1):(R-L-1);
    seg2 = Gsegment(g0, L:R, rowinds, colinds, 0);
    [V,D] = eig(seg1'*seg2);
    [M, I] = min(abs(diag(D)-1));
    g = V(:, I(1));
    g = g/sum(g);
end

function [Y, phi1_phi2, phi2_phi2]=find_phi_grammian(L, R, K, Xe, Ze, Nprime, C, g0)
    g = find_g(L, R, g0);
    matr = Gsegment(g, K + ((L-R+1):(R-L-1)), (-R+1):(2*K+R-2), 0:(K+R-2), 0);
    phi1_phi2 = C'*matr(1:(K+R-1),:);
    phi2_phi2 = matr((K+R):end,:);
    
    rhs = (Ze(1:size(phi2_phi2,1),:))'*phi2_phi2*(Ze(1:size(phi2_phi2,2),:)) + Xe'*phi1_phi2*(Ze(1:size(phi1_phi2,2),:)) + (Xe'*phi1_phi2*(Ze(1:size(phi1_phi2,2),:)))';
    ls = eye(Nprime^2) - kron(Xe',Xe');
    rs = reshape(rhs, [Nprime^2, 1]);
    Y  = reshape(ls\rs, [Nprime, Nprime]);
end

function Ypsi=find_psi_grammian(Xo, Zo, phi1_phi1, phi1_phi2, phi2_phi2)
    Ypsi = (Xo(1:size(phi1_phi1,1),:))'*phi1_phi1*Xo(1:size(phi1_phi1,2),:) + (Zo(1:size(phi2_phi2,1),:))'*phi2_phi2*(Zo(1:size(phi2_phi2,2),:)) + Xo'*phi1_phi2*(Zo(1:size(phi1_phi2,2),:)) + (Xo'*phi1_phi2*(Zo(1:size(phi1_phi2,2),:)))';
end

function [L,U]=lu_nopivot(A)
    n = size(A,1);
    L = zeros(n); U = zeros(n);
    for k=1:n
        L(:,k) = A(:,k) /A(k,k);
        U(k,:) = A(k,:);
        A = A - L(:,k)*U(k,:);
    end
end

function A=luimpl_banded(A)
    d = (size(A,1)-1)/2;
    n = size(A,2);
    for k=1:(n-1)
        % s = min(k + 1 + d, n);
        A((d+2):end,k) = A((d+2):end,k) /A(d+1,k);
        A((d+2):end,(k+1):end) = A((d+2):end,(k+1):end) - A((d+2):end,k)*A(d+1,(k+1):end);
    end
end