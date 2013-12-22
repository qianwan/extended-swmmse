function [mmse, omega, ols] = updateLPSWMmseMatrix(K, Q, M, I, N, H, U, W, S, T, maxNumCand)
    mmse = zeros(maxNumCand * M, K * I * maxNumCand * M);
    omega = zeros(K * I, maxNumCand);
    ols = zeros(K * I, K * Q);
    fMik = zeros(K * Q * M, K * Q * M);
    flag = 1;
    for ik = 1 : K * I
        rowOffset = (ik - 1) * N + 1 : ik * N;
        u = U(rowOffset, 1);
        if norm(u, 2) == 0
            continue;
        end
        Tik = T(ik, T(ik, :) ~= 0);
        Mik = zeros(length(Tik), length(Tik));
        for jm = 1 : K * I
            rowOffset = (jm - 1) * N + 1 : jm * N;
            h = H(rowOffset, Tik);
            u = U(rowOffset, 1);
            hpu = h' * u;
            Mik = Mik + W(jm) * hpu * hpu';
        end
        offset = (ik - 1) * maxNumCand * M + 1 : (ik - 1) * maxNumCand * M + length(Tik);
        mmse(1 : length(Tik), offset) = Mik;
        Sik = S(ik, S(ik, :) ~= 0);
        for s = 1 : length(Sik)
            offset = (s - 1) * M + 1 : s * M;
            omega(ik, s) = 2 * max(real(eig(Mik(offset, offset))));
        end
        rowOffset = (ik - 1) * N + 1 : ik * N;
        h = H(rowOffset, Tik);
        u = U(rowOffset, 1);
        assert(rank(Mik) == size(Mik, 1), 'rank deficient, how to obtain ols?');
        vik = W(ik) * Mik \ (h' * u);
        % if flag == 1
        %     flag = 2;
        %     fprintf(2, 'mmse rho=%f, w=%f, u=%f, v=%f\n', max(eig(Mik)), W(ik), norm(u), norm(vik));
        % end
        for s = 1 : length(Sik)
            ols(ik, Sik(s)) = norm(vik((s - 1) * M + 1 : s * M, 1));
            assert(norm(ols(ik, Sik(s))) ~= 0, 'ols should be positive');
        end
    end
    return
