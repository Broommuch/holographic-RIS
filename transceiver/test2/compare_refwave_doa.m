function compare_refwave_doa( ...
    SNRdB, MC, ...
    theta_true, phi_true, alpha_true, ...
    Ny, Nz, dy, dz, lambda, ...
    S_pilot)

L = length(theta_true);
K = Ny * Nz;
Tp = size(S_pilot,2);

b_ref = ones(K,1);    % with reference wave
b_noref = zeros(K,1);% without reference wave

MSE_theta_ref = zeros(length(SNRdB),1);
MSE_phi_ref   = zeros(length(SNRdB),1);
MSE_theta_noref = zeros(length(SNRdB),1);
MSE_phi_noref   = zeros(length(SNRdB),1);

for isnr = 1:length(SNRdB)

    sigma2 = 10^(-SNRdB(isnr)/10);

    err_th_ref = 0; err_ph_ref = 0;
    err_th_nr  = 0; err_ph_nr  = 0;

    for mc = 1:MC

        %% ===== generate pilot measurements =====
        I_ref   = zeros(K,Tp);
        I_noref = zeros(K,Tp);

        for t = 1:Tp
            y = zeros(K,1);
            for l = 1:L
                v = steering2D(theta_true(l),phi_true(l),Ny,Nz,dy,dz,lambda);
                y = y + alpha_true(l)*v*S_pilot(l,t);
            end

            I_ref(:,t)   = abs(y + b_ref).^2   + sqrt(sigma2)*randn(K,1);
            I_noref(:,t) = abs(y + b_noref).^2 + sqrt(sigma2)*randn(K,1);
        end

        %% ===== GN initialization (same for fairness) =====
        theta0 = theta_true + 5*pi/180*randn(1,L);
        phi0   = phi_true   + 5*pi/180*randn(1,L);
        alpha0 = alpha_true .* (1+0.2*(randn(L,1)+1j*randn(L,1)));

        x0 = pack_param(theta0,phi0,alpha0);

        %% ===== with reference wave =====
        x = x0;
        for iter = 1:25
            [r,J] = residual_jacobian_doa(x,I_ref,S_pilot,b_ref,...
                                          Ny,Nz,dy,dz,lambda);
            dx = -(J.'*J)\(J.'*r);
            x = x + dx;
            if norm(dx) < 1e-5, break; end
        end
        [th,ph,~] = unpack_param(x,L);
        err_th_ref = err_th_ref + mean((th-theta_true.').^2);
        err_ph_ref = err_ph_ref + mean((ph-phi_true.').^2);

        %% ===== without reference wave =====
        x = x0;
        for iter = 1:25
            [r,J] = residual_jacobian_doa(x,I_noref,S_pilot,b_noref,...
                                          Ny,Nz,dy,dz,lambda);
%             dx = -(J.'*J)\(J.'*r);
            mu = 1e-6 * trace(J.'*J) / size(J,2);   % 自适应阻尼
            dx = -(J.'*J + mu*eye(size(J,2))) \ (J.'*r);

            x = x + dx;
            if norm(dx) < 1e-5, break; end
        end
        [th,ph,~] = unpack_param(x,L);
        err_th_nr = err_th_nr + mean((th-theta_true.').^2);
        err_ph_nr = err_ph_nr + mean((ph-phi_true.').^2);
    end

    MSE_theta_ref(isnr)   = err_th_ref / MC;
    MSE_phi_ref(isnr)     = err_ph_ref / MC;
    MSE_theta_noref(isnr) = err_th_nr  / MC;
    MSE_phi_noref(isnr)   = err_ph_nr  / MC;

    fprintf('SNR %d dB done\n', SNRdB(isnr));
end

%% ===== plot =====
figure;
semilogy(SNRdB,MSE_theta_ref,'-o','LineWidth',1.5); hold on;
semilogy(SNRdB,MSE_theta_noref,'--o','LineWidth',1.5);
semilogy(SNRdB,MSE_phi_ref,'-s','LineWidth',1.5);
semilogy(SNRdB,MSE_phi_noref,'--s','LineWidth',1.5);
grid on;

xlabel('SNR (dB)');
ylabel('MSE');
legend( ...
    'Elevation (with ref)', ...
    'Elevation (no ref)', ...
    'Azimuth (with ref)', ...
    'Azimuth (no ref)', ...
    'Location','southwest');

title('Impact of Reference Wave on DOA Estimation');

end
