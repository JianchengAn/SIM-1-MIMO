clc;
clearvars;
close all;
Thickness = 0.05; %% Thickness of TX-SIM and RX-SIM
Pt = 10^(20/10); %% Transmit power
Sigma2 = 10^(-110/10); %% Average noise power at the receiver
c = 3*10^8; %% Speed of light
f0 = 28*10^9; %% Radio frequency
lambda = c/f0; %% Wavelength
N_max = 10; %% Number of meta-atoms on each row
PL = -20*log10(4*pi/lambda)-35*log10(250); %% Pathloss in dB
pathloss = 10^(PL/10); %% Pathloss
M = 100; %% Number of meta-atoms on each layer of TX-SIM
N = 100; %% Number of meta-atoms on each layer of RX-SIM
d_element_spacing = lambda/2; %% Element spacing
S = 4; %% Number of data streams
MonteCarlo = 10; %% Number of independent experiments
Max_L = 10; %% The maximum number of metasurface layers in TX-SIM
K = 10;  %% The number of metasurface layers in RX-SIM
NMSE = zeros(MonteCarlo,1);
Capacity = zeros(MonteCarlo,1);
NMSE_average = zeros(Max_L,1);
Capacity_average = zeros(Max_L,1);
for ii = 1:Max_L
    L = ii; %% The number of metasurface layers in TX-SIM
    tic
    Derivative_transmit_phase_shift = zeros(M,L);
    Derivative_receive_phase_shift = zeros(N,K);
    W_T = zeros(M,M);
    Corr_T = zeros(M,M);
    U_R = zeros(N,N);
    C_single_stream = zeros(S,1);
    Corr_R = zeros(N,N);
    Num_initialization = (max(L,K)*10); %% Initialization times
    Error_old_set = zeros(Num_initialization,1);
    phase_transmit_set = zeros(M,L,Num_initialization);
    phase_receive_set = zeros(N,K,Num_initialization);
    d_layer_spacing_transmit = Thickness/L; %% Adjacent layer spacing in TX-SIM
    d_layer_spacing_receive = Thickness/K; %% Adjacent layer spacing in RX-SIM
    W_T_1 = zeros(M,S);
    U_R_1 = zeros(S,N);
    Temp1 = zeros(S,1);
    Temp2 = zeros(S,1);
    %% Calculate inter-layer transmission coefficient matrix W_T and channel correlation matrix Corr_T associated with TX-SIM
    for mm1 = 1:M
        m_z = ceil(mm1/N_max); %% Eq. (3)
        m_x = mod(mm1-1,N_max)+1; %% Eq. (3)
        for mm2 = 1:M
            n_z = ceil(mm2/N_max); %% Eq. (3)
            n_x = mod(mm2-1,N_max)+1; %% Eq. (3)
            d_temp  = sqrt(  (m_x-n_x)^2 +  (m_z-n_z) ^2 )*d_element_spacing; %% Eq. (1)
            d_temp2 = sqrt(d_layer_spacing_transmit^2 + d_temp^2); %% Eq. (5)
            W_T(mm2,mm1) = lambda/4/pi/d_temp2*exp(-1i*2*pi*d_temp2/lambda); %% old model
            %% lambda^2/4*(d_layer_spacing_transmit/d_temp2/d_temp2*(1/2/pi/d_temp2-1i/lambda))*exp(1i*2*pi*d_temp2/lambda); %% new model
            Corr_T(mm2,mm1) = sinc(2*d_temp/lambda); %% Eq. (14)
        end
    end
    %% Calculate inter-layer transmission coefficient matrix U_R and channel correlation matrix Corr_R associated with RX-SIM
    for nn1 = 1:N
        m_z = ceil(nn1/N_max); %% Eq. (4)
        m_x = mod(nn1-1,N_max)+1; %% Eq. (4)
        for nn2 = 1:N
            n_z = ceil(nn2/N_max); %% Eq. (4)
            n_x = mod(nn2-1,N_max)+1; %% Eq. (4)
            d_temp  = sqrt( (m_x-n_x)^2 + (m_z-n_z)^2 )*d_element_spacing; %% Eq. (2)
            d_temp2 = sqrt(d_layer_spacing_receive^2 + d_temp^2); %% Eq. (6)
            U_R(nn2,nn1) = lambda/4/pi/d_temp2*exp(-1i*2*pi*d_temp2/lambda); %% old model
            Corr_R(nn2,nn1) = sinc(2*d_temp/lambda); %% Eq. (15)
        end
    end
    %% The channel from transmitter to the first layer of TX-SIM
    for mm = 1:M
        m_z = ceil(mm/N_max);
        m_x = mod(mm-1,N_max)+1;
        for nn = 1:S
            d_transmit = sqrt(d_layer_spacing_transmit^2 + ...
                ( (m_x-(1+N_max)/2)*d_element_spacing )^2 + ...
                ( (m_z-(1+N_max)/2)*d_element_spacing - (nn-(1+S)/2)*lambda/2 )^2 ); %% Eq. (7)
            W_T_1(mm,nn) = lambda/4/pi/d_transmit*exp(-1i*2*pi*d_transmit/lambda);
        end
    end
    %% The channel from the last layer of RX-SIM to the receiver
    for mm = 1:N
        m_z = ceil(mm/N_max);
        m_x = mod(mm-1,N_max)+1;
        for nn = 1:S
            d_receive = sqrt(d_layer_spacing_receive^2 +...
                ( (m_x-(1+N_max)/2)*d_element_spacing  )^2 +...
                ( (m_z-(1+N_max)/2)*d_element_spacing - (nn-(1+S)/2)*lambda/2 )^2 ); %% Eq. (8)
            U_R_1(nn,mm) = lambda/4/pi/d_receive*exp(-1i*2*pi*d_receive/lambda);
        end
    end
    rng(1)
    for jj = 1:MonteCarlo
        G_independent = sqrt(1/2)*(randn(N,M)+1i*randn(N,M));
        G = sqrt(pathloss)*(Corr_R)^(1/2)*G_independent*(Corr_T)^(1/2); %% HMIMO channel
        [G_left, G_svd, G_right] = svd(G); %% SVD of HMIMO channel
        H_true = G_svd(1:S,1:S); %% Target channel
        H_true_vec = vec(H_true);
        Norm_H = norm(H_true_vec)^2; %% Norm of the target end-to-end channel
        h_diag = diag(H_true);
        %% Power allocation using water-filling algorithm
        if S == 1
            PA_WF = Pt;
        else
            [ PA_WF ] = WF( Pt, Sigma2, h_diag );
        end
        %% Random initialization
        for tt = 1:Num_initialization
            phase_transmit = randn(M,L)+1i*randn(M,L);
            phase_transmit = phase_transmit./abs(phase_transmit); %% TX-SIM phase shifts
            phase_receive = randn(N,K)+1i*randn(N,K);
            phase_receive = phase_receive./abs(phase_receive); %% RX-SIM phase shifts
            %% Calculate TX-SIM response
            P = diag(phase_transmit(:,1))*W_T_1;
            for l=1:L-1
                P = diag(phase_transmit(:,l+1))*W_T*P;
            end
            %% Calculate RX-SIM response
            Q = U_R_1*diag(phase_receive(:,1));
            for k = 1:K-1
                Q = Q*U_R*diag(phase_receive(:,k+1));
            end
            H_SIM = Q*G*P; %% Practical SIM-aided end-to-end channel
            H_SIM_vec = vec(H_SIM);
            Factor = (H_SIM_vec'*H_SIM_vec)\H_SIM_vec'*H_true_vec; %% Compensation factor
            Error_old_set(tt) = norm(Factor*H_SIM_vec-H_true_vec)^2/Norm_H; %% Error
            phase_transmit_set(:,:,tt) = phase_transmit;
            phase_receive_set(:,:,tt) = phase_receive;

        end
        [~,d_max] = min(Error_old_set); %% Select the best one from the random codebook
        Error_old = Error_old_set(d_max);
        phase_transmit = phase_transmit_set(:,:,d_max);
        phase_phase_transmit = angle(phase_transmit);
        phase_receive = phase_receive_set(:,:,d_max);
        phase_phase_receive = angle(phase_receive);
        %% Calculate TX-SIM response
        P = diag(phase_transmit(:,1))*W_T_1;
        for l=1:L-1
            P = diag(phase_transmit(:,l+1))*W_T*P;
        end
        %% Calculate RX-SIM response
        Q = U_R_1*diag(phase_receive(:,1));
        for k = 1:K-1
            Q = Q*U_R*diag(phase_receive(:,k+1));
        end
        H_SIM = Q*G*P; %% The end-to-end channel including TX-SIM and RX-SIM
        H_SIM_vec = vec(H_SIM);
        Factor = (H_SIM_vec'*H_SIM_vec)\H_SIM_vec'*H_true_vec; %% Compensation factor
        step = 0.1; %% learning rate
        Error_new = 10000; %% A preset value as large as possible
        while abs(Error_new-Error_old) >= Error_old * 0.001
            %% Calculate partial derivative values associated with TX-SIM phase shifts Eq. (23)
            for ll = 1:L
                for mm = 1:M
                    X_left = W_T_1;
                    for ll_left = 1:ll-1
                        X_left = W_T*diag(phase_transmit(:,ll_left))*X_left;
                    end
                    X_right = Q*G;
                    for ll_right = 1:(L-ll)
                        X_right = X_right*diag(phase_transmit(:,L+1-ll_right))*W_T;
                    end
                    for ss1 = 1:S
                        temp1 = X_right(:,mm)*X_left(mm,ss1);
                        Temp1(ss1) = 2*imag((Factor*phase_transmit(mm,ll)*temp1)'*(Factor*H_SIM(:,ss1)-H_true(:,ss1)));
                    end
                    Derivative_transmit_phase_shift(mm,ll) = sum(Temp1);
                end
            end
            %% Calculate partial derivative values associated with RX-SIM phase shifts Eq. (24)
            for kk = 1:K
                for nn = 1:N
                    Y_left = U_R_1;
                    for kk_left = 1:kk-1
                        Y_left = Y_left*diag(phase_receive(:,kk_left))*U_R;
                    end
                    Y_right = G*P;
                    for kk_right = 1:(K-kk)
                        Y_right = U_R*diag(phase_receive(:,K+1-kk_right))*Y_right;
                    end
                    for ss1 = 1:S
                        Y = Y_left(ss1,nn)*Y_right(nn,:);
                        Temp2(ss1) = 2*imag((Factor*H_SIM(ss1,:)-H_true(ss1,:))*(Factor*phase_receive(nn,kk)*Y)');
                    end
                    Derivative_receive_phase_shift(nn,kk) = sum(Temp2);
                end
            end
            Derivative_transmit_phase_shift = pi*Derivative_transmit_phase_shift/max(max(Derivative_transmit_phase_shift)); %% Eq. (27)
            phase_phase_transmit = phase_phase_transmit-step*Derivative_transmit_phase_shift; %% Update phase shifts of TX-SIM
            phase_transmit = exp(1i*phase_phase_transmit);
            Derivative_receive_phase_shift = pi*Derivative_receive_phase_shift/max(max(Derivative_receive_phase_shift)); %% Eq. (28)
            phase_phase_receive = phase_phase_receive-step*Derivative_receive_phase_shift; %% Update phase shifts of RX-SIM
            phase_receive = exp(1i*phase_phase_receive);
            step = step*0.5; %% Update learning rate
            %% Calculate TX-SIM response
            P = diag(phase_transmit(:,1))*W_T_1;
            for l=1:L-1
                P = diag(phase_transmit(:,l+1))*W_T*P;
            end
            %% Calculate RX-SIM response
            Q = U_R_1*diag(phase_receive(:,1));
            for k = 1:K-1
                Q = Q*U_R*diag(phase_receive(:,k+1));
            end
            H_SIM = Q*G*P; %% The end-to-end channel including TX-SIM and RX-SIM
            H_SIM_vec = vec(H_SIM);
            Factor = (H_SIM_vec'*H_SIM_vec)\H_SIM_vec'*H_true_vec; %% Compensation factor
            Error_old = Error_new;
            Error_new = norm(Factor*H_SIM-H_true)^2/Norm_H; %% Update residual error
        end
        NMSE(jj) = Error_new;
        for pp = 1:S
            C_single_stream(pp) = log2(1+PA_WF(pp)*abs(Factor*H_SIM(pp,pp))^2/ ...
                (Sigma2+(abs(Factor*H_SIM(pp,:)).^2*PA_WF-PA_WF(pp)*abs(Factor*H_SIM(pp,pp))^2)));
        end
        Capacity(jj) = sum(C_single_stream); %% Capacity under the current setups
    end
    NMSE_average(ii) = mean(NMSE); %% Calculate NMSE
    Capacity_average(ii) = mean(Capacity); %% Ergodic capacity
    toc
end
figure
imagesc(abs(Q*G*P))
figure
plot(NMSE_average)
NMSE_K_10 = NMSE_average;
save NMSE_K_10 NMSE_K_10
figure;
plot(Capacity_average)
Capacity_K_10 = Capacity_average;
save Capacity_K_10 Capacity_K_10