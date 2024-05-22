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
MonteCarlo = 50; %% Number of independent experiments
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
        [ PA_WF ] = WF( Pt, Sigma2, h_diag );
        Capacity(jj) = sum(log2(1 + PA_WF.*h_diag.^2/Sigma2));
    end
    Capacity_average(ii) = mean(Capacity); %% Ergodic capacity
    toc
end

figure
plot(Capacity_average)
Capacity_max = Capacity_average;
save Capacity_max Capacity_max