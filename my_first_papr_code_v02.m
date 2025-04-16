clear all;
close all;
clc;

% papr small => secrecy cap low (artificial noise variance is low) and the error prob gap with ZF  will be
% small (not lots of power is lost for secrecy )


%% Parameters
Mt=70;       % Number of Tx antennas
Mr=10;       % Number of legitimate users
Mre=10;    % Number of eavesdroppers
M = 16;     % QAM size
theta = 0.9; % Mr/(Mr+Mre);

SNRdB = 0:2:20;

len = length(SNRdB);
SNR = zeros(1,length(SNRdB));

SERm_papr = zeros(1,len);
SERm_an = zeros(1,len);
SERm_zf = zeros(1,len);

SERm_papr_func = zeros(1,len);
SERm_an_func = zeros(1,len);
SERm_zf_func = zeros(1,len);

sec_cap_zf = zeros(1,len);
sec_cap_zf_papr = zeros(1,len);
sec_cap_zf_an = zeros(1,len);

hist_int = 1000;
y = randi(1000,1);
%% M-QAM modulator parameters
L = 4;
N = 128;
iter = 1e3; % Number of monte carlo interations
NiterVect = 1; % Number of iterations for PAPR reduction

AN_store = zeros(Mt, iter);
Vk_store = zeros(Mt, iter);

x_store1 = zeros(Mt, iter);
x_store2 = zeros(Mt, iter);
x_store_useful1 = zeros(Mt, iter);
x_store_AN1= zeros(Mt, iter);
x_store_useful2 = zeros(Mt, iter);
x_store_AN2= zeros(Mt, iter);

for kk = 1:length(SNRdB)

    SER_papr = 0;
    SER_an = 0;
    SER_zf = 0;

    SER_papr_func=0;
    SER_an_func=0;
    SER_zf_func=0;

    bob_cap1 = zeros(1,iter);
    bob_cap2 = zeros(1,iter);

    eve_cap_zf = zeros(1,iter);
    eve_cap_zf_papr = zeros(1,iter);
    eve_cap_zf_an = zeros(1,iter);

    for monte=1:iter

        [SNRdB(kk) monte]
        rand('seed',101*monte);

        %% Input data

        sent_bits = randi([0, 1],Mr*log2(M),1);
        sent_syms=  8*sent_bits(1:4:end)+4*sent_bits(2:4:end) +2*sent_bits(3:4:end)+sent_bits(4:4:end);
        qam_syms = qammod(sent_syms, M)/sqrt(10);
        sn = reshape(qam_syms,Mr,1);


        %% Channel
        H = (1/sqrt(2))*(randn(Mr,Mt)+1i*randn(Mr,Mt));  % Main channel for the legitimate receiver, Bob.

        H_eve = (1/sqrt(2))*(randn(Mre,Mt)+1i*randn(Mre,Mt)); % Wiretap channel for the eavesdropper, Eve.

        %%  MASSIVE MIMO precoding and PAPR reduction

        %% ZF precoder
        F = H'*inv(H*H');
        nf  = Mr/(Mt-Mr);

        %% Null space matrix
        V = eye(Mt) - H'*inv(H*H')*H;
        nV = Mt-Mr;

        %% Zero Forcing Precoded Signal

        x_zf_power = sqrt(1/nf)*F*sn; % with equal power allocation

        x_zf = single(sqrt(1/nf)*F*sn);
        % signal before PAPR reduction

        d= F*sn; % 10 symbols to transmit and the size of d is 100

        %% PAPR for ZF Precoded Signal
        PAPR_ZF(monte) = max(abs(x_zf.^2))/mean(abs(x_zf.^2)); % ZF precoded signal PAPR

        %% Gaussian Artificial Noise
        k = 1/sqrt(2)*(randn(Mt,1)+1i*randn(Mt,1));
        %Vk = V*k;
        %Vk_store(:,monte) = Vk;

        x_zf_an = single(sqrt(theta/nf)*d + sqrt((1-theta)/nV)*V*k); %ZF Signal plus Random AN

        %x_store1(:,monte) = x_zf_an;
        %x_store_useful1(:,monte) = sqrt(theta/nf)*d;
        %x_store_AN1(:,monte) = sqrt((1-theta)/nV)*V*k;


        %% PAPR for ZF Precoded Signal plus Gaussian AN
        PAPR_ZF_AN(monte) = max(abs(x_zf_an.^2))/mean(abs(x_zf_an.^2));


        %% PROPOSED ALGORITM
        x_zf_gd = sqrt(1/nf)*F*sn;
        AN =0;
        e = zeros(Mt,1);
        dx = zeros(Mt,1);
        e_store = 0;
        Le = 2*max(var(V));

        for iIter=1:NiterVect

            lambda = mean(abs(x_zf_gd).^2);

            lambda_used = sqrt(10^0)*sqrt(lambda);

            xclip = x_zf_gd;

            xclip(abs(x_zf_gd)>lambda_used)=lambda_used*exp(1i*angle(x_zf_gd(abs(x_zf_gd)>lambda_used))); %why not the sqrt of lambda or abs(x_gd)^2

            z = xclip - x_zf_gd;

            e =  e - (2/Le)*(V'*(V*e - z));

            p = sum(abs(V*e).*abs(z))/sum(abs(V*e).^2);

            x_zf_gd = x_zf_gd + p*V*e;

            AN = AN + p*V*e;


            e_store = e_store + p*e;

        end

        %% PAPR for Precoded Signal after GD algorithm
        PAPR_ZF_GD(monte) = max(abs(x_zf_gd).^2)/mean(abs(x_zf_gd).^2);

        %% Combined PAPR aware AN For ZF Precoded Signal Normalized
        AN =AN/std(AN);
        %AN_store(:,monte) = AN;

        %% ZF Precoded Signal plus PAPR aware AN
        x_zf_gd_papr = single(sqrt(theta/nf)*d + sqrt((1-theta)/nV)*V*AN);

        %x_store2(:,monte) = x_zf_gd_papr;
        %x_store_useful2(:,monte) = sqrt(theta/nf)*d;
        %x_store_AN2(:,monte) = sqrt((1-theta)/nV)*AN;


        %% PAPR for ZF Precoded Signal plus PAPR aware AN (Proposal 1 for ZF)
        PAPR_ZF_GD_PAPR_AN(monte) = max(abs(x_zf_gd_papr).^2)/mean(abs(x_zf_gd_papr).^2);


        %% Received Signal
        H = single(H);
        H_eve = single(H_eve);

        SNR =10^(SNRdB(kk)/10);
        sigma2n=1/SNR;

        wn = single(sqrt(sigma2n/2)*(randn(Mr,1)+1i*randn(Mr,1)));
        wn_eve = single(sqrt(sigma2n/2)*(randn(Mre,1)+1i*randn(Mre,1)));

        %% Received Signal with PAPR Reducing AN
        rec_papr = H*x_zf_gd_papr + wn;                      % Received signal for Bob

        rcvd_sym_papr = qamdemod(sqrt((10*nf)/theta)*rec_papr,M);

        %% Received Signal with Regular AN
        rec_an = H*x_zf_an + wn;                      % Received signal for Bob

        rcvd_sym_an = qamdemod(sqrt(10*nf/theta)*rec_an,M);

        %% Received Signal with only ZF Precoding
        rec_zf = H*x_zf + wn;
        rcvd_sym_zf = qamdemod(sqrt(10*nf)*rec_zf,M);


        %% Symbol Error Rate
        H_ser = reshape(H, [1, 10, 70]);
        my_ser = calculate_ser(real(x_zf_gd.'), imag(x_zf_gd.'), real(x_zf_gd_papr.'), imag(x_zf_gd_papr.'), real(x_zf_an.'), imag(x_zf_an.'), real(H_ser), imag(H_ser), sent_syms.', SNRdB(kk));

        SER_papr = SER_papr + sum(rcvd_sym_papr~=sent_syms)/(Mr);
        SER_papr_func = SER_papr_func + my_ser(2);

        SER_an = SER_an + sum(rcvd_sym_an~=sent_syms)/(Mr);
        SER_an_func = SER_an_func + my_ser(3);

        SER_zf = SER_zf + sum(rcvd_sym_zf~=sent_syms)/(Mr);
        SER_zf_func = SER_zf_func + my_ser(1);

        %% Secrecy Capacity

        B1 = H_eve * (F * F') * H_eve' ; %for ZF Precoding

        C1 = eye(Mre) + H_eve * (V * V') * H_eve' * SNR*(1-theta)/nV; %AN for ZF Precoded signal


        bob_cap1(monte) = log2(det(1 + SNR/nf*1));

        bob_cap2(monte) = log2(1 + (theta/nf)*SNR);

        eve_cap_zf(monte) = 1/Mre*real(log2(det(eye(Mre) + SNR/nf*B1)));
        eve_cap_zf_papr(monte) = 1/Mre*real(log2(det(eye(Mre) + SNR/nf*theta*H_eve*F*C1^(-1)*F'*H_eve'))) ; % - log2(det(eye(Mre) + SNR*C1)));
        eve_cap_zf_an(monte) = eve_cap_zf_papr(monte) ; %real(log2(det(eye(Mre) + SNR*B2 + SNR*C2)) - log2(det(eye(Mre) + SNR*C2)));
    end

    %% Secrecy Capacity
    cs_zf = bob_cap1 - eve_cap_zf;
    cs_zf_papr = bob_cap2 - eve_cap_zf_papr;
    cs_zf_an = bob_cap2 - eve_cap_zf_an;

    sec_cap_zf(kk) = max(mean(cs_zf),0);
    sec_cap_zf_papr(kk) = max(mean(cs_zf_papr),0);
    sec_cap_zf_an(kk) = max(mean(cs_zf_an),0);


    %% Bob with PAPR reducing noise
    SERm_papr(kk) = SER_papr/iter;
    SERm_papr_func(kk) = SER_papr_func/iter;

    %% Bob with legacy artificil noise
    SERm_an(kk) = SER_an/iter;
    SERm_an_func(kk) = SER_an_func/iter;

    %% Bob with only ZF Precoding
    SERm_zf(kk) = SER_zf/iter;
    SERm_zf_func(kk) = SER_zf_func/iter;

end

%% PAPR Plot
[h1,t1]=hist(PAPR_ZF,hist_int);
h1=h1/length(PAPR_ZF);
t11_1dB=10*log10(t1);

[h3,t3]=hist(PAPR_ZF_AN,hist_int);
h3=h3/length(PAPR_ZF_AN);
t13_1dB=10*log10(t3);

[h5,t5]=hist(PAPR_ZF_GD,hist_int);
h5=h5/length(PAPR_ZF_GD);
t15_1dB=10*log10(t5);

[h7,t7]=hist(PAPR_ZF_GD_PAPR_AN,hist_int);
h7=h7/length(PAPR_ZF_GD_PAPR_AN);
t17_1dB=10*log10(t7);

y1_1dB = 1 - cumsum(h1);
y3_1dB = 1-cumsum(h3);
y5_1dB = 1-cumsum(h5);
y7_1dB = 1-cumsum(h7);

%% PAPR Plot
figure
semilogy(t11_1dB, y1_1dB,'--rs','LineWidth',1,'MarkerSize',10,'MarkerIndices',1:100:length(y5_1dB))
hold all
semilogy(t17_1dB, y7_1dB,'-k*','LineWidth',1,'MarkerSize',10,'MarkerIndices',1:100:length(y7_1dB))
hold all
semilogy(t13_1dB, y3_1dB,'-.bo','LineWidth',1,'MarkerSize',10,'MarkerIndices',1:100:length(y3_1dB))
hold all
grid on
axis([0 12 1e-3 1])
legend(' ZF Precoded Signal','ZF Precoded Signal plus PAPR Reducing AN',...
    'ZF Precoded Signal plus Random AN','location','northeast','FontSize',10)
xlabel('PAPR(dB)')
ylabel('CCDF(PAPR)')


%% Secrecy Capacity Plot
figure
plot(SNRdB, sec_cap_zf,'-rs', 'linewidth', 1, 'MarkerSize',10);
hold on;
plot(SNRdB, sec_cap_zf_papr,'-k*', 'linewidth', 1, 'MarkerSize',10);
hold on;
plot(SNRdB, sec_cap_zf_an,'-bo', 'linewidth', 1, 'MarkerSize',10);
hold off;
grid on;
legend('ZF Precoded Signal','ZF Precoded Signal plus PAPR Reducing AN',...
    'ZF Precoded Signal plus Random AN','location','northwest','FontSize',10)
xlabel('SNR(dB)')
ylabel('Secrecy Rate (b/s/Hz)')

%% SER Plot
figure
semilogy(SNRdB, SERm_zf,'-rs', 'linewidth', 1, 'MarkerSize',10);
hold on;
semilogy(SNRdB, SERm_papr, '-k*', 'linewidth', 1, 'MarkerSize',10);
hold on
semilogy(SNRdB, SERm_an,'-bo', 'linewidth', 1, 'MarkerSize',10);
hold off;
legend('ZF Precoded Signal','ZF Precoded Signal plus PAPR Reducing AN','ZF Precoded Signal plus Random AN',...
    'Location','southwest','FontSize',10)
grid on;
xlabel('SNR(dB)')
ylabel('SER')

%% SER Plot
figure
semilogy(SNRdB, SERm_zf_func,'-rs', 'linewidth', 1, 'MarkerSize',10);
hold on;
semilogy(SNRdB, SERm_papr_func, '-k*', 'linewidth', 1, 'MarkerSize',10);
hold on
semilogy(SNRdB, SERm_an_func,'-bo', 'linewidth', 1, 'MarkerSize',10);
hold off;
legend('ZF Precoded Signal','ZF Precoded Signal plus PAPR Reducing AN','ZF Precoded Signal plus Random AN',...
    'Location','southwest','FontSize',10)
grid on;
xlabel('SNR(dB)')
ylabel('SER')