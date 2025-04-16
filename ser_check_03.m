close all;
clc;

%%
SNRdB = 0:2:20;
len = length(SNRdB);

SERm_alg_in = zeros(1,len);
SERm_alg_out = zeros(1,len);
SERm_nn_out = zeros(1,len);

SERm_alg_in_test = zeros(1,len);
SERm_alg_out_test  = zeros(1,len);
SERm_nn_out_test  = zeros(1,len);

iter = 1e4;

for kk = 1:length(SNRdB)

    SER_alg_in = 0;
    SER_alg_out = 0;
    SER_nn_out = 0;

    SER_alg_in_test = 0;
    SER_alg_out_test = 0;
    SER_nn_out_test = 0;

    snr_db = SNRdB(kk);

    for monte=1:iter

        [SNRdB(kk) monte]
        rand('seed',101*monte);

        alg_in_real = total_alg_in_real(monte,:);
        alg_in_imag = total_alg_in_imag(monte,:);
        alg_out_real = total_alg_out_real(monte,:);
        alg_out_imag = total_alg_out_imag(monte,:);
        nn_out_real = total_nn_out_real(monte,:);
        nn_out_imag = total_nn_out_imag(monte,:);
        main_chan_real = total_main_channels_real(monte,:);
        main_chan_imag = total_main_channels_imag(monte,:);
        sym = total_symbols(monte,:);

        main_chan_real_test = total_nn_out_real(monte,:);
        main_chan_imag_test = total_nn_out_imag(monte,:);


        ser = calculate_ser( alg_in_real, alg_in_imag, alg_out_real, alg_out_imag, nn_out_real, nn_out_imag, main_chan_real, main_chan_imag, sym, snr_db);

        ser_test = calculate_ser( alg_in_real, alg_in_imag, alg_out_real, alg_out_imag, nn_out_real, nn_out_imag, main_chan_real_test, main_chan_imag_test, sym, snr_db);


        %% Symbol Error Rate

    
        SER_alg_in_init =  ser(1);
        SER_alg_in = SER_alg_in + SER_alg_in_init;

        SER_alg_out_init = ser(2);
        SER_alg_out = SER_alg_out + SER_alg_out_init;

        SER_alg_nn_init = ser(3);
        SER_nn_out = SER_nn_out + SER_alg_nn_init;

        SER_alg_in_init_test =  ser_test(1);
        SER_alg_in_test  = SER_alg_in_test  + SER_alg_in_init_test;

        SER_alg_out_init_test  = ser_test(2);
        SER_alg_out_test  = SER_alg_out_test  + SER_alg_out_init_test ;

        SER_alg_nn_init_test  = ser_test(3);
        SER_nn_out_test  = SER_nn_out_test  + SER_alg_nn_init_test ;

    end

    %% Bob with only ZF Precoding
    SERm_alg_in(kk) = SER_alg_in/iter;
    SERm_alg_in_test(kk) = SER_alg_in_test/iter;

    %% Bob with PAPR reducing noise
    SERm_alg_out(kk) = SER_alg_out/iter;
    SERm_alg_out_test(kk) = SER_alg_out_test/iter;

    %% Bob with PAPR reducing noise
    SERm_nn_out(kk) = SER_nn_out/iter;
    SERm_nn_out_test(kk) = SER_nn_out_test/iter;

end


%% SER Plot
figure
semilogy(SNRdB, SERm_alg_in,':rs', 'linewidth', 1, 'MarkerSize',10);
hold on;
semilogy(SNRdB, SERm_alg_out, '-k*', 'linewidth', 1, 'MarkerSize',10);
hold on;
semilogy(SNRdB, SERm_nn_out,'--md', 'linewidth', 1, 'MarkerSize',10);
hold off;
legend('SER - Alg Input','SER - Alg Output','SER - NN Output','Location','southwest','FontSize',10)
grid on;
axis([0 12 1e-5 1])
xlabel('SNR(dB)')
ylabel('SER')

%% SER Plot
figure
semilogy(SNRdB, SERm_alg_in_test ,':rs', 'linewidth', 1, 'MarkerSize',10);
hold on;
semilogy(SNRdB, SERm_alg_out_test , '-k*', 'linewidth', 1, 'MarkerSize',10);
hold on;
semilogy(SNRdB, SERm_nn_out_test ,'--md', 'linewidth', 1, 'MarkerSize',10);
hold off;
legend('SER - Alg Input','SER - Alg Output','SER - NN Output','Location','southwest','FontSize',10)
grid on;
axis([0 20 1e-5 1])
xlabel('SNR(dB)')
ylabel('SER - test')