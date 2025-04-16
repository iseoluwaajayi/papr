%close all;
clc;

%%
SNRdB = 0:2:20;
len = length(SNRdB);

SERm_alg_in = zeros(1,len);
SERm_alg_out = zeros(1,len);
SERm_nn_out = zeros(1,len);

iter = 1e5;

for kk = 1:length(SNRdB)

    SER_alg_in = 0;
    SER_alg_out = 0;
    SER_nn_out = 0;

    snr_db = SNRdB(kk);

    for monte=1:iter

        [SNRdB(kk) monte]
        rand('seed',101*monte);

        alg_in_real = real(algorithm_input(monte,:));
        alg_in_imag = imag(algorithm_input(monte,:));
        alg_out_real = real(algorithm_output(monte,:));
        alg_out_imag = imag(algorithm_output(monte,:));
        nn_out_real = real(algorithm_output(monte,:));
        nn_out_imag = imag(algorithm_output(monte,:));
        main_chan_real = real(main_channels(monte,:,:));
        main_chan_imag = imag(main_channels(monte,:,:));
        sym = symbols_store(monte,:);

        ser = calculate_ser( alg_in_real, alg_in_imag, alg_out_real, alg_out_imag, nn_out_real, nn_out_imag, main_chan_real, main_chan_imag, sym, snr_db);

        %% Symbol Error Rate

        SER_alg_in_init =  ser(1);
        SER_alg_in = SER_alg_in + SER_alg_in_init;

        SER_alg_out_init = ser(2);
        SER_alg_out = SER_alg_out + SER_alg_out_init;

        SER_alg_nn_init = ser(3);
        SER_nn_out = SER_nn_out + SER_alg_nn_init;

    end

    %% Bob with only ZF Precoding
    SERm_alg_in(kk) = SER_alg_in/iter;

    %% Bob with PAPR reducing noise
    SERm_alg_out(kk) = SER_alg_out/iter;

    %% Bob with PAPR reducing noise
    SERm_nn_out(kk) = SER_nn_out/iter;

end

%% SER Plot

figure
semilogy(SNRdB, SERm_alg_in,'-rs', 'linewidth', 1, 'MarkerSize',10);
hold on;
semilogy(SNRdB, SERm_alg_out, '-k*', 'linewidth', 1, 'MarkerSize',10);
hold on
semilogy(SNRdB, SERm_nn_out,'-bo', 'linewidth', 1, 'MarkerSize',10);
hold off;
legend('SER - Alg Input','SER - Alg Output','SER - NN Output','Location','southwest','FontSize',10)
grid on;
axis([0 12 1e-5 1])
xlabel('SNR(dB)')
ylabel('SER')