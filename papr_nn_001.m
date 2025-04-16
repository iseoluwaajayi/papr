clear all;
close all;
clc;

%% Parameters
Mt=70;       % Number of Tx antennas
Mr=10;       % Number of legitimate users
Mre =10;    % Number of eavesdroppers
M = 16;     % QAM size
theta = 0.9;

SNRdB = 10;

SNR = zeros(1,length(SNRdB));
hist_int = 1000;

%% M-QAM modulator parameters
iter = 1e5; % Number of monte carlo interations
NiterVect = 40; % Number of iterations for PAPR reduction

%%
symbols_store = single(zeros(iter, Mr));
nf_store = single(zeros(iter,1));
main_channels = single(zeros(iter, Mr, Mt));
wiretap_channels = single(zeros(iter, Mr, Mt));
algorithm_input = single(zeros(iter, Mt));
algorithm_output = single(zeros(iter, Mt));

for kk = 1:length(SNRdB)

    for monte=1:iter

        [SNRdB(kk) monte]
        rand('seed',101*monte);

        %% Input data
        sent_bits = randi([0, 1],Mr*log2(M),1);
        sent_syms=  8*sent_bits(1:4:end)+4*sent_bits(2:4:end) +2*sent_bits(3:4:end)+sent_bits(4:4:end);
        qam_syms = qammod(sent_syms, M)/sqrt(10);
        sn = reshape(qam_syms,Mr,1);

        symbols_store(monte,:) = single(sent_syms);


        %% Channel
        H = (1/sqrt(2))*(randn(Mr,Mt)+1i*randn(Mr,Mt));  % Main channel for the legitimate receiver, Bob.
        main_channels(monte,:,:) = single(H);

        H_eve = (1/sqrt(2))*(randn(Mre,Mt)+1i*randn(Mre,Mt)); % Wiretap channel for the eavesdropper, Eve.
        wiretap_channels(monte,:,:) = single(H_eve);

        %% ZF precoder
        F = H'*inv(H*H');
        nf  = Mr/(Mt-Mr);


        %% Null space matrix
        V = eye(Mt) - H'*inv(H*H')*H;
        nV = Mt - Mr;

        %% Zero Forcing Precoded Signal

        x_zf = single(sqrt(1/nf)*F*sn); % signal before PAPR reduction

        d= F*sn; % 10 symbols to transmit and the size of d is 100

        %% PAPR for ZF Precoded Signal
        PAPR_ZF(monte) = max(abs(x_zf.^2))/mean(abs(x_zf.^2)); % ZF precoded signal PAPR

        %% Gaussian Artificial Noise
        k = 1/sqrt(2)*(randn(Mt,1)+1i*randn(Mt,1));
        
        Vk = V*k/norm(V*k);

        x_zf_an = single(sqrt(theta/nf)*d + sqrt((1-theta)/nV)*V*k); %ZF Signal plus Random AN


        %% PAPR for ZF Precoded Signal plus Gaussian AN
        PAPR_ZF_AN(monte) = max(abs(x_zf_an.^2))/mean(abs(x_zf_an.^2));


        %% PROPOSED ALGORITM
        x_zf_gd = sqrt(1/nf)*F*sn;
        AN =0;
        e = zeros(Mt,1);
        dx = zeros(Mt,1);
        e_store = 0;
        Le = 2*max(var(V));

        algorithm_input(monte,:) = single(x_zf_gd);

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

        %% Combined PAPR aware AN For ZF Precoded Signal Normalized
        AN =AN/std(AN);

        %% ZF Precoded Signal plus PAPR aware AN
        x_zf_gd_papr = sqrt(theta/nf)*d + sqrt((1-theta)/nV)*V*AN;

        algorithm_output(monte,:) = single(x_zf_gd_papr);

    end
end

%%
save('symbols_store_single.mat', 'symbols_store')
save('main_channels_single.mat', 'main_channels')
save('wiretap_channels_single.mat', 'wiretap_channels')
save('algorithm_input_single.mat', 'algorithm_input')
save('algorithm_output_single.mat', 'algorithm_output')

%%
%save('main_channels.mat', 'main_channels',  '-v7.3');

%%
%save('wiretap_channels.mat', 'wiretap_channels', '-v7.3')

