clear all;
close all;
clc;

%% Parameters
Mt=70;       % Number of Tx antennas
Mr=10;       % Number of legitimate users
Mre =10;    % Number of eavesdroppers

%% M-QAM modulator parameters
iter = 2e5; % Number of monte carlo interations


%%
main_channels = single(zeros(iter, Mr, Mt));
wiretap_channels = single(zeros(iter, Mr, Mt));


for monte=1:iter

    [monte]
    rand('seed',101*monte);


    %% Channel
    H = (1/sqrt(2))*(randn(Mr,Mt)+1i*randn(Mr,Mt));  % Main channel for the legitimate receiver, Bob.
    main_channels(monte,:,:) = single(H);

    H_eve = (1/sqrt(2))*(randn(Mre,Mt)+1i*randn(Mre,Mt)); % Wiretap channel for the eavesdropper, Eve.
    wiretap_channels(monte,:,:) = single(H_eve);

end


%%
save('main_channels_for_autoencoder.mat', 'main_channels')

%%
save('wiretap_channels_for_autoencoder.mat', 'wiretap_channels')

%%
%save('main_channels_for_autoencoder.mat', 'main_channels',  '-v7.3');

%%
%save('wiretap_channels_for_autoencoder.mat', 'wiretap_channels', '-v7.3')

