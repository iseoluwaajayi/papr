function ser = calculate_ser_2( nn_out_real, nn_out_imag, index_selection)

index_selection = uint32(index_selection);

algorithm_input_ = evalin('base', 'algorithm_input');
algorithm_output_ = evalin('base', 'algorithm_output');
main_channels_ = evalin('base', 'main_channels');
symbols_store_ = evalin('base', 'symbols_store');

alg_in = algorithm_input_(index_selection,:);
alg_out = algorithm_output_(index_selection,:);
main_chan = main_channels_(index_selection,:,:);
sym = symbols_store_(index_selection,:);
nn_out = complex(nn_out_real, nn_out_imag);


batch = size(main_chan, 1);
algorithm_input_size = size(alg_in,2);
algorithm_output_size = size(alg_out,2);
nn_output_size = size(nn_out,2);
Mr = 10;
Mt = 70;
M = 16; 
theta = 0.9;
normFactor = (sqrt(sum(abs(alg_in.').^2)))';

SNR =10^(20/10);
sigma2n=1/SNR;

wn = single(sqrt(sigma2n/2)*(randn(batch, Mr)+1i*randn(batch, Mr)));

algorithm_input_reshaped = reshape(alg_in, [batch, 1, algorithm_input_size]);
algorithm_output_reshaped = reshape(alg_out, [batch, 1, algorithm_output_size]);
nn_output_reshaped = reshape(nn_out, [batch, 1, nn_output_size]);

output_1 = sum(main_chan .* algorithm_input_reshaped, 3) + wn;
output_1 = output_1.*normFactor*sqrt((10*Mr)/theta)*sqrt(theta/Mr + (1-theta)/Mt);
output_1_sym = qamdemod(output_1,M);

output_2 = sum(main_chan .* algorithm_output_reshaped, 3) + wn;
output_2 = output_2.*normFactor*sqrt((10*Mr)/theta)*sqrt(theta/Mr + (1-theta)/Mt);
output_2_sym = qamdemod(output_2,M);

output_3 = sum(main_chan .* nn_output_reshaped, 3) + wn;
output_3 = output_3.*normFactor*sqrt((10*Mr)/theta)*sqrt(theta/Mr + (1-theta)/Mt);
output_3_sym = qamdemod(output_3,M);

ser_1 = sum((output_1_sym~=sym),2)/(Mr); % ser for algorithm input
ser_2 = sum((output_2_sym~=sym),2)/(Mr); % ser for algorithm output
ser_3 = sum((output_3_sym~=sym),2)/(Mr); % ser for neural network output

ser = cat(2, ser_1, ser_2, ser_3);

end