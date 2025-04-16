function ser = calculate_ser( alg_in_real, alg_in_imag, alg_out_real, alg_out_imag, nn_out_real, nn_out_imag, main_chan_real, main_chan_imag, sym)

snr_db = 10;

alg_in = single(complex(alg_in_real, alg_in_imag));
alg_out = single(complex(alg_out_real, alg_out_imag));
nn_out = single(complex(nn_out_real, nn_out_imag));
main_chan = single(complex(main_chan_real, main_chan_imag));


batch = size(main_chan, 1);
algorithm_input_size = size(alg_in,2);
algorithm_output_size = size(alg_out,2);
nn_output_size = size(nn_out,2);

Mr = 10;
Mt = 70;
nf  = (Mr/(Mt-Mr));
M = 16; 
theta = 0.9;


SNR =10^(snr_db/10);
sigma2n=1/SNR;

wn = single(sqrt(sigma2n/2)*(randn(batch, Mr)+1i*randn(batch, Mr)));

algorithm_input_reshaped = reshape(alg_in, [batch, 1, algorithm_input_size]);
algorithm_output_reshaped = reshape(alg_out, [batch, 1, algorithm_output_size]);
nn_output_reshaped = reshape(nn_out, [batch, 1, nn_output_size]);


output_1 = sum(main_chan .* algorithm_input_reshaped, 3) + wn;
output_1_sym = qamdemod(sqrt(10*nf)*output_1,M);

%output_1 = sum(main_chan .* algorithm_input_reshaped, 3) + wn;
%output_1 = output_1.*normFactor*sqrt((10*Mr)/theta)*sqrt(theta/Mr + (1-theta)/Mt);
%output_1_sym = qamdemod(output_1,M);

output_2 = sum(main_chan .* algorithm_output_reshaped, 3) + wn;
output_2_sym = qamdemod(sqrt((10*nf)/theta)*output_2,M);

output_3 = sum(main_chan .* nn_output_reshaped, 3) + wn;
output_3_sym = qamdemod(sqrt((10*nf)/theta)*output_3,M);

ser_1 = sum((output_1_sym~=sym),2)/(Mr); % ser for algorithm input
ser_2 = sum((output_2_sym~=sym),2)/(Mr); % ser for algorithm output
ser_3 = sum((output_3_sym~=sym),2)/(Mr); % ser for neural network output

ser = cat(2, ser_1, ser_2, ser_3);

end