len = length(total_alg_in_imag);

power_in = zeros(1,len);
power_out = zeros(1,len);
power_nn = zeros(1,len);

tic
for iter=1:len

    alg_in = total_alg_in_real(iter,:) + 1i*total_alg_in_imag(iter,:);
    power_in(iter) =  sum(abs(alg_in).^2);

    alg_out = total_alg_out_real(iter,:) + 1i*total_alg_out_imag(iter,:);
    power_out(iter) =  sum(abs(alg_out).^2);

    nn_out = total_nn_out_real(iter,:) + 1i*total_nn_out_imag(iter,:);
    power_nn(iter) =  sum(abs(nn_out).^2);

end
toc