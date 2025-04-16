clc
close all

%% PAPR Plot

nn_output = complex(total_alg_in_real, total_alg_in_imag);

papr_nn= max(abs(nn_output.^2))/mean(abs(nn_output.^2));

hist_int = 10;

[h1,t1]=hist(papr_nn,hist_int);
h1=h1/length(papr_nn);
t11_1dB=10*log10(t1);

y1_1dB = 1 - cumsum(h1);

figure
semilogy(t11_1dB, y1_1dB,':rs','LineWidth',1,'MarkerSize',10,'MarkerIndices',1:100:length(y1_1dB))
grid on
axis([0 12 1e-3 1])
legend('ZF Precoded Signal','location','northeast','FontSize',10)
xlabel('PAPR(dB)')
ylabel('CCDF(PAPR)')