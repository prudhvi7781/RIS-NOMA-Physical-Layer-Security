clc
clear all
close all

snr_dB = 0:5:50;    % Total power in dB
c = 2.2;            % Pathloss component
No = 1;             % Noise power
iterations = 400;    % Number of iterations for Monte Carlo simulation
N = 30;            % Number of channel realizations

% Distances
d_sr = 50;      % Distance from Source to RIS
d_ru1 = 10;          % Distance from RIS to Destination
d_ru2 = 20;    % Distance from Source to E1
d_su1=70;
d_su2=15;
d_se1 = 30;      % Distance from RIS to E1
d_se2 = 25;          % Distance from Jammer to E1
d_se3 = 15;      % Distance from Jammer to Destination
d_ju1=70;      % Distance from Source to Destination
d_ju2=80;
d_je1=2;
d_je2=3;
d_je3=4;
a2_values = 0.76; 
a1_values = 0.24;

for a_index = 1:length(a1_values)
    a1 = a1_values(a_index);
    a2 = a2_values(a_index);
    cap_nj1 = zeros(iterations, length(snr_dB));
    cap_nj2 = zeros(iterations, length(snr_dB));

    cap_uj1 = zeros(iterations, length(snr_dB));
    cap_uj2 = zeros(iterations, length(snr_dB));

    cap_kj1 = zeros(iterations, length(snr_dB));
    cap_kj2 = zeros(iterations, length(snr_dB));

    for jj = 1:iterations
         disp(['Iteration: ', num2str(jj)]);

         for ii = 1:length(snr_dB)
             Ps = 10^(snr_dB(ii) / 10);  % Power conversion from dB to linear scale
             Pt = 0.9 * Ps;              % Power for users
             Pj = 0.1* Ps;              % Power for jamming

             %channel coefficients
             hsr = sqrt(d_sr^(-c)) * (randn(N, 1) + 1i * randn(N, 1));
             hru1 = sqrt(d_ru1^(-c)) * (randn(N, 1) + 1i * randn(N, 1));
             hru2 = sqrt(d_ru2^(-c)) * (randn(N, 1) + 1i * randn(N, 1));
             hsu1 = sqrt(d_su1^(-c)) * (randn(N, 1) + 1i * randn(N, 1));
             hsu2 = sqrt(d_su2^(-c)) * (randn(N, 1) + 1i * randn(N, 1));
             hse1 = sqrt(d_se1^(-c)) * (randn(N, 1) + 1i * randn(N, 1));
             hse2 = sqrt(d_se2^(-c)) * (randn(N, 1) + 1i * randn(N, 1));
             hse3 = sqrt(d_se3^(-c)) * (randn(N, 1) + 1i * randn(N, 1));         
             hju1 = sqrt(d_ju1^(-c)) * (randn(N, 1) + 1i * randn(N, 1));         
             hju2 = sqrt(d_ju2^(-c)) * (randn(N, 1) + 1i * randn(N, 1));
             hje1 = sqrt(d_je1^(-c)) * (randn(N, 1) + 1i * randn(N, 1));
             hje2 = sqrt(d_je2^(-c)) * (randn(N, 1) + 1i * randn(N, 1));
             hje3 = sqrt(d_je3^(-c)) * (randn(N, 1) + 1i * randn(N, 1));

             %phase shift due to RIS           
             Phi1 = exp(-1i * angle(hsr .* hru1));
             Phi2 = exp(-1i * angle(hsr .* hru2));

             %snr of users
             snr1_nj = ((abs(hsr .* hru1 .* Phi1).^2 * (a1 * Pt))+(abs(hru1).^2 * (a1 * Pt)))./No
             snr2_nj = ((abs(hsr .* hru2 .* Phi2).^2 * (a2 * Pt))+(abs(hru2).^2 * (a2 * Pt)))./ (No +((abs(hsr .* hru1 .* Phi1).^2 * (a1 * Pt))+(abs(hru1).^2 * (a1 * Pt)))) 
             
             snr1_uj=((abs(hsr .* hru1 .* Phi1).^2 * (a1 * Pt))+(abs(hru1).^2 * (a1 * Pt)))./(No+(abs(hju1).^2 * Pj))
             snr2_uj=((abs(hsr .* hru2 .* Phi2).^2 * (a2 * Pt))+(abs(hru2).^2 * (a2 * Pt)))./(No+(abs(hju2).^2 * Pj)+((abs(hsr .* hru1 .* Phi1).^2 * (a1 * Pt))+(abs(hru1).^2 * (a1 * Pt))))

             snr1_kj=((abs(hsr .* hru1 .* Phi1).^2 * (a1 * Pt))+(abs(hru1).^2 * (a1 * Pt)))./(No)
             snr2_kj=((abs(hsr .* hru2 .* Phi2).^2 * (a2 * Pt))+(abs(hru2).^2 * (a2 * Pt)))./(No+((abs(hsr .* hru1 .* Phi1).^2 * (a1 * Pt))+(abs(hru1).^2 * (a1 * Pt))))
             %snr of eavesdroppers   
             snre_nj = (abs(hse1).^2 * (a2 * Pt))./No;
             
             snre_uj=(abs(hse1).^2 * (a2 * Pt))./(No+(abs(hje1).^2 * Pj));

             snre_kj=(abs(hse1).^2 * (a2 * Pt))./(No+abs(hje1).^2*Pj);
             %secrecy rate calculations

             sr_nj1 = max(0.5 * (log2(1 + mean(snr1_nj)) - log2(1 + mean(snre_nj))), 0)
             sr_nj2 = max(0.5 * (log2(1 + mean(snr2_nj)) - log2(1 + mean(snre_nj))), 0)

             sr_uj1 = max(0.5 * (log2(1 + mean(snr1_uj)) - log2(1 + mean(snre_uj))), 0)
             sr_uj2 = max(0.5 * (log2(1 + mean(snr2_uj)) - log2(1 + mean(snre_uj))), 0)

             sr_kj1 = max(0.5 * (log2(1 + mean(snr1_kj)) - log2(1 + mean(snre_kj))), 0)
             sr_kj2 = max(0.5 * (log2(1 + mean(snr2_kj)) - log2(1 + mean(snre_kj))), 0)

             cap_nj1(jj, ii) = sr_nj1;
             cap_nj2(jj, ii) = sr_nj2;

             cap_uj1(jj, ii) = sr_uj1;
             cap_uj2(jj, ii) = sr_uj2;
             
             cap_kj1(jj, ii) = sr_kj1;
             cap_kj2(jj, ii) = sr_kj2;


             asnr1_nj = ((abs(hsr .* hru1 .* Phi1).^2 * (a1 * Pt))+(abs(hru1).^2 * (a1 * Pt)))./No
             asnr2_nj = ((abs(hsr .* hru2 .* Phi2).^2 * (a2 * Pt))+(abs(hru2).^2 * (a2 * Pt)))./ (No +((abs(hsr .* hru1 .* Phi1).^2 * (a1 * Pt))+(abs(hru1).^2 * (a1 * Pt)))) 
             
             asnr1_uj=((abs(hsr .* hru1 .* Phi1).^2 * (a1 * Pt))+(abs(hru1).^2 * (a1 * Pt)))./(No+(abs(hju1).^2 * Pj))
             asnr2_uj=((abs(hsr .* hru2 .* Phi2).^2 * (a2 * Pt))+(abs(hru2).^2 * (a2 * Pt)))./(No+(abs(hju2).^2 * Pj)+((abs(hsr .* hru1 .* Phi1).^2 * (a1 * Pt))+(abs(hru1).^2 * (a1 * Pt))))

             asnr1_kj=((abs(hsr .* hru1 .* Phi1).^2 * (a1 * Pt))+(abs(hru1).^2 * (a1 * Pt)))./(No)
             asnr2_kj=((abs(hsr .* hru2 .* Phi2).^2 * (a2 * Pt))+(abs(hru2).^2 * (a2 * Pt)))./(No+((abs(hsr .* hru1 .* Phi1).^2 * (a1 * Pt))+(abs(hru1).^2 * (a1 * Pt))))
             %snr of eavesdroppers   
             asnre_nj =[(abs(hse1).^2 * (a2 * Pt))./No]+[(abs(hse2).^2 * (a2 * Pt))./No]+[(abs(hse3).^2 * (a2 * Pt))./No];
             
             asnre_uj=[(abs(hse1).^2 * (a2 * Pt))./(No+(abs(hje1).^2 * Pj))]+[(abs(hse3).^2 * (a2 * Pt))./(No+(abs(hje3).^2 * Pj))]+[(abs(hse2).^2 * (a2 * Pt))./(No+(abs(hje2).^2 * Pj))];

             asnre_kj=[(abs(hse1).^2 * (a2 * Pt))./(No+abs(hje1).^2*Pj)]+[(abs(hse2).^2 * (a2 * Pt))./(No+abs(hje2).^2*Pj)]+[(abs(hse3).^2 * (a2 * Pt))./(No+abs(hje3).^2*Pj)];
             %secrecy rate calculations

             asr_nj1 = max(0.5 * (log2(1 + mean(asnr1_nj)) - log2(1 + mean(asnre_nj))), 0)
             asr_nj2 = max(0.5 * (log2(1 + mean(asnr2_nj)) - log2(1 + mean(asnre_nj))), 0)

             asr_uj1 = max(0.5 * (log2(1 + mean(asnr1_uj)) - log2(1 + mean(asnre_uj))), 0)
             asr_uj2 = max(0.5 * (log2(1 + mean(asnr2_uj)) - log2(1 + mean(asnre_uj))), 0)

             asr_kj1 = max(0.5 * (log2(1 + mean(asnr1_kj)) - log2(1 + mean(asnre_kj))), 0)
             asr_kj2 = max(0.5 * (log2(1 + mean(asnr2_kj)) - log2(1 + mean(asnre_kj))), 0)

             acap_nj1(jj, ii) = asr_nj1;
             acap_nj2(jj, ii) = asr_nj2;

             acap_uj1(jj, ii) = asr_uj1;
             acap_uj2(jj, ii) = asr_uj2;
             
             acap_kj1(jj, ii) = asr_kj1;
             acap_kj2(jj, ii) = asr_kj2;

        end
    end

    % Averaging over iterations
    C_RIS1_nojam1(a_index, :) = mean(cap_nj1, 1);
    C_RIS1_nojam2(a_index, :) = mean(cap_nj2, 1);

    C_RIS1_unknownjam1(a_index, :) = mean(cap_uj1, 1);
    C_RIS1_unknownjam2(a_index, :) = mean(cap_uj2, 1);

    C_RIS1_knownjam1(a_index, :) = mean(cap_kj1, 1);
    C_RIS1_knownjam2(a_index, :) = mean(cap_kj2, 1);

    aC_RIS1_nojam1(a_index, :) = mean(acap_nj1, 1);
    aC_RIS1_nojam2(a_index, :) = mean(acap_nj2, 1);

    aC_RIS1_unknownjam1(a_index, :) = mean(acap_uj1, 1);
    aC_RIS1_unknownjam2(a_index, :) = mean(acap_uj2, 1);

    aC_RIS1_knownjam1(a_index, :) = mean(acap_kj1, 1);
    aC_RIS1_knownjam2(a_index, :) = mean(acap_kj2, 1);
end


% Plotting the results
set(0, 'DefaultTextFontName', 'Times New Roman');
set(0, 'DefaultAxesFontName', 'Times New Roman');

fig = figure;
hold on;
grid on

% Colors for each power allocation scheme
colors = ['g', 'r', 'b', 'm', 'k', 'c'];

plot(snr_dB, C_RIS1_nojam2, '-ms', 'DisplayName', 'user 2 sr nojamming in non collabrative case', 'LineWidth', 2.2, 'Color', 'g');
plot(snr_dB, C_RIS1_unknownjam2, '-m^', 'DisplayName', 'user 2 sr unknownjamming in non collabrative case', 'LineWidth', 2.2, 'Color', 'g');
plot(snr_dB, C_RIS1_knownjam2, '-mx', 'DisplayName', 'user 2 sr knownjamming in non collabrative case', 'LineWidth', 2.2, 'Color', 'g');


plot(snr_dB, aC_RIS1_nojam2, '-ms', 'DisplayName', 'user 2 sr nojamming in collabrative case', 'LineWidth', 2.2, 'Color', 'r');
plot(snr_dB, aC_RIS1_unknownjam2, '-m^', 'DisplayName', 'user 2 sr unknownjamming in collabrative case', 'LineWidth', 2.2, 'Color', 'r');
plot(snr_dB, aC_RIS1_knownjam2, '-mx', 'DisplayName', 'user 2 sr knownjamming in collabrative case', 'LineWidth', 2.2, 'Color', 'r');

xlabel('SNR (dB)', 'FontWeight', 'bold');
ylabel('Secrecy Rate (bits/sec/Hz)', 'FontWeight', 'bold');
title('secrecy rate vs snr at user 2 in both collabrative and non collabrative cases')
grid on;
legend('Location', 'best', 'NumColumns', 2, 'FontSize', 8);

% Save the plot
exportgraphics(fig, 'myplot.jpeg', 'Resolution', 800);

















