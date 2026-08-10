clc
clear all
close all

snr_dB = 40;    % Total power in dB
c = 2.2;            % Pathloss component
No = 1;             % Noise power
iterations = 400;    % Number of iterations for Monte Carlo simulation
N = 30;            % Number of channel realizations

% Distances
d_sr = 50;      % Distance from Source to RIS
d_ru1 = 10;          % Distance from RIS to Destination
d_ru2 = 20;    % Distance from Source to E1
d_su1=60;
d_su2=70;
d_se1 = 30;      % Distance from RIS to E1
d_se2 = 40;          % Distance from Jammer to E1
d_se3 = 50;      % Distance from Jammer to Destination
d_ju1=50;      % Distance from Source to Destination
d_ju2=30;
d_je1=0:3:30
d_je2=0:3:30
d_je3=0:3:30;

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

         for ii = 1:length(d_je1)
             Ps = 10^(snr_dB/ 10);  % Power conversion from dB to linear scale
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
             hje1 = sqrt((d_je1(ii))^(-c)) * (randn(N, 1) + 1i * randn(N, 1));
             hje2 = sqrt(d_je2(ii)^(-c)) * (randn(N, 1) + 1i * randn(N, 1));
             hje3 = sqrt(d_je3(ii)^(-c)) * (randn(N, 1) + 1i * randn(N, 1));

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
             snre_nj =[(abs(hse1).^2 * (a2 * Pt))./No]+[(abs(hse2).^2 * (a2 * Pt))./No]+[(abs(hse3).^2 * (a2 * Pt))./No];
             
             snre_uj=[(abs(hse1).^2 * (a2 * Pt))./(No+(abs(hje1).^2 * Pj))]+[(abs(hse3).^2 * (a2 * Pt))./(No+(abs(hje3).^2 * Pj))]+[(abs(hse2).^2 * (a2 * Pt))./(No+(abs(hje2).^2 * Pj))];

             snre_kj=[(abs(hse1).^2 * (a2 * Pt))./(No+abs(hje1).^2*Pj)]+[(abs(hse2).^2 * (a2 * Pt))./(No+abs(hje2).^2*Pj)]+[(abs(hse3).^2 * (a2 * Pt))./(No+abs(hje3).^2*Pj)];
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

        end
    end

    % Averaging over iterations
    C_RIS1_nojam1(a_index, :) = mean(cap_nj1, 1);
    C_RIS1_nojam2(a_index, :) = mean(cap_nj2, 1);

    C_RIS1_unknownjam1(a_index, :) = mean(cap_uj1, 1);
    C_RIS1_unknownjam2(a_index, :) = mean(cap_uj2, 1);

    C_RIS1_knownjam1(a_index, :) = mean(cap_kj1, 1);
    C_RIS1_knownjam2(a_index, :) = mean(cap_kj2, 1);
end


% Plotting the results
set(0, 'DefaultTextFontName', 'Times New Roman');
set(0, 'DefaultAxesFontName', 'Times New Roman');

fig = figure;
hold on;
grid on

% Colors for each power allocation scheme
colors = ['g', 'r', 'b', 'm', 'k', 'c'];

plot(d_je1, C_RIS1_nojam2, '-m^', 'DisplayName', 'user 2 sr nojamming', 'LineWidth', 2.2, 'Color', 'g');

plot(d_je1, C_RIS1_unknownjam2, '-m^', 'DisplayName', 'user 2 sr unknownjamming', 'LineWidth', 2.2, 'Color', 'r');
 
plot(d_je1, C_RIS1_knownjam2, '-m^', 'DisplayName', 'user 2 sr knownjamming', 'LineWidth', 2.2, 'Color', 'k');

xlabel('d_je', 'FontWeight', 'bold');
ylabel('Secrecy Rate (bits/sec/Hz)', 'FontWeight', 'bold');
title("secrecy rate of user 2 vs d_je1 in collabrative case")
grid on;
legend('Location', 'best', 'NumColumns', 2, 'FontSize', 8);

% Save the plot
exportgraphics(fig, 'myplot.jpeg', 'Resolution', 800);




