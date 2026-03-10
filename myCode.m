% Servo Valf Statik Özelliklerinin Test Verileri ile Modellenmesi
% Verilerin Yüklenmesi

%%

clear all
format compact
format long

psi2Pa = 6894.75729;
lpm2m3Ps = 1/60e3; 
m3Ps2lpm = 1/ lpm2m3Ps; 
bar2Pa = 1e5; 
Pa2Bar = 1/1e5; 
addpath('C:\_PhD\4_Tez\3_Tests\1_TestData\_lib');
% constants; 

%%
load("servoData.mat");
close all; 
%% 

% Flow Gain Test
figure
plot(flowGainCurrent_A*1e3,flowGainFlow_lpm,"-")
xlabel("Current (mA)")
ylabel("Flow (lpm)")
grid on

% Pressure Gain Test
figure
plot(pressureGainCurrent_A*1e3,pressureGainPressure_bar,"-")
xlabel("Current (mA)")
ylabel("Pressure (bar)")
grid on

% Leakage Test
figure
plot(leakageCurrent_filtered_A*1e3,leakageFlow_filtered_lpm)
hold on 
xlabel("Current (mA)")
ylabel("Flow (lpm)")
grid on

%% 1. Sisteme Ait Giriş Verileri (Katalog ve Test Koşulları)
close all 

Ps = 2980 * psi2Pa; % 3000psi yapıldı  % [bar] Besleme basıncı (Örnek değer, kendi test verinizle güncelleyin)
Pr = 0;             % [bar] Dönüş hattı basıncı

% Katalogdan veya Testten Alınacak Veriler:
Q_max = max(flowGainFlow_lpm * lpm2m3Ps); 
I_max = 25e-3; 
DeltaP_FlowGain = 1032 * psi2Pa; % Flow gain testleri 1000 psi deltaP'de yapıldı

% Flow gain hesaplama 
K = (sqrt(2) * Q_max) / (sqrt(DeltaP_FlowGain) * I_max) % [m3Ps/sqrt(Pa).A]

Q_tare  = 0.64 * lpm2m3Ps; % Test verisinden elde edildi
Qs_test =  leakageFlow_filtered_lpm * lpm2m3Ps - Q_tare; % Test verisinden elde edildi.
Qs_at_x0 = max(Qs_test)

x0 = Qs_at_x0 / (sqrt(2)*K*sqrt(Ps - Pr));

PL_at_x0 = 193 * bar2Pa;

k = (1/2) * sqrt((Ps - Pr + PL_at_x0)/(Ps - Pr - PL_at_x0)) - 1; 
%% 4. Model Simülasyonu ve Görselleştirme
I_range = linspace(-25e-3, 25e-3, 10e3); % Akım aralığı [A]

% Sürekli f(x) fonksiyonu
f_x = (1 + abs(I_range)/x0).^2 .* (1 + k * abs(I_range)/x0).^2;

% Sızıntı Debisi Modeli (Qs) 
Qs_model =  2 * K * sqrt(Ps - Pr) * (x0 + abs(I_range)) .* (1 + f_x).^(-1/2); 

% Yük Basıncı Modeli (PL) 
PL_model = ((f_x - 1) ./ (f_x + 1)) * (Ps - Pr) .* sign(I_range);

Qs_model = Qs_model + Q_tare; 

% Grafik Çizimi
figure;
subplot(2,1,1);

plot(I_range*1e3, Qs_model*m3Ps2lpm, 'r-', 'LineWidth', 4.5); hold on;
plot(leakageCurrent_filtered_A*1e3, leakageFlow_filtered_lpm, '.','Color',[0.6 0.6 0.6],'LineWidth',0.2); % Test verinizi buraya ekleyin

title('Sızıntı Debisi (Internal Leakage)');
xlabel('Akım (mA)'); ylabel('Debi (l/min)'); grid on;

subplot(2,1,2);
plot(I_range*1e3, PL_model*1/bar2Pa, 'r-', 'LineWidth', 1.5); hold on;
plot(pressureGainCurrent_A*1e3, pressureGainPressure_bar, 'bo'); % Test verinizi buraya ekleyin
title('Basınç Hassasiyeti (Pressure Sensitivity)');
xlabel('Akım (mA)'); ylabel('Yük Basıncı (bar)'); grid on;

%% 
clc 
format long
fprintf('------------------------------------- \n');
fprintf('K değeri \t: \t %1.8f \n',K);
fprintf('Qs_at_x0 \t: \t %1.8f \n',Qs_at_x0);
fprintf('x0 \t \t \t: \t %1.8f \n',x0);
fprintf('k \t \t \t:  \t %1.8f \n',k);
fprintf('------------------------------------- \n');

%% Dara Akış Modelleme 
K_tare = Q_tare/sqrt(Ps - Pr); 
