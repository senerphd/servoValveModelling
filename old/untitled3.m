%% Tamburrano & Merritt Tabanlı Fiziksel Köprü Modeli
% Bu model, sızıntı bölgesindeki akışı LAMINAR, ana akışı TÜRBÜLANSLI kabul eder.
% Eryılmaz modelindeki "Volkan Ağzı" (çökme) sorununu fiziksel olarak yok eder.

clearvars -except flowGainedited pressureGainedited leakage500edited psi2Pa lpm2m3Ps m3Ps2lpm bar2Pa
clc; close all;

%% 1. Veri Hazırlığı ve Sabitler
if ~exist('psi2Pa','var'), psi2Pa = 6894.75729; end
if ~exist('lpm2m3Ps','var'), lpm2m3Ps = 1/60e3; end
if ~exist('bar2Pa','var'), bar2Pa = 1e5; end
m3Ps2lpm = 1/lpm2m3Ps;

% Verileri Çek
I_leak = leakage500edited.CurrentA; 
Q_leak = leakage500edited.Flowlpm; 
I_press = pressureGainedited.CurrentA; 
P_press_bar = pressureGainedited.Pressurebar; 

% Sistem Basınçları
Ps = 3000 * psi2Pa; 
Pr = 0;             

%% 2. Parametrelerin Belirlenmesi
% A. Türbülanslı Kazanç (K_turb): Flow Gain testinden
Q_max_test = max(flowGainedited.Flowlpm * lpm2m3Ps);
I_max = 25e-3;
DeltaP_FG = 1000 * psi2Pa; 
K_turb = Q_max_test / (I_max * sqrt(DeltaP_FG)); 

% B. Laminar Sızıntı İletkenliği (C_null): Null sızıntısından
Q_tare = 0.64 * lpm2m3Ps; % Ofset (Pilot flow - Tamburrano Q_leak,I)
Q_leak_peak_net = max(Q_leak * lpm2m3Ps) - Q_tare; 

% Null noktasında (x=0) basınç düşümü Ps/2'dir (Tamburrano Eq. 16)
% Q_peak = C_null * Ps (Toplam sızıntı)
C_null = Q_leak_peak_net / Ps; 

% C. Sızıntı Kapanma Hızı (x_lap): EN ÖNEMLİ AYAR
% Bu parametre sızıntı grafiğinin genişliğini belirler.
% Tamburrano makalesindeki "Overlap" kavramının elektriksel karşılığıdır.
x_lap = 0.00035; % [Amper] <-- Burayı değiştirerek genişliği ayarlayın!

fprintf('--- Model Parametreleri ---\n');
fprintf('K_turb : %.4e (Ana Akış Gücü)\n', K_turb);
fprintf('C_null : %.4e (Sızıntı Tepesi)\n', C_null);
fprintf('x_lap  : %.4e (Sızıntı Genişliği)\n', x_lap);

%% 3. Simülasyon (Fiziksel Denge Çözümü)
I_sim = linspace(-25e-3, 25e-3, 2000);
Qs_model = zeros(size(I_sim));
PL_model = zeros(size(I_sim));

for i = 1:length(I_sim)
    x = I_sim(i);
    x_abs = abs(x);
    
    % --- İLETKENLİKLER (Conductances) ---
    % 1. Kapanan Hat (Laminar Sızıntı): Spool hareket ettikçe direnci artar.
    % G_close = C_null * (x_lap / (x_lap + x))
    G_close = C_null * (x_lap / (x_lap + x_abs)); 
    
    % 2. Açılan Hat (Türbülanslı + Laminar Sızıntı):
    % Açılan tarafta hem ana orifis açılır hem de null sızıntısı devam eder.
    % Bu "paralel" yapı, x=0'da sürekliliği sağlar.
    G_open_turb = K_turb * x_abs; % Türbülanslı kısım
    
    % --- NODE A (Besleme -> A -> Tank) BASINÇ DENGESİ ---
    % Tamburrano Eq. 5 ve Eq. 11-12 temel alınarak:
    % P_L hesabı için Açılan ve Kapanan alanların oranını (Alpha) kullanıyoruz.
    
    Alpha = (G_open_turb + G_close) / G_close;
    
    % Basınç Kazancı (Pressure Gain)
    if x >= 0
        PL = (Alpha^2 - 1) / (Alpha^2 + 1) * (Ps - Pr);
    else
        PL = -1 * (Alpha^2 - 1) / (Alpha^2 + 1) * (Ps - Pr);
    end
    
    PL_model(i) = PL;
    
    % --- SIZINTI HESABI (Qs) ---
    % Tamburrano Eq. 6: Toplam Sızıntı = Q2 + Q4
    % Bizim modelde: Qs = Q(P->A) + Q(P->B)
    % Pa ve Pb basınçlarını bulalım
    Pa = (Ps + Pr + PL)/2;
    Pb = (Ps + Pr - PL)/2;
    
    % P->A Akışı (Giriş):
    % Hem laminar sızıntı hem türbülanslı ana akış katkısı (Eq. 1)
    Q_PA = G_close * (Ps - Pa) + K_turb * x_abs * (x>=0) * sqrt(Ps - Pa);
    
    % P->B Akışı (Giriş):
    % x>0 iken B tarafı kapanıyor (Sadece laminar sızıntı)
    Q_PB = G_close * (Ps - Pb) + K_turb * x_abs * (x<0) * sqrt(Ps - Pb);
    
    Qs_model(i) = Q_PA + Q_PB;
end

% Tare Flow Ekle
Qs_total = Qs_model + Q_tare;

%% 4. Görselleştirme
figure('Color','w', 'Position', [100 100 1000 700]);

% Sızıntı
subplot(2,1,1);
plot(I_leak*1e3, Q_leak, 'k.', 'Color',[0.6 0.6 0.6]); hold on;
plot(I_sim*1e3, Qs_total*m3Ps2lpm, 'r-', 'LineWidth', 2.5);
title(['Tamburrano Fiziksel Modeli (Volkansız Sivri Tepe) | x\_lap: ' num2str(x_lap)]);
ylabel('Debi (lpm)'); xlabel('Akım (mA)'); grid on;
xlim([-10 10]);
legend('Test Verisi', 'Fiziksel Model');

% Basınç
subplot(2,1,2);
plot(I_press*1e3, P_press_bar, 'b.', 'MarkerSize',4); hold on;
plot(I_sim*1e3, PL_model/bar2Pa, 'r-', 'LineWidth', 2.5);
title('Basınç Kazancı (Doğal Fiziksel Geçiş)');
ylabel('Basınç (bar)'); xlabel('Akım (mA)'); grid on;
xlim([-2 2]);
legend('Test Verisi', 'Fiziksel Model');