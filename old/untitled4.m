%% 1. Sisteme Ait Giriş Verileri
Ps = 3000 * psi2Pa; 
Pr = 0;             
I_max = 25e-3; 
close all 

% Akış Kazancı (K) Hesabı
Q_max = max(flowGainFlow_lpm * lpm2m3Ps); 
DeltaP_FlowGain = 1000 * psi2Pa; 
K = (sqrt(2) * Q_max) / (sqrt(DeltaP_FlowGain) * I_max);

% Tare Flow (Ofset)
Q_tare = 0.64 * lpm2m3Ps; 

% Pressure Gain Verisi (Unique hale getirilmiş - Lookup için şart)
[uniq_I, uniq_idx] = unique(pressureGainCurrent_A);
uniq_P = pressureGainPressure_bar(uniq_idx) * bar2Pa;

% Test Sızıntı Tepe Noktası (Net)
Q_leak_net_target = (1.32 * lpm2m3Ps) - Q_tare; 

%% 2. OTOMATİK OPTİMİZASYON (x0 Bulucu)
fprintf('Optimizasyon Başlıyor... En iyi x0 aranıyor...\n');

% Başlangıç Tahmini (Sizin son denediğiniz değer)
initial_guess_x0 = 3.5e-4; 

% Cost Function Tanımı (Hatayı hesaplayan fonksiyon)
% Bu fonksiyon aşağıda tanımlıdır, script sonunda yer alır.
cost_func = @(x) calculate_error(x, leakageCurrent_A, leakageFlow_lpm, ...
                                 Ps, Pr, K, Q_tare, uniq_I, uniq_P, ...
                                 Q_leak_net_target, lpm2m3Ps);

% MATLAB fminsearch ile en iyi x0'ı bul
options = optimset('Display','iter','TolX',1e-8);
best_x0 = fminsearch(cost_func, initial_guess_x0, options);

fprintf('------------------------------------------------\n');
fprintf('BULUNAN EN İYİ x0: %.4e A (%.4f mA)\n', best_x0, best_x0*1e3);
fprintf('------------------------------------------------\n');

%% 3. Bulunan x0 ile Final Model Parametrelerini Hesapla
% (Simülasyon ve Çizim İçin)

% 1. k (Sızıntı Katsayısı)
PL_at_x0 = interp1(uniq_I, uniq_P, best_x0, 'linear', 'extrap');
limit_val = 0.65 * (Ps - Pr);
if PL_at_x0 < limit_val, PL_at_x0 = limit_val + 1e5; end
k_opt = 0.5 * sqrt((Ps - Pr + PL_at_x0)/(Ps - Pr - PL_at_x0)) - 1; 

% 2. Correction Factor (C_leak)
Q_model_raw_peak = 2 * K * sqrt(Ps - Pr) * best_x0;
C_leak_opt = Q_leak_net_target / Q_model_raw_peak;

fprintf('Final k Değeri: %.4f\n', k_opt);
fprintf('Final C_leak  : %.4f\n', C_leak_opt);

%% 4. Model Simülasyonu ve Görselleştirme
I_range = linspace(-25e-3, 25e-3, 1000); 

% Sürekli f(x) fonksiyonu
f_x = (1 + abs(I_range)/best_x0).^2 .* (1 + k_opt * abs(I_range)/best_x0).^2;

% Sızıntı Debisi Modeli (Qs)
Qs_model = 2 * (K * C_leak_opt) * sqrt(Ps - Pr) * (best_x0 + abs(I_range)) .* (1 + f_x).^(-1/2);

% Yük Basıncı Modeli (PL) 
PL_model = ((f_x - 1) ./ (f_x + 1)) * (Ps - Pr) .* sign(I_range);

% --- GRAFİK ÇİZİMİ ---
figure('Color', 'w', 'Position', [100 100 1000 700]);

% Grafik 1: Sızıntı
subplot(2,1,1);
plot(leakageCurrent_A*1e3, leakageFlow_lpm, 'ro', 'MarkerSize', 3, 'MarkerFaceColor', 'r'); hold on;
plot(I_range*1e3, (Qs_model*m3Ps2lpm) + (Q_tare*m3Ps2lpm), 'b-', 'LineWidth', 2.5); 
title(['OTOMATİK FIT: Sızıntı Debisi (x0 = ' num2str(best_x0*1e3, '%.4f') ' mA)']);
xlabel('Akım (mA)'); ylabel('Debi (l/min)'); grid on; 
legend('Test Verisi', 'Optimize Edilmiş Model', 'Location','South');
xlim([-15 15]);

% Grafik 2: Basınç
subplot(2,1,2);
plot(pressureGainCurrent_A*1e3, pressureGainPressure_bar, 'bo', 'MarkerSize', 3); hold on;
plot(I_range*1e3, PL_model*1/bar2Pa, 'r-', 'LineWidth', 2.5); 
title('Basınç Hassasiyeti (Model Tutarlılık Kontrolü)');
xlabel('Akım (mA)'); ylabel('Yük Basıncı (bar)'); grid on;


%% --- YARDIMCI FONKSİYONLAR ---
% Bu fonksiyon x0 denemeleri yaparak hatayı hesaplar
function error_val = calculate_error(x_trial, I_test, Q_test_lpm, Ps, Pr, K, Q_tare, I_press, P_press, Q_target, unit_conv)
    
    % Negatif veya anlamsız x0 gelirse cezalandır
    if x_trial <= 1e-6 || x_trial > 1e-3
        error_val = 1e9; 
        return;
    end

    % 1. Bu x0 için k hesapla
    PL_val = interp1(I_press, P_press, x_trial, 'linear', 'extrap');
    limit = 0.65 * (Ps - Pr);
    if PL_val < limit
        PL_val = limit + 1e5; % Limit düzeltmesi
    end
    k_trial = 0.5 * sqrt((Ps - Pr + PL_val)/(Ps - Pr - PL_val)) - 1;

    % 2. Bu x0 için Yükseklik Katsayısı (C_leak) hesapla
    Q_raw = 2 * K * sqrt(Ps - Pr) * x_trial;
    C_corr = Q_target / Q_raw;

    % 3. Modeli TEST NOKTALARINDA oluştur (I_test vektörü kullanılarak)
    f_x_trial = (1 + abs(I_test)/x_trial).^2 .* (1 + k_trial * abs(I_test)/x_trial).^2;
    Qs_model_net = 2 * (K * C_corr) * sqrt(Ps - Pr) * (x_trial + abs(I_test)) .* (1 + f_x_trial).^(-1/2);
    
    % Toplam Model (Net + Tare)
    Qs_model_total_lpm = (Qs_model_net + Q_tare) / unit_conv;

    % 4. Hatayı Hesapla (RMS Error)
    % Farkın karesini alıp topluyoruz. Amacımız bunu 0 yapmak.
    residuals = Q_test_lpm - Qs_model_total_lpm;
    error_val = sum(residuals.^2); 
end