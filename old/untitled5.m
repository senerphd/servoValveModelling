%% OPTİMİZASYON BÖLÜMÜ (k > 0 Kısıtı Eklenmiş)
fprintf('Optimizasyon başlatılıyor (k > 0 Kısıtı ile)...\n');
close all 

% 1. Optimizasyon için Veri Hazırlığı
subset_range_1 = 1:25000; 
subset_range_2 = 25001:54900; 

% Sızıntı verisi
I_leak_raw = leakageCurrent_filtered_A(subset_range_1);
Q_leak_raw = leakageFlow_filtered_lpm(subset_range_1);

% Basınç verisi
I_press_raw = pressureGainCurrent_A(subset_range_2);
P_press_raw = pressureGainPressure_bar(subset_range_2);

% 2. Akım Limiti ve Filtreleme
limit_I = 3e-3; % Arama aralığı (Modelin geçerli olduğu dar bölge)

idx_leak = abs(I_leak_raw) <= limit_I;
I_leak_opt = I_leak_raw(idx_leak);
Q_leak_opt = Q_leak_raw(idx_leak) * lpm2m3Ps; 

idx_press = abs(I_press_raw) <= limit_I;
I_press_opt = I_press_raw(idx_press);
P_press_opt = P_press_raw(idx_press) * bar2Pa; 

% Basınç verisi Unique hale getirilir (Enterpolasyon için)
[I_press_unique, unique_idx] = unique(I_press_opt);
P_press_unique = P_press_opt(unique_idx);

% 3. Amaç Fonksiyonu (Cost Function)
objectiveFunc = @(x) calculateTotalError(x, ...
    I_leak_opt, Q_leak_opt, ...
    I_press_unique, P_press_unique, ... 
    I_press_opt, P_press_opt, ...       
    Ps, Pr, K, Q_leak_net_max, Q_tare);

% 4. Optimizasyonu Çalıştır
options = optimset('Display','iter','TolX',1e-12,'MaxFunEvals',200);
% Alt sınırı çok küçük yapmayın, fiziksel olarak x0 0 olamaz.
[x0_opt, min_error] = fminbnd(objectiveFunc, 0.00001, 0.003, options);

fprintf('Optimizasyon Tamamlandı.\n');
fprintf('Optimize edilmiş x0: %.6f Amper (%.3f mA)\n', x0_opt, x0_opt*1e3);

%% 5. Modelin Kurulması ve Kontrolü
x0 = x0_opt; 

% x0'a karşılık gelen PL değerini bul
PL_at_x0 = interp1(I_press_unique, P_press_unique, x0, 'linear', 'extrap');
PL_at_x0 = abs(PL_at_x0); 

% k hesabı
term = (Ps - Pr + PL_at_x0) / (Ps - Pr - PL_at_x0);
if term < 0, term = 0; end
k = (1/2) * sqrt(term) - 1;

fprintf('Bu x0 değerindeki Deneysel Yük Basıncı (PL): %.2f Bar\n', PL_at_x0 / bar2Pa);
fprintf('Model Parametresi k: %.4f (k > 0 olmalı)\n', k);

% --- K_sc (Scaling Factor) Otomatik Hesaplama ---
Calculated_Peak = (2 * K * sqrt(Ps - Pr) * x0); % Modelin teorik zirvesi
Raw_Scale_Factor = Q_leak_net_max / Calculated_Peak; % Veriye oranla
Manual_Correction = (1.3/1.1); % Sizin manuel düzeltmeniz (isteğe bağlı)
K_sc = Raw_Scale_Factor * Manual_Correction; 
fprintf('K_sc (Scaling Factor): %.4f\n', K_sc);
% ------------------------------------------------

% Simülasyon
I_range = linspace(-25e-3, 25e-3, 10e3); 
f_x = (1 + abs(I_range)/x0).^2 .* (1 + k * abs(I_range)/x0).^2;

Qs_model = K_sc * 2 * K * sqrt(Ps - Pr) * (x0 + abs(I_range)) .* (1 + f_x).^(-1/2); 
Qs_model = Qs_model + Q_tare; 

PL_model = ((f_x - 1) ./ (f_x + 1)) * (Ps - Pr) .* sign(I_range);

%% 6. Grafik Çizimi
figure('Name', 'Constraint Optimized Model (k>0)', 'Color', 'w');

% Grafik 1: Sızıntı
subplot(2,1,1);
plot(I_range*1e3, Qs_model*m3Ps2lpm, 'r-', 'LineWidth', 2.5); hold on;
plot(I_leak_raw*1e3, Q_leak_raw, '.','Color',[0.6 0.6 0.6],'LineWidth',0.5); 
xline(limit_I*1e3, '--k'); xline(-limit_I*1e3, '--k');
title(['Sızıntı Debisi (x0 = ' num2str(x0*1e3, '%.3f') ' mA)']);
xlabel('Akım (mA)'); ylabel('Debi (l/min)'); grid on; xlim([-25 25]);

% Grafik 2: Basınç
subplot(2,1,2);
plot(I_range*1e3, PL_model/bar2Pa, 'r-', 'LineWidth', 2); hold on;
plot(I_press_raw*1e3, P_press_raw, 'b.', 'MarkerSize', 4); 
xline(limit_I*1e3, '--k'); xline(-limit_I*1e3, '--k');
title(['Basınç Hassasiyeti (k = ' num2str(k, '%.3f') ')']);
xlabel('Akım (mA)'); ylabel('Yük Basıncı (bar)'); grid on; xlim([-25 25]);

%% Yardımcı Fonksiyon (DÜZELTİLMİŞ)
function total_error = calculateTotalError(x0_try, I_leak, Q_leak_data, I_press_uni, P_press_uni, I_press_data, P_press_data, Ps, Pr, K, Q_net_max, Q_tare)
    
    % 1. PL değerini çek
    try
        PL_try = interp1(I_press_uni, P_press_uni, x0_try, 'linear', 'extrap');
        PL_try = abs(PL_try); 
    catch
        total_error = 1e9; return; 
    end
    
    % 2. k değerini hesapla
    % Matematiksel limit kontrolü (karekök içi negatif olmamalı)
    if (Ps - Pr - PL_try) <= 0
         total_error = 1e9; return; % Payda sıfır veya negatif -> Ceza
    end
    
    term = (Ps - Pr + PL_try)/(Ps - Pr - PL_try);
    k_try = 0.5 * sqrt(term) - 1;
    
    % --- KRİTİK KISIT ---
    % Eğer k <= 0 ise bu x0 değeri geçersizdir.
    % Algoritmayı buradan uzaklaştırmak için çok büyük bir hata döndür.
    if k_try <= 0.001 
        total_error = 1e12; % Çok büyük ceza puanı
        return; % Fonksiyondan çık, hesaplamaya devam etme
    end
    % --------------------
    
    % 3. Modelleri Hesapla
    f_x_leak = (1 + abs(I_leak)/x0_try).^2 .* (1 + k_try * abs(I_leak)/x0_try).^2;
    

    
    Q_model_net = 2 * K * sqrt(Ps - Pr) * (x0_try + abs(I_leak)) .* (1 + f_x_leak).^(-1/2);
    Q_model_total = Q_model_net + Q_tare; 
    
    f_x_press = (1 + abs(I_press_data)/x0_try).^2 .* (1 + k_try * abs(I_press_data)/x0_try).^2;
    PL_model = ((f_x_press - 1) ./ (f_x_press + 1)) * (Ps - Pr) .* sign(I_press_data);
    
    % 4. Hata Hesabı
    err_leak = sum(((Q_model_total - Q_leak_data) / max(Q_leak_data)).^2);
    err_press = sum(((PL_model - P_press_data) / Ps).^2);
    
    % Ağırlıklar: %90 Sızıntı, %10 Basınç (Sızıntıyı daha iyi oturtmak için)
    total_error = 10 * err_leak + 1 * err_press;
end