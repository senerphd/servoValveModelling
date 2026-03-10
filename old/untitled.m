%% OPTİMİZASYON BÖLÜMÜ: Nonlinear Least Squares (Trust-Region-Reflective)
fprintf('Optimizasyon başlatılıyor (Yöntem: lsqnonlin | Ağırlık: %%80 Sızıntı, %%20 Basınç)...\n');

% 1. Veri Hazırlığı (İlk 25.000 veri + 5mA Filtresi)
subset_range_1 = 1:25000; 
subset_range_2 = 25001:54900; 
% Sızıntı verisini kes
I_leak_raw = leakageCurrent_filtered_A(subset_range_1);
Q_leak_raw = leakageFlow_filtered_lpm(subset_range_1);

% Basınç verisini kes
I_press_raw = pressureGainCurrent_A(subset_range_2);
P_press_raw = pressureGainPressure_bar(subset_range_2);


% +/- 5mA Filtreleme
limit_I = 5e-3; 

% Sızıntı verisi (Optimizasyon Vektörü için)
idx_leak = abs(I_leak_raw) <= limit_I;
I_leak_opt = I_leak_raw(idx_leak);
Q_leak_opt = Q_leak_raw(idx_leak) * lpm2m3Ps; % [m3/s]

% Basınç verisi (Optimizasyon Vektörü için)
idx_press = abs(I_press_raw) <= limit_I;
I_press_opt = I_press_raw(idx_press);
P_press_opt = P_press_raw(idx_press) * bar2Pa; % [Pa]

% Basınç verisini enterpolasyon için hazırla (Unique hale getir)
% x0 değiştiğinde PL'yi bu tablodan çekeceğiz
[I_press_unique, unique_idx] = unique(I_press_opt);
P_press_unique = P_press_opt(unique_idx);

% 2. Optimizasyon Ayarları (Trust-Region-Reflective)
% Başlangıç tahmini
x0_init = 0.0005; 
lb = 0.00001; % Alt sınır (0 olamaz)
ub = 0.005;   % Üst sınır (5mA)

options = optimoptions('lsqnonlin', ...
    'Algorithm', 'trust-region-reflective', ... % İstediğiniz algoritma
    'Display', 'iter-detailed', ...             % Detaylı gösterim
    'FunctionTolerance', 1e-8, ...
    'StepTolerance', 1e-8, ...
    'MaxFunctionEvaluations', 1000);

% 3. Amaç Fonksiyonunun Tanımlanması
% lsqnonlin bizden kareleri alınacak bir "Hata Vektörü" (Residuals) ister.
% Ağırlıklandırmayı vektör elemanlarını çarparak yapıyoruz.
objectiveFunc = @(x) calculateResiduals(x, ...
    I_leak_opt, Q_leak_opt, ...
    I_press_unique, P_press_unique, ...
    I_press_opt, P_press_opt, ...
    Ps, Pr, K, Q_leak_net_max, Q_tare);

% 4. Çözücü Çalıştırma
[x0_opt, resnorm, residual, exitflag, output] = lsqnonlin(objectiveFunc, x0_init, lb, ub, options);

fprintf('\n--- Sonuçlar ---\n');
fprintf('Optimize edilmiş x0: %.6f Amper (%.3f mA)\n', x0_opt, x0_opt*1e3);
fprintf('Kalan Hata Karesi Normu (Resnorm): %.4e\n', resnorm);

%% 5. Sonuçların Hesaplanması ve Görselleştirme
x0 = x0_opt; 

% Optimize x0 için final PL ve k değerleri
PL_at_x0 = interp1(I_press_unique, P_press_unique, x0, 'linear', 'extrap');
PL_at_x0 = abs(PL_at_x0); 

term = (Ps - Pr + PL_at_x0) / (Ps - Pr - PL_at_x0);
if term < 0, term = 0; end
k = (1/2) * sqrt(term) - 1;

fprintf('Final Model Parametreleri:\n   k = %.4f\n   PL(@x0) = %.2f Bar\n', k, PL_at_x0/bar2Pa);

% Çizim için simülasyon (Tam aralık)
I_range = linspace(-25e-3, 25e-3, 10e3); 
f_x = (1 + abs(I_range)/x0).^2 .* (1 + k * abs(I_range)/x0).^2;

% Ölçekleme
% Hedef_Tepe = Q_leak_net_max; 
% Hesaplanan_Tepe = (2 * K * sqrt(Ps - Pr) * x0);
% Skala = Hedef_Tepe / Hesaplanan_Tepe;

% Modeller
Qs_model = 2 * K * sqrt(Ps - Pr) * (x0 + abs(I_range)) .* (1 + f_x).^(-1/2); 
PL_model = ((f_x - 1) ./ (f_x + 1)) * (Ps - Pr) .* sign(I_range);

% Grafik
figure('Name', 'Nonlinear Least Squares Optimization (80/20 Weight)', 'Color', 'w');

% Grafik 1: Sızıntı
subplot(2,1,1);
plot(I_range*1e3, Qs_model*m3Ps2lpm + Q_tare*m3Ps2lpm, 'r-', 'LineWidth', 3); hold on; % Dönüşüm kontrolü
plot(I_leak_raw*1e3, Q_leak_raw, '.','Color',[0.6 0.6 0.6],'LineWidth',0.5); 
xline(5, '--k'); xline(-5, '--k');
title(['Sızıntı Debisi (Weight: %80, x0 = ' num2str(x0*1e3, '%.3f') ' mA)']);
xlabel('Akım (mA)'); ylabel('Debi (l/min)'); grid on; xlim([-25 25]);

% Grafik 2: Basınç
subplot(2,1,2);
plot(I_range*1e3, PL_model/bar2Pa, 'r-', 'LineWidth', 2); hold on;
plot(I_press_raw*1e3, P_press_raw, 'b.', 'MarkerSize', 4); 
xline(5, '--k'); xline(-5, '--k');
title(['Basınç Hassasiyeti (Weight: %20, k = ' num2str(k, '%.3f') ')']);
xlabel('Akım (mA)'); ylabel('Yük Basıncı (bar)'); grid on; xlim([-25 25]);

%% Yardımcı Fonksiyon: Residual Vector Calculation
function residuals = calculateResiduals(x0_try, I_leak, Q_leak_data, I_press_uni, P_press_uni, I_press_data, P_press_data, Ps, Pr, K, Q_net_max, Q_tare)
    
    % 1. Fiziksel Parametrelerin Hesabı
    % x0'a karşılık gelen PL'yi bul
    PL_try = interp1(I_press_uni, P_press_uni, x0_try, 'linear', 'extrap');
    PL_try = abs(PL_try); 
    
    % k katsayısı
    term = (Ps - Pr + PL_try)/(Ps - Pr - PL_try);
    if term < 0 || isnan(term), term = 0; end
    k_try = 0.5 * sqrt(term) - 1;
    
    % 2. Sızıntı Modeli Hesabı (Vektörel)
    f_x_leak = (1 + abs(I_leak)/x0_try).^2 .* (1 + k_try * abs(I_leak)/x0_try).^2;
    Calc_Peak = (2 * K * sqrt(Ps - Pr) * x0_try);
    Scale = Q_net_max / Calc_Peak;
    
    Q_model_net = Scale * (1.3/1.1) * 2 * K * sqrt(Ps - Pr) * (x0_try + abs(I_leak)) .* (1 + f_x_leak).^(-1/2);
    Q_model_total = Q_model_net + Q_tare;
    
    % 3. Basınç Modeli Hesabı (Vektörel)
    f_x_press = (1 + abs(I_press_data)/x0_try).^2 .* (1 + k_try * abs(I_press_data)/x0_try).^2;
    PL_model = ((f_x_press - 1) ./ (f_x_press + 1)) * (Ps - Pr) .* sign(I_press_data);
    
    % 4. Hata Vektörlerinin Oluşturulması ve Normalizasyon
    % Normalizasyon: Hataları aynı mertebeye (0-1 arasına) çekmek şarttır.
    % Q (Debi) hatası normalize: (Model - Data) / Max_Data
    Q_residuals = (Q_model_total - Q_leak_data) / max(Q_leak_data); 
    
    % P (Basınç) hatası normalize: (Model - Data) / Ps
    P_residuals = (PL_model - P_press_data) / Ps;
    
    % 5. Ağırlıklandırma
    % lsqnonlin fonksiyonu sum(residuals.^2) işlemini yapar.
    % Biz %80 Etki istiyorsak: 0.8 * J_leak + 0.2 * J_press olmalı.
    % J = sum( (w * r)^2 ) mantığıyla, Residual vektörünü karekök(Ağırlık) ile çarpmalıyız.
    
    % Not: Veri sayısı (N) farklıysa N'e bölmek de gerekebilir ama 
    % burada veri setlerini aynı filtreden geçirdiğimiz için N'ler yakındır.
    
    w_leak = sqrt(0.980);
    w_press = sqrt(0.020);
    
    % Tek bir uzun vektör olarak birleştirme
    residuals = [w_leak * Q_residuals; w_press * P_residuals];
    
    % NaN kontrolü (Algoritmanın çökmemesi için)
    residuals(isnan(residuals)) = 0;
end