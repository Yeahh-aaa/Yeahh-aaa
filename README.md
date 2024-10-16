% Model Gravitasi Sederhana dengan Lipatan atau Patahan
clear; clc; close all;

% Parameter model
nx = 100; % jumlah titik di sepanjang sumbu x
nz = 50;  % jumlah titik di sepanjang sumbu z
x = linspace(-5000, 5000, nx); % jarak horizontal (meter)
z = linspace(0, 5000, nz); % kedalaman (meter)
[X, Z] = meshgrid(x, z);

% Model densitas batuan (kg/m^3)
rho_bg = 2670; % densitas batuan latar belakang (umum untuk batuan granit)
rho_fault = 2800; % densitas batuan pada zona patahan (batuan lebih padat)
rho_fold = 2500; % densitas batuan pada zona lipatan (batuan lebih ringan)

% Buat model struktur
model = rho_bg * ones(nz, nx);

% Membuat patahan (Fault)
fault_pos = round(nx/2); % posisi patahan di tengah
fault_depth = round(nz/3); % kedalaman patahan
model(fault_depth:end, fault_pos:end) = rho_fault;

% Membuat lipatan (Fold)
fold_amplitude = 500; % amplitudo lipatan (meter)
fold_wavelength = 3000; % panjang gelombang lipatan (meter)
fold_depth = round(nz/2); % kedalaman mulai lipatan
for i = 1:nx
    fold_height = fold_amplitude * sin(2*pi*(x(i))/fold_wavelength);
    fold_height_idx = find(z >= fold_depth - fold_height, 1, 'first');
    if ~isempty(fold_height_idx)
        model(fold_height_idx:end, i) = rho_fold;
    end
end

% Plot model densitas
figure;
imagesc(x, z, model);
colorbar;
title('Model Struktur Gravitasi (Lipatan dan Patahan)');
xlabel('Jarak (m)');
ylabel('Kedalaman (m)');
set(gca, 'YDir', 'reverse');

% Hitung anomali gravitasi (menggunakan metode forward modeling sederhana)
G = 6.67430e-11; % konstanta gravitasi (m^3 kg^-1 s^-2)
anomali_gravitasi = zeros(1, nx); % inisialisasi vektor anomali gravitasi

% Perhitungan anomali gravitasi dengan integral (2D forward modeling)
dz = z(2) - z(1); % resolusi vertikal (meter)
dx = x(2) - x(1); % resolusi horizontal (meter)

for ix = 1:nx
    for iz = 1:nz
        r = sqrt(x(ix)^2 + z(iz)^2); % jarak dari pengukuran ke elemen densitas
        if r ~= 0
            anomali_gravitasi(ix) = anomali_gravitasi(ix) + ...
                2 * G * (model(iz, ix) - rho_bg) * dz / r; % sum over z
        end
    end
end

% Plot anomali gravitasi
figure;
plot(x, anomali_gravitasi * 1e8); % konversi ke mGal
title('Anomali Gravitasi (mGal)');
xlabel('Jarak (m)');
ylabel('Anomali Gravitasi (mGal)');
grid on;
