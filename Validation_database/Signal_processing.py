# Procesamiento de base de datos de senales EMG
# Este codigo unicamente se usa como una validacion controlada del procesamiento de las senales EMG
# EL codigo sirve como borrador, tras la validacion de la base de datos se transportara a C/C++ para su uso en la ESP32
# La base de datos fue obtenida de Physionet: https://physionet.org/content/emgdb/1.0.0/

import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft, fftfreq
from scipy.signal import butter, lfilter, iirnotch
import os

# ================================================
# 3. FILTRO DIGITAL Y LOGICA DE VENTANA MUSCULAR
# ================================================

# Funcion de procesamiento de senal EMG 
def emg_processing_pipeline(raw_signal, fs, low=20, high=450, notch_freq=60):
    nyq = 0.5 * fs
    
    # Filtro Notch para eliminar ruido de 60 Hz
    b_notch, a_notch = iirnotch(notch_freq / (0.5 * fs), 4)
    sig_notch = lfilter(b_notch, a_notch, raw_signal)

    # Filtrado de banda (20-450 Hz)
    lowcut = low / nyq
    highcut = high / nyq
    b_band, a_band = butter(4, [lowcut, highcut], btype='band')
    sig_filtered = lfilter(b_band, a_band, sig_notch)

    # Envolvente RMS (ventana de 200 ms)
    window_samples = int(0.2 * fs)
    envelope = np.sqrt(np.convolve(sig_filtered**2, np.ones(window_samples)/window_samples, mode='same'))
    
    return sig_filtered, envelope

# ====================================================
# 1. CARGA Y VISUALIZACION DE LA SENAL EMG CRUDA
# ====================================================

# Cargar y procesar la senal EMG desde la base de datos
base_path = os.path.dirname(__file__)                       
file_path = os.path.join(base_path, "examples-of-electromyograms-1.0.0","emg_neuropathy.txt")

# Cargar datos desde el archivo de texto, ignorando filas con datos faltantes o no numericos
data = np.genfromtxt(file_path, encoding='latin-1', invalid_raise=False)
data = data[~np.isnan(data).any(axis=1)]

# Separar columnas
time = data[:, 0]
signal = data[:, 1]

# Remover componente DC (centrar la senal)
signal_centered = signal - np.mean(signal)      

# Calcular frecuencia de muestreo (fs) 
fs = 1 / np.mean(np.diff(time))
n = len(signal)

# ==============================================
# 2. VALIDACION DE LA CALIDAD DE LA SENAL Y FFT
# ==============================================

# Calcular FFT y espectro de frecuencia 
yf = fft(signal_centered)
xf = fftfreq(n, 1/fs)[:n//2]
magnitude = 2.0/n * np.abs(yf[0:n//2])

# Visualizacion
plt.figure(figsize=(12, 6))

# Grafica de la senal en el tiempo (para ver calidad de cruda)
plt.subplot(2, 1, 1)
plt.plot(time, signal, color='tab:blue')
plt.title("Señal EMG Cruda (Dominio del Tiempo)")
plt.xlabel("Tiempo (s)")
plt.ylabel("Amplitud (mV)")
plt.grid(True)

# Grafica del Espectro de Potencia (FFT)
plt.subplot(2, 1, 2)
plt.plot(xf, magnitude, color='tab:red')
plt.title("Espectro de Frecuencia (FFT)")
plt.xlabel("Frecuencia (Hz)")
plt.ylabel("Magnitud")
plt.xlim(0, 500)  # El EMG util llega hasta 500Hz
plt.grid(True)

plt.tight_layout()
plt.show()

# =============================================================
# 3. PROCESAMIENTO EMG: FILTRADO, ENVOLVENTE Y VALIDACION FFT
# =============================================================

# Procesamiento de la senal EMG utilizando el pipeline definido
sig_f, env = emg_processing_pipeline(signal_centered, fs)

# Visualizacion de la senal procesada   
plt.figure(figsize=(12, 8))

# Grafica de la senal filtrada
plt.subplot(4, 1, 1)
plt.plot(time, signal_centered, color='lightgray', lw=0.5, label='Cruda')
plt.plot(time, sig_f, color='tab:blue', lw=0.8, label='Filtrada (20-450Hz)')
plt.title("Filtrado Digital")
plt.xlabel("Tiempo (s)")
plt.ylabel("Amplitud (mV)")
plt.grid(True)
plt.legend()

# Grafica de la senal rectificada (valor absoluto)
plt.subplot(4, 1, 2)
plt.plot(time, np.abs(sig_f), color='orange', lw=0.5)
plt.title("Senal Rectificada (Valor Absoluto)")
plt.xlabel("Tiempo (s)")
plt.ylabel("Amplitud (mV)")
plt.grid(True)

# Grafica de la envolvente RMS (ventana muscular)
plt.subplot(4, 1, 3)
plt.plot(time, env, color='red', lw=2)
plt.fill_between(time, env, color='red', alpha=0.2)
plt.title("Envolvente RMS (Ventana Muscular)")
plt.xlabel("Tiempo (s)")
plt.ylabel("Amplitud (mV)")
plt.grid(True)

# Validacion FFT de la senal filtrada 
yf_f = fft(sig_f)
mag_f = 2.0/n * np.abs(yf_f[0:n//2])

plt.subplot(4, 1, 4)
plt.plot(xf, mag_f, color='purple')
plt.title("Validacion: Espectro de la Senal Filtrada")
plt.xlabel("Frecuencia (Hz)")
plt.xlim(0, 500)
plt.grid(True)

plt.tight_layout()
plt.show()

# =====================================================================
# 4. SIMULACION DE CALIBRACION DINAMICA Y LOGICA DE DETECCION BINARIA
# =====================================================================

# Simulacion de calibracion: reposos (primeros 2s) y CVM (pico maximo de la senal) 
baseline_segment = env[time < 2.0]
rest_level = np.mean(baseline_segment)
max_voluntary_contraction = np.max(env)

# Umbral dinamico: Baseline + 20% del rango entre util
threshold = rest_level + 0.2 * (max_voluntary_contraction - rest_level)

# Logica de deteccion de activacion binaria
binary_detection = (env > threshold).astype(int)

# Visualizacion de la deteccion binaria
plt.figure(figsize=(12, 6))

# Envolvente y umbral de calibracion
plt.subplot(2, 1, 1)
plt.plot(time, env, color='red', lw=1.5, label='Envolvente RMS')
plt.fill_between(time, env, color='red', alpha=0.2)
plt.axhline(y=threshold, color='green', linestyle='--', label='Umbral Calibrado')
plt.axhline(y=rest_level, color='black', linestyle=':', label='Baseline (Reposo)')
plt.title("Envolvente y Calibracion")
plt.ylabel("Amplitud (mV)")
plt.legend()

# Salida de control (binaria)
plt.subplot(2, 1, 2)
plt.fill_between(time, 0, binary_detection, color='green', alpha=0.3)
plt.plot(time, binary_detection, color='darkgreen', label='Detección (1/0)')
plt.title("Salida de Control (Logica Binaria)")
plt.xlabel("Tiempo (s)")
plt.ylim(-0.1, 1.2)
plt.legend()

plt.tight_layout()
plt.show()