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

# ==========================================================
# 4. CALIBRACION DE UMBRALES Y LOGICA DE DETECCION BINARIA
# ==========================================================

# Funcion para calibrar umbral y aplicar logica de deteccion con debouncing
def calibrate_threshold(envelope, fs, rest_time=2.0):
    # Calibracion del umbral utilizando los primeros 2 segundos de la senal (reposo)
    baseline_segment = env[time < 2.0]
    rest_level = np.mean(baseline_segment)
    max_voluntary_contraction = np.max(env)
    threshold = rest_level + 0.2 * (max_voluntary_contraction - rest_level)

    # Logica de deteccion binaria con debouncing
    raw_detection = (env > threshold).astype(int)
    clean_detection = np.zeros_like(raw_detection)

    # Parametro de estabilidad (debouncing)
    hold_samples = int(0.15 * fs)       # 150ms de tolerancia
    counter = 0
    
    for i in range(len(raw_detection)):
        if raw_detection[i] == 1:
            clean_detection[i] = 1
            counter = hold_samples      # Reiniciamos el contador de deteccion
        else:
            if counter > 0:
                clean_detection[i] = 1  # Mantenemos la deteccion activa durante el periodo de hold
                counter -= 1
            else:
                clean_detection[i] = 0  # No hay deteccion activa
    return threshold, raw_detection, clean_detection

# ================================================
# 1. CARGA Y VISUALIZACION DE LA SENAL EMG CRUDA
# ================================================

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

# ===============================================
# 2. VALIDACION DE LA CALIDAD DE LA SENAL Y FFT
# ===============================================

# Calcular FFT y espectro de frecuencia 
yf = fft(signal_centered)
xf = fftfreq(n, 1/fs)[:n//2]
magnitude = 2.0/n * np.abs(yf[0:n//2])

plt.figure(figsize=(12, 6))

# Grafica de la senal en el tiempo (para ver calidad de cruda)
plt.subplot(2, 1, 1)
plt.plot(time, signal, color='tab:blue')
plt.title("Señal EMG cruda")
plt.xlabel("Tiempo (s)")
plt.ylabel("Amplitud (mV)")
plt.grid(True)

# Grafica del Espectro de Potencia (FFT)
plt.subplot(2, 1, 2)
plt.plot(xf, magnitude, color='tab:red')
plt.title("Espectro de frecuencia (FFT)")
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
plt.title("Filtrado digital")
plt.xlabel("Tiempo (s)")
plt.ylabel("Amplitud (mV)")
plt.grid(True)
plt.legend()

# Grafica de la senal rectificada (valor absoluto)
plt.subplot(4, 1, 2)
plt.plot(time, np.abs(sig_f), color='orange', lw=0.5)
plt.title("Senal rectificada (valor absoluto)")
plt.xlabel("Tiempo (s)")
plt.ylabel("Amplitud (mV)")
plt.grid(True)

# Grafica de la envolvente RMS (ventana muscular)
plt.subplot(4, 1, 3)
plt.plot(time, env, color='red', lw=2)
plt.fill_between(time, env, color='red', alpha=0.2)
plt.title("Envolvente RMS (ventana muscular)")
plt.xlabel("Tiempo (s)")
plt.ylabel("Amplitud (mV)")
plt.grid(True)

# Validacion FFT de la senal filtrada 
yf_f = fft(sig_f)
mag_f = 2.0/n * np.abs(yf_f[0:n//2])

plt.subplot(4, 1, 4)
plt.plot(xf, mag_f, color='purple')
plt.title("Espectro de la senal filtrada")
plt.xlabel("Frecuencia (Hz)")
plt.xlim(0, 500)
plt.grid(True)

plt.tight_layout()
plt.show()

# =====================================================================
# 4. SIMULACION DE CALIBRACION DINAMICA Y LOGICA DE DETECCION BINARIA
# =====================================================================

threshold, raw_detection, clean_detection = calibrate_threshold(env, fs)

# Grafica de la envolvente con el umbral
plt.figure(figsize=(12, 6))
plt.subplot(2, 1, 1)
plt.plot(time, env, color='red', label='Envolvente RMS')
plt.axhline(y=threshold, color='green', linestyle='--', label='Umbral')
plt.title("Envolvente con Umbral")
plt.legend()

# Grafica de la deteccion binaria (con y sin debouncing)
plt.subplot(2, 1, 2)
plt.plot(time, raw_detection, color='gray', alpha=0.5, label='Detección con ruido (Original)')
plt.plot(time, clean_detection, color='darkgreen', lw=2, label='Detección Limpia (Control)')
plt.fill_between(time, 0, clean_detection, color='green', alpha=0.2)
plt.title("Salida de Control Filtrada")
plt.ylim(-0.1, 1.2)
plt.legend()

plt.tight_layout()
plt.show()