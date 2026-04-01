# Procesamiento de base de datos de senales EMG
# Este codigo unicamente se usa como una validacion controlada del procesamiento de las senales EMG
# EL codigo sirve como borrador, tras la validacion de la base de datos se transportara a C/C++ para su uso en la ESP32
# La base de datos fue obtenida de Physionet: https://physionet.org/content/emgdb/1.0.0/

import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft, fftfreq
from scipy.signal import butter, filtfilt, iirnotch
import os

base_path = os.path.dirname(__file__)                       
file_path = os.path.join(base_path, "examples-of-electromyograms-1.0.0","emg_healthy.txt")

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
plt.xlim(0, 500)  # El EMG útil llega hasta 500Hz
plt.grid(True)

plt.tight_layout()
plt.show()