# Procesamiento de base de datos de senales EMG
# Este codigo unicamente se usa como una validacion controlada del procesamiento de las senales EMG
# EL codigo sirve como borrador, tras la validacion de la base de datos se transportara a C/C++ para su uso en la ESP32
# La base de datos fue obtenida de Physionet: https://physionet.org/content/emgdb/1.0.0/

import numpy as np
import matplotlib.pyplot as plt
import os

base_path = os.path.dirname(__file__)                       
file_path = os.path.join(base_path, "examples-of-electromyograms-1.0.0","emg_healthy.txt")

data = np.genfromtxt(file_path, encoding='latin-1', invalid_raise=False)

# Separar columnas
x = data[:, 0]
y = data[:, 1]

# Graficar
plt.plot(x, y)
plt.xlabel('Time (sec)')
plt.ylabel('Amplitude (mV)')
plt.title('EMG Healthy')
plt.xlim(0, 1.5)
plt.show()