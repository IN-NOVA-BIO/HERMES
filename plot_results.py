import pandas as pd
import matplotlib.pyplot as plt

# Cargar los datos desde el CSV generado por HERMES (con codificacion utf-16 para evitar problemas de caracteres)
try:
    df = pd.read_csv('resultados.csv', encoding='utf-16')
except:
    # Si falla, intenta la normal
    df = pd.read_csv('resultados.csv')

# Configuracion de la grafica
fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(12, 8), sharex=True)
fig.suptitle('Analisis de Senal EMG - Proyecto HERMES', fontsize=16)

# Subplot 1. Senal cruda vs filtrada
ax1.plot(df['time_s'], df['raw_signal'], label='Cruda', alpha=0.5, color='gray')
ax1.plot(df['time_s'], df['filtered_signal'], label='Filtrada (Notch + BP)', color='blue')
ax1.set_ylabel('mV')
ax1.grid(True, alpha=0.3)

# Subplot 2. Envolvente RMS y umbral
ax2.plot(df['time_s'], df['envelope'], label='Envolvente RMS', color='orange', linewidth=2)
ax2.axhline(y=df['envelope'].mean() * 1.2, color='red', linestyle='--', label='Umbral estatico') 
ax2.set_ylabel('RMS')
ax2.legend(loc='upper right')
ax2.grid(True, alpha=0.3)

# Subplot 3. Deteccion binaria (salida para el teclado)
ax3.fill_between(df['time_s'], 0, df['detection'], color='green', alpha=0.3, label='Estado activado')
ax3.plot(df['time_s'], df['detection'], color='green', linewidth=1)
ax3.set_ylabel('Deteccion (0/1)')
ax3.set_xlabel('Tiempo (s)')
ax3.set_ylim(-0.1, 1.1)
ax3.legend(loc='upper right')
ax3.grid(True, alpha=0.3)

plt.tight_layout()
plt.show()