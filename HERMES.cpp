#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cstring>
#include <vector>       
#include <algorithm>    
#include <numeric>      

#define M_PI 3.14159265358979323846

// DECLARACION DE PARAMETROS DEL SISTEMA
// Parametro de muestreo
static float FS = 0.0f;

// Parametros de filtrado
static constexpr float NOTCH_FREQ = 60.0f;      // Hz - frecuencia de la red electrica
static constexpr float NOTCH_Q = 4.0f;          // Factor de calidad del notch
static constexpr float BAND_LOW = 20.0f;        // Hz - limite inferior EMG
static constexpr float BAND_HIGH = 450.0f;      // Hz - limite superior EMG
static constexpr int BUTTER_ORDER = 4;          // Orden del filtro Butterworth

// Parametro de la ventana de segmentacion
static constexpr float RMS_WINDOW_SIZE = 0.2f;  // 200 ms

// Parametros de deteccion binaria (debouncing)
// Se pueden ajustar para mejorar la deteccion de contracciones musculares y la calibracion propuesta
static constexpr float HOLD_TIME_SEC = 0.15f;   // 150 ms de hold time
static constexpr float REST_CALIB_SEC = 2.0f;   // 2 segundos para calibrar 
static constexpr float THRESHOLD_SCALE = 0.2f;  // 20% del rango MVC-reposo

// ESTRUCTURAS DE FILTROS IIR DE SEGUNDO ORDEN (BIQUAD)
struct Biquad{
    // Coeficientes del filtro
    float b0, b1, b2; // Coeficientes de la parte de alimentacion (feedforward)
    float a1, a2;     // Coeficientes de la parte de retroalimentacion (feedback)

    // Estados del filtro (memoria de las muestras anteriores)
    float w1 = 0.0f, w2 = 0.0f; 

    // Procesa una muestra de entrada y devuelve la salida filtrada
    // OPTIMIZACION: inline para reducir la sobrecarga de la llamada a la función
    inline float process(float x) {
        // Calculo de la salida del filtro usando la estructura de biquad
        float y = b0 * x + w1;      // Salida actual basada en la entrada y el estado anterior
        w1 = b1 * x - a1 * y + w2;  // Actualiza w1 para la proxima muestra
        w2 = b2 * x - a2 * y;       // Actualiza w2 para la proxima muestra
        return y;
    }

    // Resetea los estados del filtro 
    void reset() { w1 = 0.0f; w2 = 0.0f; }
};

// CASCADA DE BIQUADS PARA FILTROS DE ORDEN SUPERIOR
struct BiquadCascade {
    static constexpr int MAX_STAGES = 2;    // Orden 4 -> 2 biquads
    Biquad stages[MAX_STAGES];              // Array de biquads para la cascada
    int numStages = 0;                      // Numero de etapas en la cascada

    // Procesa una muestra de entrada a traves de todas las etapas de biquad en la cascada
    // Devuelve la salida filtrada despues de pasar por cada etapa secuencialmente
    float process(float x) {
        float y = x;
        for (int i = 0; i < numStages; ++i) {
            y = stages[i].process(y);       // Procesa la salida de cada etapa secuencialmente
        }
        return y;   
    }

    // Resetea los estados de todas las etapas de biquad en la cascada
    void reset() {
        for (int i = 0; i < numStages; i++) {
            stages[i].reset();
        }
    }
};

// DISENO DE FILTRO NOTCH
Biquad designNotch(float fs, float freq, float Q) {
    float w0 = 2.0f * M_PI * freq / fs;     // Frecuencia angular normalizada
    float alpha = sin(w0) / (2.0f * Q);     // Calculo de los coeficientes del filtro notch usando la formula estandar
    float cos_w0 = cos(w0);                 // Calculo del coseno de w0 para los coeficientes
    float a0 = 1.0f + alpha;                // Coeficiente de normalización

    Biquad bq;
    bq.b0 = 1.0f / a0;                      
    bq.b1 = -2.0f * cos_w0 / a0;            
    bq.b2 = 1.0f / a0;                     
    bq.a1 = -2.0f * cos_w0 / a0;            
    bq.a2 = (1.0f - alpha) / a0;            
    return bq;
}

// DISENO DE FILTRO BUTTERWORTH BANDPASS DE ORDEN 4 (2 SECCIONES BIQUAD)
void designButterworthBandpass(float lowHz, float highHz, float fs, BiquadCascade& cascade) {
    cascade.numStages = 2;      // Orden 4 -> 2 biquads
    
    // Seccion 1. Paso alto de 2 orden en lowHz
    // Transformacion bilineal para convertir el filtro analogico a digital
    {
        float wc = 2.0f * M_PI * lowHz / fs;    // Frecuencia angular normalizada para lowHz
        float k = std::tan(wc / 2.0f);          // Calculo de k para la transformacion bilineal
        float k2 = k * k;                       // k al cuadrado para los coeficientes
        float sqrt2 = std::sqrt(2.0f);          // Raiz cuadrada de 2 para los coeficientes de Butterworth
        float a0 = 1.0f + sqrt2 * k + k2;       // Coeficiente de normalizacion para la primera seccion

        // Coeficientes para la primera seccion biquad (lowHz)
        Biquad& bq = cascade.stages[0];
        // Paso alto: el numerador tiene +1, -2, +1 
        bq.b0 = 1.0f / a0;      
        bq.b1 = -2.0f / a0;     
        bq.b2 = 1.0f / a0;      
        bq.a1 = (2.0f * (k2 - 1.0f)) / a0;   
        bq.a2 = (1.0f - sqrt2 * k + k2) / a0;
    }

    // Seccion 2. Paso bajo de 2 orden en highHz
    // Transformacion bilineal para convertir el filtro analogico a digital
    {
        float wc = 2.0f * M_PI * highHz / fs;   // Frecuencia angular normalizada para highHz
        float k = std::tan(wc / 2.0f);          // Calculo de k para la transformacion bilineal
        float k2 = k * k;                       // k al cuadrado para los coeficientes
        float sqrt2 = std::sqrt(2.0f);          // Raiz cuadrada de 2 para los coeficientes de Butterworth
        float a0 = 1.0f + sqrt2 * k + k2;       // Coeficiente de normalizacion para la segunda seccion

        // Coeficientes para la segunda seccion biquad (highHz)
        Biquad& bq = cascade.stages[1];
        // Paso bajo: el numerador tiene k2, 2*k2, k2
        bq.b0 = k2 / a0;
        bq.b1 = 2.0f * k2 / a0;
        bq.b2 = k2 / a0;
        bq.a1 = (2.0f * (k2 - 1.0f)) / a0;
        bq.a2 = (1.0f - sqrt2 * k + k2) / a0;
    }
}

// BUFFER CIRCULAR PARA LA ENVOLVENTE RMS
struct RmsBuffer {
    static constexpr int MAX_SIZE = 2048;   // Maximo tamano del buffer para la ventana de RMS 
    float buf[MAX_SIZE] = {};               // Buffer circular para almacenar las muestras de la ventana de RMS
    int head = 0;                           // Indice del head del buffer circular
    int winSize = 0;                        // Tamaño de la ventana de RMS en muestras
    double sumSq = 0.0;                     // Acumulador de la suma de los cuadrados

    // Inicializa el buffer circular para la ventana de RMS con el tamaño especificado
    void init(int WindowSamples) {
        winSize = WindowSamples;                        // Establece el tamaño de la ventana de RMS en muestras
        head = 0;                                       // Resetea el indice del head del buffer circular
        sumSq = 0.0;                                    // Resetea el acumulador de la suma de los cuadrados
        std::memset(buf, 0, sizeof(float) * winSize);   // Limpia el buffer circular para la ventana de RMS
    }

    float update(float sample) {
        float newSq = sample * sample;          // Calcula el cuadrado de la nueva muestra
        sumSq -= buf[head];                     // Resta el valor antiguo del cuadrado que se va a eliminar de la ventana
        sumSq += newSq;                         // Agrega el nuevo valor del cuadrado a la suma total
        buf[head] = newSq;                      // Almacena el nuevo valor del cuadrado en el buffer
        head = (head + 1) % winSize;            // Avanza el indice del head del buffer circular
        return std::sqrt(sumSq / winSize);      // Devuelve la raiz cuadrada de la media de los cuadrados
    }
};