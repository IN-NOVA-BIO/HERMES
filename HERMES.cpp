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
static constexpr float THRESHOLD_SCALE = 0.25f; // 25% del rango MVC-reposo

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
        for (int i = 0; i < numStages; ++i) 
            y = stages[i].process(y);       // Procesa la salida de cada etapa secuencialmente
        return y;   
    }

    // Resetea los estados de todas las etapas de biquad en la cascada
    void reset() {
        for (int i = 0; i < numStages; i++)
            stages[i].reset();
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

    // Actualiza el buffer circular con una nueva muestra y devuelve el valor RMS actualizado
    float update(float sample) {
        float newSq = sample * sample;          // Calcula el cuadrado de la nueva muestra
        sumSq -= buf[head];                     // Resta el valor antiguo del cuadrado que se va a eliminar de la ventana
        sumSq += newSq;                         // Agrega el nuevo valor del cuadrado a la suma total
        buf[head] = newSq;                      // Almacena el nuevo valor del cuadrado en el buffer
        head = (head + 1) % winSize;            // Avanza el indice del head del buffer circular
        return std::sqrt(sumSq / winSize);      // Devuelve la raiz cuadrada de la media de los cuadrados
    }

    // Resetea el buffer circular y el acumulador de la suma de los cuadrados
    void reset() {
        head = 0;                                       // Resetea el indice del head del buffer circular
        sumSq = 0.0;                                    // Resetea el acumulador de la suma de los cuadrados
        std::memset(buf, 0, sizeof(float) * winSize);   // Limpia el buffer circular para la ventana de RMS
    }

};

// ESTADO GLOBAL DEL PIPELINE
struct EMGPipeline {
    Biquad notchFilter;             // Filtro notch 60 Hz
    BiquadCascade bandpassFilter;   // Bandpass 20-450 Hz
    RmsBuffer rmsBuffer;            // Envolvente RMS

    // Estado de calibracion 
    bool calibrated = false;        // Calibracion del sistema
    float threshold = 0.0f;          // Umbral de deteccion binaria 
    float restLevel = 0.0f;         // Nivel de reposo 
    float maxMVC = 0.0f;            // Nivel de contraccion maxima voluntaria

    // Estado del debouncing
    int holdCounter = 0;            // Muestras restantes del hold activo
    int holdSamples = 0;            // Duracion del hold time en muestras

    // Contadores para calibracion offline (primeros 2 s)
    double calibSum = 0.0;          // Acumulador para la suma de las muestras
    int calibCount = 0;             // Contador de muestras para la calibracion
    int calibSamples = 0;           // Total de muestras en el periodo

    // Inicializa el pipeline de procesamiento de EMG con la fs 
    void init(float fs) {
        // Diseno de filtros segun fs
        notchFilter = designNotch(fs, NOTCH_FREQ, NOTCH_Q);                 // Filtro notch para eliminar ruido de 60 Hz
        designButterworthBandpass(BAND_LOW, BAND_HIGH, fs, bandpassFilter); // Filtro bandpass para limitar el rango

        // Inicializacion del buffer de RMS con el tamaño de la ventana en muestras
        int winSamples = (int)(RMS_WINDOW_SIZE * fs);       // Calcula el tamano de la ventana en muestras
        rmsBuffer.init(winSamples);                         // Inicializa el buffer con el tamano calculado  

        // Calcula los parametros de calibracion y debouncing en muestras
        holdSamples = (int)(HOLD_TIME_SEC * fs);            // Calcula la duracion del hold time en muestras
        calibSamples = (int)(REST_CALIB_SEC * fs);          // Calcula el numero de muestras para la calibracion

        // Resetea el estado de calibracion y debouncing
        calibrated = false;         // Resetea el estado de calibracion
        calibSum = 0.0;             // Resetea el acumulador de la suma para la calibracion
        calibCount = 0;             // Resetea el contador de muestras para la calibracion
        holdCounter = 0;            // Resetea el contador de hold time para el debouncing
        maxMVC = 0.0f;              // Resetea el nivel de contraccion maxima voluntaria
    }

    // Resetea el estado de todo el pipeline
    void softReset() {
        notchFilter.reset();        
        bandpassFilter.reset();     
        rmsBuffer.reset();          
        calibrated = false;         
        calibSum = 0.0;             
        calibCount = 0;      
        maxMVC = 0.0f;
        holdCounter = 0;       
    }
};

// PROCESAMIENTO DE UNA MUESTRA (PIPELINE COMPLETO)
int processSample(EMGPipeline& p, float rawSample, float dcOffset, float& outFiltered, float& outEnvelope) {
    // Paso 1. Eliminacion de offset DC
    float x = rawSample - dcOffset;                 // Resta el offset DC de la muestra cruda

    // Paso 2. Filtrado Notch 60 Hz
    float xNotch = p.notchFilter.process(x);        // Procesa la muestra a traves del filtro notch

    // Paso 3. Filtrado Bandpass 20-450 Hz
    float xBand = p.bandpassFilter.process(xNotch); // Procesa la muestra a traves del filtro bandpass
    outFiltered = xBand;                            // Salida del filtro bandpass 

    // Paso 4. Calculo de la envolvente RMS
    float env = p.rmsBuffer.update(xBand);          // Actualiza el buffer de RMS con la muestra filtrada
    outEnvelope = env;                              // Salida de la envolvente RMS
    
    // Paso 5. Calibracion offline (primeros 2 segundos)
    if (!p.calibrated) {
        p.calibSum += env;
        p.calibCount++;

        // Al acumular suficientes muestras para la calibracion, calcula el nivel de reposo y el umbral de deteccion   
        if(p.calibCount >= p.calibSamples) {
            // Calcula el nivel de reposo como el promedio de las muestras durante el periodo de calibracion
            p.restLevel = (float)(p.calibSum / p.calibCount);
            // Establece el umbral de deteccion como un porcentaje por encima del nivel de reposo 
            p.threshold = p.restLevel * (1.0f + THRESHOLD_SCALE); 
            p.calibrated = true;            
        }
        return -1; // Indica que el sistema aun no esta calibrado   
    }
    
    // Paso 6. Deteccion binaria con debouncing
    int rawDetection = (env > p.threshold) ? 1 : 0; // Deteccion binaria sin debouncing
    int cleanDetection;                             // Deteccion binaria con debouncing

    if(rawDetection == 1) {
        cleanDetection = 1;                 // Salida de deteccion limpia (contraccion detectada)
        p.holdCounter = p.holdSamples;      // Reinicia el contador de hold time
        if(env > p.maxMVC) p.maxMVC = env;  // Actualiza el nivel de CMV si se detecta un nuevo maximo
    } else if(p.holdCounter > 0) {
        cleanDetection = 1;                 // Mantiene la deteccion activa durante el hold time
        p.holdCounter--;                    // Decrementa el contador de hold time
    } else {
        cleanDetection = 0;                 // Salida de deteccion limpia (reposo detectado)
    }

    return cleanDetection;                  // Devuelve la deteccion binaria con debouncing
}

// CARGA DE DATOS DESDE ARCHIVO (MODO OFFLINE)
struct SignalData {
    std::vector<float> time;        // Vector para almacenar los tiempos de cada muestra
    std::vector<float> signal;      // Vector para almacenar las muestras de la señal EMG
};

bool loadSignalFromFile(const std::string& path, SignalData& data) {
    std::ifstream file(path); // Abre el archivo para lectura
    if (!file.is_open()) {
        std::cerr << "[ERROR] No se pudo abrir el archivo " << path << std::endl;
        return false; 
    }

    std::string line;   // Variable para almacenar cada linea leida del archivo
    // Lee el archivo linea por linea y extrae el tiempo y la señal de cada linea
    while(std::getline(file, line)) { 
        std::istringstream iss(line);
        float t, s;
        if (!(iss >> t >> s)) continue;                 // Ignora lineas con formato incorrecto
        if (std::isnan(t) || std::isnan(s)) continue;   // Ignora filas con valores no numericos (NaN)
        data.time.push_back(t);                         // Almacena el tiempo en el vector
        data.signal.push_back(s);                       // Almacena la señal en el vector
    }

    // Verifica que se hayan cargado datos correctamente
    return !data.time.empty();                          
}