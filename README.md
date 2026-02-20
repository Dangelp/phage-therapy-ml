# 🦠 Predicción de Susceptibilidad Fágica mediante Machine Learning

![Python](https://img.shields.io/badge/Python-3.8%2B-blue?style=flat-square&logo=python)
![Scikit-Learn](https://img.shields.io/badge/Scikit--Learn-Machine%20Learning-orange?style=flat-square&logo=scikit-learn)
![Biopython](https://img.shields.io/badge/Biopython-Genómica-green?style=flat-square)
![Status](https://img.shields.io/badge/Estado-Completado-success?style=flat-square)

> **Herramienta bioinformática in-silico para predecir la susceptibilidad de cepas clínicas de *Escherichia coli* a bacteriófagos utilizando firmas genómicas libres de alineamiento (K-mers) y algoritmos de Random Forest Multi-Output.**

## 📋 Resumen del Proyecto
La Resistencia a los Antimicrobianos (RAM) es una de las mayores amenazas para la salud pública global. La **fagoterapia** (uso de virus bacteriófagos) emerge como una alternativa terapéutica de alta precisión. Sin embargo, el emparejamiento entre una bacteria clínica y el fago correcto depende actualmente de cribados fenotípicos in-vitro (spot tests) que son lentos, costosos y difíciles de escalar.

Este proyecto propone una **solución computacional** para predecir el éxito terapéutico en minutos. A partir de secuencias genómicas bacterianas crudas (`.fna`), el pipeline extrae firmas genómicas y entrena un modelo capaz de recetar un cóctel óptimo de bacteriófagos, optimizando la toma de decisiones clínicas.

## 🔬 Metodología

El flujo de trabajo bioinformático consta de tres pilares fundamentales:

1. **Extracción de Características (Alignment-Free):** Procesamiento de 912 genomas completos de *E. coli* para extraer frecuencias de tetranucleótidos (k-mers, $k=4$). Esto genera un vector matemático de 256 dimensiones que captura la huella filogenética y el uso de codones sin necesidad de alineamientos computacionalmente costosos.
2. **Selección del Cóctel Terapéutico:** Implementación de un **Algoritmo Greedy** (Set Cover Problem) para identificar el "Cóctel de Oro" de 20 fagos que maximiza la cobertura de lisis sobre la población bacteriana.
3. **Aprendizaje Automático (Machine Learning):** Entrenamiento de un modelo `RandomForestClassifier` acoplado a un `MultiOutputClassifier`. El sistema predice simultáneamente la eficacia de los 20 fagos del cóctel frente a una nueva firma genómica de entrada.

## 📊 Resultados Clave

* **Alta Capacidad Predictiva:** El modelo alcanzó una **Exactitud Global (Global Accuracy) del 82.32%** en el conjunto de prueba independiente.
* **Interpretabilidad Biológica (Biomarcadores):** El análisis de importancia de características (Feature Importance) reveló que el modelo utiliza motivos genéticos específicos (ej. `TACG`, indicativo de sitios CpG) para predecir la resistencia. Esto sugiere que el algoritmo detecta indirectamente la presencia de **Sistemas de Restricción-Modificación (R-M)** bacterianos.

## 🚀 Estructura del Repositorio

```text
phage-therapy-ml/
├── data/
│   ├── raw/                 # Muestra de genomas crudos (.fna)
│   └── processed/           # Matrices generadas (X_Kmers_Matrix.csv)
├── scripts/
│   └── TFM_Modelo_Predictivo_Fagos.py  # Pipeline principal documentado
├── results/
│   └── figures/             # Gráficas de rendimiento y feature importance
├── requirements.txt         # Dependencias del entorno de Python
└── README.md                # Documentación del proyecto
