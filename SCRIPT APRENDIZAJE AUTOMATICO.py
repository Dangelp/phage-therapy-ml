#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 18 19:59:44 2026

@author: david
"""


#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
===============================================================================
Modelo de Aprendizaje Automático para Prescripción de Fagoterapia
AUTOR: David Angel Perez
FECHA: Febrero 2026
DESCRIPCIÓN: 
    Este script implementa un pipeline bioinformático completo que:
    1. Procesa genomas bacterianos crudos (.fna) extrayendo firmas de K-mers (k=4).
    2. Cruza la información genómica con una base de datos de interacciones fago-bacteria.
    3. Selecciona un 'Cocktail Terapéutico' óptimo de 20 fagos (Algoritmo Greedy).
    4. Entrena un modelo Random Forest Multi-Output para predecir la susceptibilidad.
    5. Identifica marcadores biológicos (motivos de ADN) asociados al éxito/fracaso viral.

REQUISITOS:
    - Python 3.8+
    - Librerías: pandas, numpy, matplotlib, seaborn, scikit-learn, biopython, tqdm
    - Instalación: pip install pandas numpy matplotlib seaborn scikit-learn biopython tqdm
===============================================================================
"""

import os
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO
from itertools import product
from tqdm import tqdm
from sklearn.model_selection import train_test_split
from sklearn.multioutput import MultiOutputClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score

# ==============================================================================
# CONFIGURACIÓN DEL ENTORNO
# ==============================================================================
# RUTAS DE ARCHIVOS (Ajustar según la máquina donde se ejecute)
# Ruta a la carpeta con los genomas FASTA (.fna)
PATH_GENOMES = "/mnt/d/genomes/e.coli/filtered_1000_genomes/*.fna" 
# Ruta a la base de datos de interacciones (Tabla Plana)
PATH_INTERACTIONS = "/mnt/d/genomes/e.coli/exports del DB/PLANE TABLE TO ml.csv"
# Ruta donde se guardará/cargará la matriz procesada de K-mers
PATH_KMERS_MATRIX = "X_Kmers_Matrix.csv"

# PARÁMETROS DEL MODELO
K_SIZE = 4            # Longitud del k-mer (Tetranucleótidos -> 256 características)
N_COCKTAIL = 20       # Número de fagos a seleccionar para el cocktail terapéutico
TEST_SIZE = 0.2       # Porcentaje de datos reservados para validación (20%)
RANDOM_SEED = 42      # Semilla para reproducibilidad

# ==============================================================================
# BLOQUE 1: EXTRACCIÓN DE CARACTERÍSTICAS GENÓMICAS (K-MERS)
# ==============================================================================
def generate_kmers(k):
    """Genera todas las combinaciones posibles de nucleótidos de longitud k."""
    bases = ['A', 'C', 'G', 'T']
    return [''.join(p) for p in product(bases, repeat=k)]

def get_kmer_counts(file_path, kmer_list):
    """Lee un genoma y calcula la frecuencia normalizada de cada k-mer."""
    counts = {k: 0 for k in kmer_list}
    total_kmers = 0
    try:
        with open(file_path, "r") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                seq = str(record.seq).upper()
                for i in range(len(seq) - K_SIZE + 1):
                    kmer = seq[i : i + K_SIZE]
                    if kmer in counts:
                        counts[kmer] += 1
                        total_kmers += 1
    except Exception as e:
        print(f"Error leyendo {file_path}: {e}")
        return None, 0
    
    # Normalización (Frecuencia Relativa)
    if total_kmers > 0:
        for k in counts:
            counts[k] /= total_kmers
    return counts, total_kmers

def build_genomic_matrix():
    """Construye o carga la matriz X de características genómicas."""
    if os.path.exists(PATH_KMERS_MATRIX):
        print(f"✅ Cargando matriz de K-mers existente desde: {PATH_KMERS_MATRIX}")
        return pd.read_csv(PATH_KMERS_MATRIX, index_col='assembly_accession')
    
    print("⚠️ No se encontró matriz pre-calculada. Iniciando extracción desde genomas crudos...")
    files = glob.glob(PATH_GENOMES)
    if not files:
        raise FileNotFoundError(f"No se encontraron archivos .fna en {PATH_GENOMES}")
        
    all_kmers = generate_kmers(K_SIZE)
    data_rows = []
    ids = []
    
    print(f"Procesando {len(files)} genomas...")
    for f in tqdm(files):
        filename = os.path.basename(f)
        accession = filename.replace('.fna', '') # Limpieza del ID
        kmer_profile, total = get_kmer_counts(f, all_kmers)
        if kmer_profile and total > 0:
            data_rows.append(kmer_profile)
            ids.append(accession)
            
    df = pd.DataFrame(data_rows)
    df['assembly_accession'] = ids
    df.set_index('assembly_accession', inplace=True)
    df.to_csv(PATH_KMERS_MATRIX)
    print(f"✅ Matriz generada y guardada: {df.shape}")
    return df

# ==============================================================================
# BLOQUE 2: PROCESAMIENTO Y ALINEACIÓN DE DATOS
# ==============================================================================
def load_and_align_data():
    print("\n--- 1. CARGA Y ALINEACIÓN DE DATOS ---")
    
    # 1. Cargar Features (X)
    X = build_genomic_matrix()
    
    # 2. Cargar Targets (Y) - Resultados de Laboratorio
    print(f"Cargando base de datos de interacciones desde: {PATH_INTERACTIONS}")
    df_db = pd.read_csv(PATH_INTERACTIONS)
    df_db["y"] = df_db["resultado"].map({"Susceptible": 1, "Resistente": 0})
    
    # Pivotar para crear matriz de interacciones (Filas=Bacterias, Col=Fagos)
    # fill_value=0 asume resistencia si no hay dato experimental
    Y_raw = df_db.pivot_table(index='assembly_accession', 
                              columns='virus_accesion', 
                              values='y', 
                              fill_value=0)
    
    # 3. Intersección (Alineación X-Y)
    common_ids = X.index.intersection(Y_raw.index)
    print(f"Bacterias con Genoma (X): {len(X)} | Bacterias con Fenotipo (Y): {len(Y_raw)}")
    print(f" >> Bacterias Coincidentes para Entrenamiento: {len(common_ids)}")
    
    if len(common_ids) == 0:
        raise ValueError("Error crítico: No hay coincidencia entre IDs de genomas y base de datos.")
        
    return X.loc[common_ids], Y_raw.loc[common_ids]

# ==============================================================================
# BLOQUE 3: SELECCIÓN DEL COCKTAIL (ALGORITMO GREEDY + RELLENO)
# ==============================================================================
def select_optimal_cocktail(Y_matrix, n_phages=20):
    print(f"\n--- 2. SELECCIÓN DE COCKTAIL ({n_phages} FAGOS) ---")
    selected = []
    current_matrix = Y_matrix.copy()
    
    for i in range(n_phages):
        if current_matrix.shape[0] == 0:
            print("   INFO: Cobertura total alcanzada. Rellenando con fagos de alta potencia...")
            remaining = [p for p in Y_matrix.columns if p not in selected]
            ranking = Y_matrix[remaining].sum().sort_values(ascending=False)
            extras = ranking.head(n_phages - len(selected)).index.tolist()
            selected.extend(extras)
            break
            
        # Selección Greedy: El fago que mata más bacterias remanentes
        scores = current_matrix.sum()
        best_phage = scores.idxmax()
        score = scores.max()
        
        selected.append(best_phage)
        
        # Eliminar bacterias ya cubiertas
        dead_indices = current_matrix[current_matrix[best_phage] == 1].index
        current_matrix = current_matrix.drop(dead_indices)
        
        print(f"   Fago #{i+1}: {best_phage} (Cobertura incremental: {score} cepas)")
        
    return selected

# ==============================================================================
# BLOQUE 4: ENTRENAMIENTO Y EVALUACIÓN (MACHINE LEARNING)
# ==============================================================================
def train_and_evaluate(X, Y, phage_names):
    print("\n--- 3. ENTRENAMIENTO DEL MODELO (RANDOM FOREST) ---")
    
    # División Train/Test
    X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=TEST_SIZE, random_state=RANDOM_SEED)
    
    # Modelo Multi-Output (Predice susceptibilidad para los 20 fagos simultáneamente)
    model = MultiOutputClassifier(RandomForestClassifier(n_estimators=100, random_state=RANDOM_SEED, n_jobs=-1))
    model.fit(X_train, Y_train)
    
    # Predicción
    Y_pred = model.predict(X_test)
    
    # Evaluación por fago
    accuracies = [accuracy_score(Y_test.iloc[:, i], Y_pred[:, i]) for i in range(len(phage_names))]
    global_acc = np.mean(accuracies)
    
    print(f"\nResultados de Validación:")
    print(f" >> Precisión Global Promedio: {global_acc:.2%}")
    
    # Gráfico de Precisión
    plt.figure(figsize=(12, 6))
    eje_x = list(range(1, len(phage_names) + 1))
    sns.barplot(x=eje_x, y=accuracies, palette="viridis")
    plt.title(f'Desempeño Predictivo por Fago del Cocktail (Promedio: {global_acc:.2%})')
    plt.xlabel('Fago # (Ranking)')
    plt.ylabel('Exactitud (Accuracy)')
    plt.ylim(0, 1.05)
    plt.axhline(global_acc, color='red', linestyle='--', label='Promedio Global')
    plt.legend()
    plt.show()
    
    return model, global_acc

# ==============================================================================
# BLOQUE 5: INTERPRETABILIDAD BIOLÓGICA (MOTIVOS GENÉTICOS)
# ==============================================================================
def analyze_feature_importance(model, X_features, Y_matrix, phage_names):
    print("\n--- 4. ANÁLISIS DE BIOMARCADORES (K-MERS) ---")
    
    # Buscar un fago "interesante" (ni 100% letal ni 0% letal) para ver patrones
    kill_rates = Y_matrix.mean()
    balanced_phages = kill_rates[(kill_rates > 0.2) & (kill_rates < 0.8)].index.tolist()
    
    target_phage = balanced_phages[0] if balanced_phages else phage_names[1]
    print(f"Analizando determinantes genéticos para el fago: {target_phage}")
    
    # Extraer importancias del modelo específico para ese fago
    idx = phage_names.index(target_phage)
    importances = model.estimators_[idx].feature_importances_
    
    # Ranking Top 15 K-mers
    indices = np.argsort(importances)[::-1][:15]
    top_kmers = X_features.columns[indices]
    top_scores = importances[indices]
    
    # Gráfico
    plt.figure(figsize=(12, 6))
    plt.title(f"Firmas Genómicas Determinantes para Susceptibilidad a {target_phage}")
    plt.bar(range(15), top_scores, color="#d62728")
    plt.xticks(range(15), top_kmers, rotation=45)
    plt.ylabel("Importancia Relativa (Gini)")
    plt.xlabel("Motivo Tetranucleotídico (4-mer)")
    plt.tight_layout()
    plt.show()

# ==============================================================================
# EJECUCIÓN PRINCIPAL
# ==============================================================================
if __name__ == "__main__":
    try:
        # 1. Preparar Datos
        X, Y_all = load_and_align_data()
        
        # 2. Seleccionar Cocktail
        top_phages = select_optimal_cocktail(Y_all, n_phages=N_COCKTAIL)
        Y_cocktail = Y_all[top_phages]
        
        # 3. Entrenar y Evaluar
        model, acc = train_and_evaluate(X, Y_cocktail, top_phages)
        
        # 4. Análisis Biológico
        analyze_feature_importance(model, X, Y_cocktail, top_phages)
        
        print("\n✅ EJECUCIÓN COMPLETADA EXITOSAMENTE.")
        
    except Exception as e:
        print(f"\n❌ ERROR CRÍTICO EN LA EJECUCIÓN:\n{e}")