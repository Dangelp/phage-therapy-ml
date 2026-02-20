# 🦠 Machine Learning for Precision Phage Therapy



!\[Python](https://img.shields.io/badge/Python-3.8%2B-blue)

!\[Scikit-Learn](https://img.shields.io/badge/Scikit--Learn-Machine%20Learning-orange)

!\[Biopython](https://img.shields.io/badge/Biopython-Genomics-green)

!\[Status](https://img.shields.io/badge/Status-Completed-success)



> \\\*\\\*Predicting \\\*Escherichia coli\\\* susceptibility to bacteriophages using alignment-free genomic signatures and Multi-Output Random Forests.\\\*\\\*



\## 📋 Overview

Antimicrobial Resistance (AMR) is a critical global health crisis. Bacteriophage therapy offers a highly specific alternative to traditional antibiotics. However, matching a clinical bacterial strain to the right phage currently relies on slow, in-vitro phenotypic screening (spot testing) that can take days.



This project introduces an \*\*in-silico, machine learning-driven approach\*\* to bypass manual screening. By extracting alignment-free genomic signatures (k-mers) directly from raw bacterial `.fna` assemblies, we can predict the therapeutic success of a 20-phage optimal cocktail with high accuracy, transforming a process that takes days into one that takes minutes.



\## 🔬 Methodology



Our bioinformatic pipeline consists of three core components:



1\. \*\*Alignment-Free Feature Extraction:\*\* Raw bacterial genomes (n=912) are processed to extract tetranucleotide frequencies (k=4). This generates a 256-dimensional "genomic fingerprint" capturing phylogenetic signatures and codon usage bias without requiring computationally expensive BLAST alignments or prior gene annotation.

2\. \*\*Greedy Set Cover Algorithm (Cocktail Selection):\*\* A greedy algorithmic approach selects the optimal combination of 20 bacteriophages (the "Golden Cocktail") that maximizes the lysis coverage across the \*E. coli\* cohort.

3\. \*\*Multi-Output Machine Learning:\*\*

   A `RandomForestClassifier` wrapped in a `MultiOutputClassifier` is trained to simultaneously predict the susceptibility of a clinical strain against all 20 phages in the cocktail.



\## 📊 Key Results



\* \*\*High Predictive Power:\*\* The system achieved a \*\*Global Accuracy of 82.32%\*\* on unseen test data, proving that raw genomic sequences contain sufficient biological signal to predict viral infection success.

\* \*\*Biological Interpretability:\*\* Feature importance analysis revealed that specific tetranucleotides (e.g., `TACG`) act as critical predictors for certain selective phages. This implicitly captures the presence of \*\*Restriction-Modification (R-M) systems\*\* and epigenetic signatures acting as bacterial defense mechanisms.



\## 🚀 Repository Structure



\\`\\\\`\\`text

phage-therapy-ml/

├── data/

│   ├── raw/                 # Sample raw .fna genomes

│   └── processed/           # X\_Kmers\_Matrix.csv \& Interactions DB

├── scripts/

│   └── TFM\_Modelo\_Predictivo\_Fagos.py  # Main pipeline script

├── results/

│   └── figures/             # Feature importance and accuracy plots

├── requirements.txt         # Dependencies

└── README.md

\\`\\\\`\\`



\## ⚙️ How to Run



1\. \*\*Clone the repository:\*\*

   \\`\\\\`\\`bash

   git clone https://github.com/YourUsername/phage-therapy-ml.git

   cd phage-therapy-ml

   \\`\\\\`\\`



2\. \*\*Install dependencies:\*\*

   \\`\\\\`\\`bash

   pip install -r requirements.txt

   \\`\\\\`\\`



3\. \*\*Run the predictive pipeline:\*\*

   \\`\\\\`\\`bash

   python scripts/TFM\_Modelo\_Predictivo\_Fagos.py

   \\`\\\\`\\`



\## 👨‍🔬 Author



\*\*David Angel Perez\*\* \*Microbiologist | M.Sc. Candidate in Bioinformatics\*



Bridging the gap between clinical microbiology and artificial intelligence to solve real-world public health challenges.



\* \[LinkedIn Profile](https://www.linkedin.com/in/tu-enlace-aqui/)

\* \[Contact Me](mailto:tu-correo@email.com)



---

\*This project was developed as the Master's Final Project (TFM) for the Master's Degree in Bioinformatics.\*

