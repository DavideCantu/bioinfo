# bioinfo
Progetto di Bioinformatica - Anno 2024/2025 - Tema 4

Dato un file FASTQ di reads della stessa lunghezza, una soglia di frequenza F (0 ≤ F ≤ 1) e una lunghezza k, produrre in output un report testuale dell'uso dei k-mers (sottostringhe di lunghezza k) per posizione nei reads in input. Precisamente, riportare, per ciascun k-mer con frequenza almeno F nel dataset di input, quante volte occorre in ciascuna posizione dei reads.

Selezionare poi il k-mer K che appare più volte in una posizione e visualizzare un diagramma dell’andamento del numero di occorrenze di K rispetto alle posizioni dei reads.

Produrre poi un file FASTA con i reads che contengono K esattamente nella posizione dove K occorre più volte, riportando nell’header FASTA la qualità media delle basi del read, oltre al suo identificatore.

https://github.com/bioinformatica-corso/lezioni/tree/master/laboratorio/progetto4
