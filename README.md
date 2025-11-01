# Clustering-Based Gene Variants Identification

## Project Description

This project applies clustering algorithms to DNA sequencing reads from deer to identify different variants of the Major Histocompatibility Complex (MHC) gene. The MHC plays a crucial role in the immune system by enabling the organism to recognize foreign molecules such as bacteria or viruses. Deer with a higher diversity of MHC variants can potentially respond to a broader range of pathogens.

The provided dataset contains sequencing data from 47 FASTQ files, corresponding to 41 individual deer (identified by numbers J1 to J41). Some files share the same numeric identifier but differ in suffixes (e.g., GK, S, L), indicating multiple technical or experimental replicates for the same deer. Each file contains between 200 and 3,000 extractions of the MHC gene.

The program aims to:

1. Eliminate corrupted or erroneous sequences.  
2. Correct measurement errors by averaging the sequences to identify true gene variants.  
3. Detect and distinguish multiple gene variants within a single deer, which may arise from inherited mutations preserved by natural selection.  

#### Output :
The program takes a FASTQ file as input (representing one deer) and outputs the **k main gene variants**, where the value of k is specified by the user.

#### Full Report :

To get full details about the results of our project please refer to `/FinalReport.pdf`

## Usage

### Prerequisites
- C++ compiler (g++)
- SPOA library https://github.com/rvaser/spoa (easy to use on Linux)

### Compilation

After downloading the repository, and going int hte folder, compile the code using one of the following commands:

If SPOA is already installed globally

```g++ -o code code.cpp -lspoa```

or

If SPOA is not yet installed globally, specify the paths to the headers and library
(make sure to replace ~/spoa/spoa/include and ~/spoa/build/lib with the actual paths on your system).

```g++ -o code code.cpp -I~/spoa/spoa/include -L~/spoa/build/lib -lspoa```

### Execution

Run the compiled executable with the input file as the first argument, here is an example:

```./code fastq/J29_B_CE_IonXpress_005.fastq```

---

This project was completed as part of the Bioinformatics course at FER University in Zagreb.