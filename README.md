Data Download and Preprocessing
To begin, download bulk RNA-Seq data from the following publicly available studies:

| Dataset ID | Cell Type          | Variant                | Sample Info                                                                  | Location & Time Period             |
| ---------- | ------------------ | ---------------------- | ---------------------------------------------------------------------------- | ---------------------------------- |
| GSE157103  | Leukocytes         | Original Wuhan         | 10 COVID-negative; 48 COVID-positive (hospitalized/ICU); 50 non-hospitalized | New York, USA (2020)               |
| GSE171110  | Whole blood        | French COVID           | 10 COVID-negative; 44 COVID-positive (severe)                                | Créteil, France (2020)             |
| GSE189039  | PBMCs              | Beta                   | 8 COVID-negative; 9 COVID-positive                                           | South Tyrol, Italy (2021)          |
| GSE201530  | PBMCs              | Omicron                | 8 COVID-negative; 38 COVID-positive                                          | Tyrol, Austria (2021–2022)         |
| GSE152418  | PBMCs              | Original Wuhan         | 17 COVID-negative; 16 COVID-positive                                         | Atlanta, Georgia, USA (2020)       |
| PMC8202013 | Whole blood        | Original Wuhan         | 27 COVID-negative; 103 COVID-positive                                        | Lausanne, Switzerland (2020)       |
| GSE161731  | Whole blood        | Original Wuhan         | 16 COVID-negative; 12 COVID-positive                                         | Durham, North Carolina, USA (2020) |
| GSE166190  | Whole blood        | Original Wuhan         | 11 COVID-negative; 10 COVID-positive                                         | Geneva, Switzerland (2020)         |
| GSE294888  | pDCs and DC2s      | Delta and Omicron BA.1 | 30 total samples (5 replicates each per condition)                           | Paris, France (2025)               |
| GSE239595  | NP lymphoid tissue | Omicron                | 9 total (3 COVID-negative; 6 COVID-positive)                                 | Seoul, South Korea (2022–2023)     |


File Format
Each dataset should be formatted as follows:

Samples are columns, genes are rows.

The first row must include column names.

The first column must be titled Gene Symbol, listing the gene names.

Sample names should be cleaned to use consistent naming across datasets, using (1), (2), etc., to distinguish replicates within each group.

Pipeline Steps
Once all datasets are correctly formatted, run the following scripts in order:

0000000000-part 1-DE analysis.R – to generate differentially expressed genes (DGEs).

0000000000-part 2-classification.ipynb – for classification analysis.

0000000000-part 3-GO and Pathway analysis.ipynb – for pathway and GO enrichment analysis.

Note: The script Pathway and GO.R must be located in the same directory as the notebooks for part 3 to run properly.
