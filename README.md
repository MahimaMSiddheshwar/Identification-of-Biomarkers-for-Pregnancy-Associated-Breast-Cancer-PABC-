# Identification of Biomarkers for Pregnancy Associated Breast Cancer (PABC)                                                                           


**Name:** Mahima Mahabaleshwar Siddheshwar  
**Programming Language:** R  
**Operating System:** Windows  
**Version:** R 4.3.1  
**File Format:** R Script  

### **Description**
This project focuses on identifying biomarkers for **Pregnancy Associated Breast Cancer (PABC)** using gene expression data from the **GSE31192** dataset. The analysis pipeline includes:

- **Importing** gene expression data.
- **Preprocessing** and normalization.
- **Statistical analysis** to identify differentially expressed genes (DEGs).
- **Pathway enrichment analysis** using **Gene Ontology (GO)** and **KEGG pathways** to assess the biological significance of DEGs.

The primary objective is to identify and interpret key genes associated with PABC, providing insights into the underlying biological mechanisms.

---

## **Input Files**
- **GSE31192 dataset:** Downloaded from the [Gene Expression Omnibus (GEO)](https://www.ncbi.nlm.nih.gov/geo/).

---

## **Required Files and Dependencies**
### **R Script Requirements**
- The R code is provided in a script file format. Ensure all paths are updated in the code to match your local environment.

### **Libraries/Packages**
Install the following R libraries before execution:
- **Data Acquisition:** GEOquery
- **Data Manipulation:** dplyr
- **Data Visualization:** ggplot2, pheatmap, ggrepel
- **Statistical Analysis:** limma
- **Biological Analysis:** clusterProfiler, org.Hs.eg.db

### **Software Requirements**
- **RStudio** for running and managing the R script.
- **Internet connectivity** for downloading data via the GEOquery package.

---

## **Execution Instructions**
1. **Download** the **GSE31192** dataset using the `GEOquery` package.
2. **Install** all required R libraries (see above).
3. **Open** the R script file in RStudio.
4. **Update** the hardcoded paths for data files and output directories in the script to match your local environment.
5. **Execute** each section of the script sequentially to generate results and visualizations.

---

## **Outputs**
### **Visualizations**
The project generates the following visual outputs in PNG format:
- **Box Plots** of gene expression data.
- **Heatmaps** to visualize clustering and expression patterns.
- **Volcano Plots** highlighting DEGs.
- **Scatter Plots** for comparative expression analysis.
- **MA Plots** for log fold-change versus mean expression.

### **Statistical Outputs**
- A comprehensive list of **Differentially Expressed Genes (DEGs)** with corresponding p-values and adjusted p-values.

#### **Example Files Generated:**
- `Box_plot.png`
- `Expression_Heatmap.png`
- Other visual outputs.

---

## **Notes**
- Ensure all dependencies are installed before running the scripts.
- **Internet connectivity** is required to fetch input data using GEOquery.
- Modify script thresholds and parameters for extended analysis.
- Adjust hardcoded paths to match your local file structure.

---

## **Additional Resources**
- A detailed code repository with version control is available on GitHub. Refer to the accompanying Word document for the repository link.
- For further analysis and modifications, explore documentation for the included libraries and packages.

---

## **Contact**
For any issues or queries regarding this project, please contact:  
**Mahima Mahabaleshwar Siddheshwar**


