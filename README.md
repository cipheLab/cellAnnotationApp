# Cell Annotation Tool

## Description
The Cell Annotation Tool is a Shiny web application designed to facilitate the annotation of cytometry data. This tool allows users to upload FCS files, apply necessary preprocessing transformations, and use various algorithms for cell annotation. The results can be visualized and downloaded for further analysis.

## Features
- **Upload and View Data**: Load one or more FCS files and view their contents.
- **Preprocessing**: Apply compensation and transformation to the data.
- **Annotation Algorithms**: 
  - **XGBoost**: Use a machine learning model for cell annotation.
  - **Scyan**: Apply the Scyan algorithm with a knowledge table for annotation.
  - **Scaffold**: Perform clustering and build Scaffold maps for annotation.
- **Results Visualization**: View and analyze the annotation results.
- **Download**: Export the annotated FCS files and results.

## Installation
To run the application, you need to have R and the required packages installed. 

1. Clone the repository:
   ```sh
   git clone https://github.com/Maelle83/CellAnnotation.git
   cd CellAnnotation
   ```
2. Install the required R packages:
 ```r
install.packages(c("shiny", "shinydashboard", "shinyjs", "shinybusy", "bslib", "xgboost", "dplyr", "tidyr", "flowCore", "DT", "FlowCIPHE", "plyr", "reticulate", "gtools", "igraph", "Rcpp", "openxlsx", "parallel"))
```
3. Run the Shiny application:

 ```r
shiny::runApp()
```

   ## Usage

### Upload and View Data
1. Navigate to the **Upload FCS** tab.
2. Click on the "Upload..." button to select and upload one or multiple FCS files.
3. View the contents of the uploaded files in the table.

### Preprocessing
1. Go to the **Preprocessing** tab.
2. Select the markers to transform and the type of transformation (logicle or arcsinh).
3. Click "Apply" to perform the preprocessing.

### Annotation

#### XGBoost
1. Navigate to the **XGBoost** tab.
2. Load or choose an XGBoost model.
3. Select markers used for prediction.
4. Click "Annotate Selected Files" to perform the annotation.

#### Scyan
1. Go to the **Scyan** tab.
2. Load or build a Scyan knowledge table.
3. Select markers for prediction.
4. Click "Run Scyan Annotation" to perform the annotation.

#### Scaffold
1. Navigate to the **Scaffold** tab.
2. Perform clustering or upload already clustered files.
3. Build or load a Scaffold map.
4. Click "Run Scaffold Annotation" to perform the annotation.

### Results Visualization and Download
1. Go to the **Results** tab to view annotation statistics and visualizations.
2. Click "Download ALL CSV results" to export the results.

## System Requirements
- **OS**: Linux, Windows, or macOS
- **R version**: 3.6 or higher
- **Python version**: 3.8 or higher

<details>
<summary>Dependencies versions</summary>
<br>
Paste here what 'pip list' gives you.
</details>

## Contributing
Contributions are welcome! Please feel free to submit a Pull Request or open an issue.

## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contact
If you have any questions, feel free to reach out to [maelle-marine.monier@inserm.fr].

