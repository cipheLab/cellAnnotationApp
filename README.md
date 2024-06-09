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
