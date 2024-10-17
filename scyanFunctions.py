import scyan 
import pandas as pd
import numpy as np



def readKnowledgeTable(pathKnowledgeTable):
    # Find file extension .csv, .xlsx, .xls)
    file_extension = pathKnowledgeTable.split('.')[-1]

    if file_extension == 'csv':
        
        table = pd.read_csv(pathKnowledgeTable, index_col=[0])
    elif file_extension in ['xlsx', 'xls']:
     
        table = pd.read_excel(pathKnowledgeTable, index_col=[0])

    # Define Population as index 
    if "Populations" in table.columns:
        table.set_index("Populations", inplace=True)
    
    return table

# Function to build the AnnDataObject that contain the expression matrix of the fcs file

def builAnnDataObject (data):
  
  adata=scyan.preprocess.AnnData(X=data.loc[:, data.columns != "Time"].values.astype(np.float32), var=pd.DataFrame(index=data.columns[data.columns != "Time"]), obs = data.loc[:, data.columns == "Time"])
  
  print(f"Created anndata object with {adata.n_obs} cells and {adata.n_vars} markers.\n\n-> The markers names are: {', '.join(adata.var_names)}\n-> The non-marker names are: {', '.join(adata.obs.columns)}")
 
  return adata


# Run the Scyan model

def runModel(adata, table) :
 
 model = scyan.Scyan(adata, table, prior_std=0.25, lr=0.0001)
 model.fit(accelerator="cpu", profiler="simple")

 return model


# Predict annotations

def predictAnnotation(modell):
  
  populations = model.predict()
  
  return populations

  
def test(self): 
    log_prob_th = float(-50)
    key_added: Optional[str] = "scyan_pop"
    add_levels: bool = True
   
   
   
# Predict probability for each population
    df = self.predict_proba()
    
    print(df)
     
    self.adata.obs["scyan_log_probs"] = df["max_log_prob_u"].values
    
    populations = df.iloc[:, :self.n_pops].idxmax(axis=1).astype("category")

    populations[df["max_log_prob_u"] < log_prob_th] = np.nan


    if key_added is not None:

        self.adata.obs[key_added] = pd.Categorical(
            populations, categories=self.pop_names
        )

    if add_levels and isinstance(self.table.index, pd.MultiIndex):
        utils._add_level_predictions(self, key_added)

    missing_pops = self.n_pops - len(populations.cat.categories)

    # if missing_pops:
    #     log.warning(
    #         f"{missing_pops} population(s) were not predicted. It may be due to:\n"
    #         f"  - Errors in the knowledge table (see https://mics-lab.github.io/scyan/advice/#advice-for-the-creation-of-the-table)\n"
    #         f"  - The model hyperparameters choice (see https://mics-lab.github.io/scyan/advanced/parameters/)\n"
    #         f"  - Or maybe these populations are really absent from this dataset."
    #     )
#     
#     return populations

