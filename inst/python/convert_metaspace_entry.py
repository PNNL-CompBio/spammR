import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Initialize Metaspace client
from metaspace import SMInstance


def download_ion_data(project_id, fdr_val = 0.2, database=('SwissLipids','2018-02-02')):
    
    sm = SMInstance()
    #project_id = "2024-02-15_20h37m13s"
    #fdr_val = 0.2

    # Fetch dataset details
    ds = sm.dataset(id=project_id)

    ##which database
    ##swiss lipids for lipids,
    ##core metabolome for metabolites
    result = ds.results(database=database)

    ##filter by 20% FDR here
    result = result[result.fdr <= fdr_val]

    #feature meta is here
    fmeta = pd.DataFrame(result)
    #let's reduce the features
    fmeta = fmeta[['ionFormula','ion','mz','msm','moleculeNames','moleculeIds']]

    ##rowData will be result with ion as rownames
    all_rows = []
    for i in range(len(result)):
        row = result.iloc[i]
        (sf, adduct) = row.name
        #get image
        images = ds.isotope_images(sf,adduct)

        #chris says i care about first image
        image = images[0]

        #get pixel information
        df = pd.DataFrame(image)
        x_coord = df.columns
        y_coord = df.index
        df['y_coord'] = y_coord
        long_tab = df.melt(value_name = 'intensity',value_vars=x_coord, id_vars='y_coord',var_name = 'x_coord')
        long_tab['sample_id'] = ""+long_tab['x_coord'].map(str)+"_"+long_tab['y_coord'].map(str)
        long_tab['ion'] = row.ion
        all_rows.append(long_tab)

    newtab = pd.concat(all_rows) ##full long table here
    # Save the DataFrame to a CSV file
    newtab.to_csv("metaspace_data.csv", index = False)
    fmeta.to_csv('metaspace_feature_data.tsv', index = False)
    return(newtab,fmeta)

    print("Data downloaded and saved to 'metaspace_data.csv'.")
    
