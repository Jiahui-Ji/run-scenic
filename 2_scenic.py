#load required modules
import os
import glob
import pickle
import pandas as pd
import numpy as np

from dask.diagnostics import ProgressBar

from arboreto.utils import load_tf_names
from arboreto.algo import grnboost2

from pyscenic.rnkdb import FeatherRankingDatabase as RankingDatabase
from pyscenic.utils import modules_from_adjacencies, load_motifs
from pyscenic.prune import prune2df, df2regulons
from pyscenic.aucell import aucell

import seaborn as sns





#parameter
celltype="Tcells0"

DATA_FOLDER="/data/jj1419/t-cell/tmp"
RESOURCES_FOLDER="/data/jj1419/t-cell/resources"
DATABASE_FOLDER = "/data/jj1419/t-cell/databases"
SCHEDULER="123.122.8.24:8786"
DATABASES_GLOB = os.path.join(DATABASE_FOLDER, "mm10*.mc9nr.feather")

MOTIF_ANNOTATIONS_FNAME = os.path.join(RESOURCES_FOLDER, "motifs-v9-nr.mgi-m0.001-o0.0.tbl")
HS_TFS_FNAME = os.path.join(RESOURCES_FOLDER, 'mm_mgi_tfs.txt')

sc_input=celltype+".csv"
regulon_fname=celltype+"_regu.p"
motif_fname=celltype+"_motifs.csv"
output_csv=celltype+"_AUCMatrix.csv"
output_loom=celltype+".loom"

SC_EXP_FNAME = os.path.join(RESOURCES_FOLDER, sc_input)
REGULONS_FNAME = os.path.join(DATA_FOLDER, regulon_fname)
MOTIFS_FNAME = os.path.join(DATA_FOLDER, motif_fname)



##load data
#load expression data
ex_matrix = pd.read_csv(SC_EXP_FNAME, sep=',', header=0, index_col=0).T
ex_matrix.shape
#load transcription factor
tf_names = load_tf_names(HS_TFS_FNAME)
#load database
db_fnames = glob.glob(DATABASES_GLOB)
def name(fname):
    return os.path.splitext(os.path.basename(fname))[0]
    
     
dbs = [RankingDatabase(fname=fname, name=name(fname)) for fname in db_fnames]
dbs

##co-expression modules
adjacencies = grnboost2(ex_matrix, tf_names=tf_names, verbose=True)
modules = list(modules_from_adjacencies(adjacencies, ex_matrix))

##prune co-expression modules using tf-motif enrichment analysis 
# Calculate a list of enriched motifs and the corresponding target genes for all modules.
with ProgressBar():
    df = prune2df(dbs, modules, MOTIF_ANNOTATIONS_FNAME,
                    client_or_address="custom_multiprocessing", num_workers=8)

# Create regulons from this table of enriched motifs.
regulons = df2regulons(df)

# Save the enriched motifs and the discovered regulons to disk.
df.to_csv(MOTIFS_FNAME)
with open(REGULONS_FNAME, "wb") as f:
    pickle.dump(regulons, f)

## cellular regulon enrichment matrix (AUCell scoring the regulons)
auc_mtx = aucell(ex_matrix, regulons, num_workers=4)
auc_mtx.to_csv(output_csv)

#image
img=sns.clustermap(auc_mtx, figsize=(8,8))

##save loom file
from pyscenic.export import export2loom

export2loom(ex_mtx = ex_matrix, auc_mtx = auc_mtx, regulons = [r.rename(r.name.replace('(',' (')) for r in regulons], out_fname = output_loom)





