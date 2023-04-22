import os,sys
import loompy as lp
import numpy as np
import scanpy as sc
project=sys.argv[1]
output=sys.argv[2]
x=sc.read_csv(output+"/"+project+"_count_for_pyscenic.csv")
row_attrs={"Gene":np.array(x.var_names),}
col_attrs={"CellID":np.array(x.obs_names)}
lp.create(output+"/"+project+"_count_for_pyscenic.loom",x.X.transpose(),row_attrs,col_attrs)
