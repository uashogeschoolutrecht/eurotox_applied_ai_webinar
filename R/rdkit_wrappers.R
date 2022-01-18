library(reticulate)
library(tidyverse)
library(magrittr)
library(rsvg)
reticulate::conda_list()
reticulate::py_config()
reticulate::use_condaenv("r-reticulate",
                         conda = "/home/marc/.local/share/r-miniconda/envs/r-reticulate/bin/python")
Sys.getenv("LD_LIBRARY_PATH")
Sys.getenv("RDBASE")

rdkit   = reticulate::import("rdkit")
Chem    = rdkit$Chem
AllChem = rdkit$Chem$AllChem

openpng       = function(img) { grid::grid.raster(png::readPNG(img))}
opensvg       = function(svg) { tmp = tempfile(fileext = "png"); rsvg::rsvg_png(svg,tmp); openpng(tmp) }
showmols      = function(mols){ tmp = tempfile(fileext = "svg"); drawMols(mols,tmp); opensvg(tmp) }

defaultops    = Chem$Draw$DrawingOptions() %T>% (\(opts){opts$bgColor=py_none()})
drawMols      = function(mols,fname){ Chem$Draw$MolsToGridImage(mols,useSVG = T,molsPerRow = 3L) |> writeLines(con=fname) }
molFromSmiles  = Vectorize(function(smiles){Chem$MolFromSmiles(smiles)})
molFromSmiles2 = function(smiles){molFromSmiles(smiles) |> set_names(NULL)}
bvToIndex     = \(bv){ bv$GetOnBits() |> (\(ob){map_int(0:(bv$GetNumOnBits()-1),~ob[.x])})() }
bitvecToIndex = function(bitvec){lapply(bitvec,possibly(bvToIndex,otherwise = list()))}
rdkitFP       = Vectorize(function(mol){AllChem$GetMorganFingerprintAsBitVect(mol,3L,nBits=4096L)})
smilesToFP    = purrr::compose(bitvecToIndex,rdkitFP,molFromSmiles)
tanimoto    = function(a,b){length(intersect(a,b))/length(union(a,b))}
tanimoto.pb = function(fp,t){ pblapply(fp,function(a){tanimoto(a,t)}) |> unlist()}

#molFromSmiles2(c("c1ccccc1")) |> showmols()
## dt  = vroom::vroom("~/Downloads/Result_7.tsv",col_names = F) |> set_colnames(c("hazard","val","smiles"))
#dt2 = dt |> group_by(smiles) |> mutate(mol = smilesToFP(first(smiles))) |> ungroup()
