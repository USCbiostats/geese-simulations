#!/bin/sh
#SBATCH --account=lc_pdt
#SBATCH --partition=thomas
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=12GB
#SBATCH --job-name=parameter_estimates
#SBATCH --mail-type=ALL

library(aphylo)
library(geese)
library(coda)

# Fitting partially annotated trees --------------------------------------------

partially_annotated <- readRDS("data/candidate_trees.rds")

# Parsing the data
treeids <- sort(unique(names(partially_annotated)))
data    <- vector("list", length(treeids))
names(data) <- treeids
for (tn in treeids) {

  # Preparing annotations
  tmp_trees <- partially_annotated[names(partially_annotated) == tn]
  tmp_ann   <- lapply(tmp_trees, function(a) rbind(a$tip.annotation, a$node.annotation))
  tmp_ann   <- do.call(cbind, tmp_ann)

  # Creating the aphylo version

  if (sum(names(partially_annotated) == tn) > 1)
    adata_tm <- do.call(c, tmp_trees)
  else
    adata_tm <- tmp_trees[[1L]]

  # Offspring -> parent, we need to add the root
  tmp_tree  <- rbind(tmp_trees[[1]]$tree$edge[, 2:1], c(Ntip(tmp_trees[[1]]) + 1, -1))

  # Sorting the annotations accordingly
  tmp_ann <- tmp_ann[tmp_tree[,1],,drop=FALSE]
  tmp_ann <- lapply(1:nrow(tmp_ann), function(i) tmp_ann[i ,])

  # Computing polytomies
  polyt <- table(table(tmp_tree[,2]))
  polyt <- max(as.integer(names(polyt)))

  data[[tn]] <- list(
    edges = tmp_tree,
    ann   = tmp_ann,
    dpl   = with(tmp_trees[[1L]], c(tip.type, node.type))[tmp_tree[,1]] == 0L,
    genes = with(tmp_trees[[1L]]$tree, c(tip.label, node.label))[tmp_tree[,1]],
    tree  = adata_tm,
    polytomies = polyt,
    nfuns      = length(tmp_ann[[1]]),
    max_size   = 2^(( polyt + 1 ) * length(tmp_ann[[1]]) )
  )

}

saveRDS(data, file="data/model_data.rds")

