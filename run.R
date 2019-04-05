#!/usr/local/bin/Rscript

task <- dyncli::main()

library(dplyr, warn.conflicts = FALSE)
library(purrr, warn.conflicts = FALSE)

library(monocle, warn.conflicts = FALSE)
library(dynwrap, warn.conflicts = FALSE)

#   ____________________________________________________________________________
#   Load data                                                               ####

parameters <- task$parameters
counts <- as.matrix(task$counts)

#   ____________________________________________________________________________
#   Infer trajectory                                                        ####

# just in case
if (is.factor(parameters$norm_method)) {
  parameters$norm_method <- as.character(parameters$norm_method)
}

# TIMING: done with preproc
checkpoints <- list(method_afterpreproc = as.numeric(Sys.time()))

# load in the new dataset
pd <- Biobase::AnnotatedDataFrame(data.frame(row.names = rownames(counts)))
fd <- Biobase::AnnotatedDataFrame(data.frame(row.names = colnames(counts), gene_short_name = colnames(counts)))
cds <- monocle::newCellDataSet(t(counts), pd, fd)

# estimate size factors and dispersions
cds <- BiocGenerics::estimateSizeFactors(cds)
cds <- BiocGenerics::estimateDispersions(cds)

# filter features if requested
if (parameters$filter_features) {
  disp_table <- dispersionTable(cds)
  ordering_genes <- subset(disp_table, mean_expression >= parameters$filter_features_mean_expression)
  cds <- setOrderingFilter(cds, ordering_genes)

  print(nrow(ordering_genes))
}

# if low # cells or features -> https://github.com/cole-trapnell-lab/monocle-release/issues/26
# this avoids the error "initial centers are not distinct."
if (ncol(counts) < 500 || nrow(counts) < 500) {
  parameters$auto_param_selection <- FALSE
}

# reduce the dimensionality
cds <- monocle::reduceDimension(
  cds,
  max_components = parameters$max_components,
  reduction_method = parameters$reduction_method,
  norm_method = parameters$norm_method,
  auto_param_selection = parameters$auto_param_selection
)

# order the cells
cds <- monocle::orderCells(cds)

# TIMING: done with method
checkpoints$method_aftermethod <- as.numeric(Sys.time())

# construct dimred output
dimred <- t(reducedDimS(cds))
colnames(dimred) <- paste0("Comp", seq_len(ncol(dimred)))

# construct milestone network
degs <- igraph::degree(cds@minSpanningTree)
milestone_ids <- names(which(degs != 2))

mst <- dynwrap::simplify_igraph_network(cds@minSpanningTree)
milestone_network <- mst %>% igraph::as_data_frame() %>% rename(length = weight)

# construct milestone dimred output
dimred_milestones <- t(cds@reducedDimK)[milestone_ids, , drop = FALSE]
colnames(dimred_milestones) <- paste0("Comp", seq_len(ncol(dimred_milestones)))

# rename milestones
milestone_id_map <- set_names(paste0("M", seq_along(milestone_ids)), milestone_ids)
rownames(dimred_milestones) <- milestone_id_map[rownames(dimred_milestones)]
milestone_ids <- milestone_id_map[milestone_ids] %>% set_names(NULL)
milestone_network <- milestone_network %>% mutate(from = milestone_id_map[from] %>% set_names(NULL), to = milestone_id_map[to] %>% set_names(NULL))

# wrap output
cell_ids <- rownames(dimred)
output <- lst(
  cell_ids,
  milestone_ids,
  milestone_network,
  dimred,
  dimred_milestones,
  timings = checkpoints
)

#   ____________________________________________________________________________
#   Save output                                                             ####

output <- dynwrap::wrap_data(cell_ids = cell_ids) %>%
  dynwrap::add_dimred_projection(
    milestone_network = milestone_network,
    dimred = dimred,
    dimred_milestones = dimred_milestones
  ) %>%
  dynwrap::add_timings(checkpoints)

dyncli::write_output(output, task$output)
