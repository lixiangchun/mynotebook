# Network analysis with R package `igraph`

```R
library(cgat)
library(igraph)

inputData <- readRDS(system.file("data/haptrigData.rds",  package="cgat"))

# Construct a network from data frame
g <- graph_from_data_frame(inputData$humanBioGRIDnonredundant, directed=FALSE)

# Simplify g
g <- simplify(g)

# Extract node name: V(g)$name

# Calculate connected components of network
clu <- clusters(g)

# Calculate modularity score
modularity_score <- modularity(g, clu$membership)

# Get the maximum subgraph
g <- induced_subgraph(g, names(clu$membership)[clu$membership==1])

# Extract a small subgraph
x <- induced_subgraph(g, intersect(names(clu$membership)[clu$membership==1], inputData$keggGenes$symbol))

# Perform greedy optimisation of modularity
fc <- cluster_fast_greedy(x) # alias: fastgreedy.community

# additional methods for graph clustering: cluster_walktrap 
# for more info, type `?cluster_` at R terminal.

# Plot fc as dendrogram
plot_dendrogram(fc)

```

Additional tutorial with examples can be found in [http://igraph.wikidot.com/community-detection-in-r](http://igraph.wikidot.com/community-detection-in-r).