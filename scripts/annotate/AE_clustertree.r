# thank you to Luke Zappia: 
# https://lazappi.id.au/posts/2017-07-19-building-a-clustering-tree/

# Description -------------------------------------------------------------
# Generate clustertree of clusters 
base::source(here::here("islet_cartography_scrna/scripts/misc/set_up.R"))
set.seed(1000)

# Graphs
library(igraph)

# Plotting
library(ggraph)
library(viridis)

# Data manipulation
library(tidyverse)

clusterings <- vroom::vroom(here::here("islet_cartography_scrna/data/annotate/files/leiden_clusterings.csv")) |> 
  tibble::column_to_rownames("barcode") |> 
  dplyr::rename_all(~stringr::str_replace(., "leiden_", "res."))

getEdges <- function(clusterings) {

    # Loop over the different resolutions
    transitions <- lapply(1:(ncol(clusterings) - 1), function(i) {

        # Extract two neighbouring clusterings
        from.res <- sort(colnames(clusterings))[i]
        to.res <- sort(colnames(clusterings))[i + 1]

        # Get the cluster names
        from.clusters <- sort(unique(clusterings[, from.res]))
        to.clusters <- sort(unique(clusterings[, to.res]))

        # Get all possible combinations
        trans.df <- expand.grid(FromClust = from.clusters,
                                ToClust = to.clusters)

        # Loop over the possible transitions
        trans <- apply(trans.df, 1, function(x) {
            from.clust <- x[1]
            to.clust <- x[2]

            # Find the cells from those clusters
            is.from <- clusterings[, from.res] == from.clust
            is.to <- clusterings[, to.res] == to.clust

            # Count them up
            trans.count <- sum(is.from & is.to)

            # Get the sizes of the two clusters
            from.size <- sum(is.from)
            to.size <- sum(is.to)

            # Get the proportions of cells moving along this edge
            trans.prop.from <- trans.count / from.size
            trans.prop.to <- trans.count / to.size

            return(c(trans.count, trans.prop.from, trans.prop.to))
        })

        # Tidy up the results
        trans.df$FromRes <- as.numeric(gsub("res.", "", from.res))
        trans.df$ToRes <- as.numeric(gsub("res.", "", to.res))
        trans.df$TransCount <- trans[1, ]
        trans.df$TransPropFrom <- trans[2, ]
        trans.df$TransPropTo <- trans[3, ]

        return(trans.df)
    })

    # Bind the results from the different resolutions together
    transitions <- do.call("rbind", transitions)
    
    transitions$FromClust <- factor(transitions$FromClust)
    transitions$ToClust   <- factor(transitions$ToClust)
    
    return(transitions)
}

edges <- getEdges(clusterings)
head(edges)

##   FromClust ToClust FromRes ToRes TransCount TransPropFrom TransPropTo
## 1         0       0     0.0   0.3        135     0.3600000           1
## 2         0       1     0.0   0.3        100     0.2666667           1
## 3         0       2     0.0   0.3         60     0.1600000           1
## 4         0       3     0.0   0.3         50     0.1333333           1
## 5         0       4     0.0   0.3         30     0.0800000           1
## 6         0       0     0.3   0.6          0     0.0000000           0

# Some of these columns are pretty obvious but the last three could do with an explanation. 
# TransCount is the number of cells that move along this edge. 
# TransPropFrom is the proportion of the cells in the lower resolution cluster that have made this 
# transition and TransPropTo is the proportion of cells in the higher resolution cluster that came from this edge.

# Getting the information about the nodes of the tree is easier as these just represent the clusters. 
# This function summarises the cluster information and converts it to long format.

getNodes <- function(clusterings) {
    nodes <- clusterings %>%
        gather(key = Res, value = Cluster) %>%
        group_by(Res, Cluster) %>%
        summarise(Size = n()) %>%
        ungroup() %>%
        mutate(Res = stringr::str_replace(Res, "res.", "")) %>%
        mutate(Res = as.numeric(Res), Cluster = as.numeric(Cluster)) %>%
        mutate(Node = paste0("R", Res, "C", Cluster)) %>%
        select(Node, everything())
}

nodes <- getNodes(clusterings)
head(nodes)

## # A tibble: 6 x 4
##     Node   Res Cluster  Size
##    <chr> <dbl>   <dbl> <int>
## 1   R0C0   0.0       0   375
## 2 R0.3C0   0.3       0   135
## 3 R0.3C1   0.3       1   100
## 4 R0.3C2   0.3       2    60
## 5 R0.3C3   0.3       3    50
## 6 R0.3C4   0.3       4    30

# Each node needs a unique ID which I have made by combining the resolution and cluster number. 
# We also record the number of cells in each cluster.

# Now we can build the graph we will use as the starting point for our plot. 
# Some of the possible edges between clusters will have no cells travelling along them so we filter them out. 
# We also remove edges that correspond to a small proportion (< 2%) of cells in the higher resolution cluster.

graph <- edges %>%
    # Remove edges without any cell...
    filter(TransCount > 0) %>%
    # ...or making up only a small proportion of the new cluster
    filter(TransPropTo > 0.02) %>%
    # Rename the nodes
    mutate(FromNode = paste0("R", FromRes, "C", FromClust)) %>%
    mutate(ToNode = paste0("R", ToRes, "C", ToClust)) %>%
    # Reorder columns
    select(FromNode, ToNode, everything()) %>%
    # Build a graph using igraph
    graph_from_data_frame(vertices = nodes)

print(graph)

## IGRAPH b1b93c3 DN-- 23 23 --
## + attr: name (v/c), Res (v/n), Cluster (v/n), Size (v/n),
## | FromClust (e/c), ToClust (e/c), FromRes (e/n), ToRes (e/n),
## | TransCount (e/n), TransPropFrom (e/n), TransPropTo (e/n)
## + edges from b1b93c3 (vertex names):
##  [1] R0C0  ->R0.3C0 R0C0  ->R0.3C1 R0C0  ->R0.3C2 R0C0  ->R0.3C3
##  [5] R0C0  ->R0.3C4 R0.3C1->R0.6C0 R0.3C0->R0.6C1 R0.3C4->R0.6C1
##  [9] R0.3C0->R0.6C2 R0.3C2->R0.6C3 R0.3C3->R0.6C4 R0.6C0->R0.9C0
## [13] R0.6C2->R0.9C1 R0.6C3->R0.9C2 R0.6C1->R0.9C3 R0.6C4->R0.9C4
## [17] R0.6C1->R0.9C5 R0.9C0->R1.2C0 R0.9C1->R1.2C1 R0.9C2->R1.2C2
## [21] R0.9C3->R1.2C3 R0.9C4->R1.2C4 R0.9C5->R1.2C5

# Plot our graph using the `tree` layout
pdf(here::here("islet_cartography_scrna/data/annotate/plot/clustertree.pdf"),
    width = 4, height = 4)
ggraph(graph, layout = "tree") +
    # Plot the edges, colour is the number of cells and transparency is the
    # proportion contribution to the new cluster
    geom_edge_link(arrow = arrow(length = unit(1, 'mm')),
                    end_cap = circle(2, "mm"), edge_width = 0.5,
                    aes(colour = log(TransCount), alpha = TransPropTo)) +
    # Plot the nodes, size is the number of cells
    geom_node_point(aes(colour = factor(Res),
                        size = Size)) +
    geom_node_text(aes(label = Cluster), size = 1) +
    # Adjust the scales
    scale_size(range = c(1, 5)) +
    scale_edge_colour_gradientn(colours = viridis(100)) +
    # Add legend labels
    guides(size = guide_legend(title = "Cluster Size", title.position = "top"),
            colour = guide_legend(title = "Clustering Resolution",
                                    title.position = "top"),
            edge_colour = guide_edge_colorbar(title = "Cell Count (log)",
                                                title.position = "top"),
            edge_alpha = guide_legend(title = "Cluster Prop",
                                        title.position = "top", nrow = 2)) +
    # Remove the axes as they don't really mean anything
  my_theme_void() +
  theme(legend.position = "right")
dev.off()