# InDepthMap helper script and functions


######################################

# function to display a data frame with the Cluster members.

CreateClusterDF <- function(Memberlist){
  GeneNames = names(Memberlist)
  
  Cluster = as.numeric(Memberlist)
  
  DF <- data.frame(ClusterID = paste("C", Cluster, sep = ""), GeneNames = GeneNames) %>%
    dplyr::group_by(ClusterID) %>%
    summarize(Members = paste(GeneNames, collapse = ", "))
  
  return(DF)
}


# NetworkMapper should take the desired matrix, the list of genes and other specifications as arguments
NetworkMapper <- function(GeneList, Matrix, DependencyScoreDF, CorrType, nCoDeps, CorrThreshold, ClusteringMethod, eps, NodeBorders){
  
  # intersect the search list with the available genes of the CoDepMatrix
  # redundant, since we have based the search criteria on available genes
  AvCol <- dplyr::intersect(GeneList, colnames(Matrix))
  
  
  # lets gather the columns of the CoDepMatrix of interest
  CoDepencyDF <- AvCol %>%
    
    # for each search gene..
    lapply(function(x){
      # select all co deps of the gene name, but omit NA of the triangular matrix
      
      # if triangular input is used:
      CoDeps <- c(Matrix[, x], Matrix[x,])
      
      # if triangular matrix input is used:
      # CoDepsCols <- na.omit(Matrix[, x])
      # 
      # CoDepsRows <- na.omit(Matrix[x,])
      # 
      # # combine and remove NA
      # 
      # CoDeps <- c(CoDepsCols, CoDepsRows)
      
      if(CorrType == "Positive"){
        # top n codeps
        CoDepOutput <- names(sort(CoDeps, decreasing = TRUE)[1:nCoDeps])
        
        # create DF
        df <- data.frame(target = x, CoDeps = CoDepOutput)
        
        return(df)
      }
      
      if(CorrType == "Negative"){
        # top n negative correlation codeps
        CoDepOutput <- names(sort(CoDeps)[1:nCoDeps])
        
        # create DF
        df <- data.frame(target = x, CoDeps = CoDepOutput)
        
        return(df)
      }
      
      if(CorrType == "Both"){
        # top n negative and positive correlation codeps
        AbsCoDeps = abs(CoDeps)
        
        # create a way to mark which values are negatively and positively correlated
        # not important for now.
        
        CoDepOutput <- names(sort(AbsCoDeps, decreasing = TRUE)[1:nCoDeps])
        
        df <- data.frame(target = x, CoDeps = CoDepOutput)
        
        return(df)
      }
    }) %>%
    rbindlist
  
  # filter for unique entries to generate a search list
  SearchList <- unique(c(CoDepencyDF$CoDeps, AvCol))
  
  # define corr matrix
  CorrMat <- Matrix[SearchList, SearchList]
  
  # order Na to the corner
  CorrMat <- CorrMat[order(rowSums(!is.na(CorrMat))), order(colSums(is.na(CorrMat)))]
  
  # convert to complete matrix, if triangular input
  tCorrMat <- t(CorrMat)
  # 
  CorrMat[upper.tri(CorrMat)] <- tCorrMat[upper.tri(tCorrMat)]
  
  # define adjacency matrix
  AdjMatrix <- CorrMat
  
  if(CorrType == "Positive"){
    # for positive corr, silence values below thresholg
    AdjMatrix[AdjMatrix < CorrThreshold] <- 0
  }
  
  if(CorrType == "Negative"){
    
    # convert threshold to negative
    NegCorrThreshold <- - CorrThreshold
    
    # silence values below threshold
    AdjMatrix[AdjMatrix > NegCorrThreshold] <- 0
    
    AdjMatrix <- abs(AdjMatrix)
    
  }
  
  if(CorrType == "Both"){
    
    # convert adj matrix to absolute values
    AdjMatrix <- abs(AdjMatrix)
    
    # silence values below threshold
    AdjMatrix[AdjMatrix < CorrThreshold] <- 0
    
  }
  
  # set NA diagonal values to 0
  diag(AdjMatrix) <- 0
  
  
  # create network
  network <- graph_from_adjacency_matrix(AdjMatrix,
                                         mode = "undirected",
                                         #mode = "directed",
                                         weighted = TRUE)
  
  ###### delete zero connection vertices #######
  
  # find entries without connections
  isolated = which(igraph::degree(network) == 0)
  
  # remove search targets  from list, so they are always present
  isolated <- isolated[!(names(isolated) %in% AvCol)]
  
  # remove 0 connection genes from network
  network_deleted_vertices <- delete.vertices(network, isolated)
  
  # add colours to search genes
  #V(network_deleted_vertices)$color <- ifelse(V(network_deleted_vertices)$name %in% AvCol, "red", "white")
  
  
  ################################ clustering #################################
  
  if(ClusteringMethod == TRUE) {
    
    # Do tSNE+DBSCAN clustering 
    set.seed(100)
    
    # tSNE from distance matrix
    tSNE_output <- Rtsne(as.dist(1 - CorrMat),
                         theta = 0.0,
                         perplexity = 3,
                         pca = FALSE,
                         is_distance = TRUE)
    
    # seperate out coordinates
    tSNEplotDf = data.frame(x = tSNE_output$Y[, 1], y = tSNE_output$Y[, 2], gene.name = colnames(CorrMat))
    
    # use t-sne coordinates for dbscan cluster recognition
    # here one could consider an iterative approach to cluster size, eps,
    # and check whether genes are consistently clustered together.
    dbscan_res = dbscan(x = data.frame(x = tSNEplotDf$x, y = tSNEplotDf$y)
                        , eps = eps)
    
    clusters <- dbscan_res$cluster
    
    clusters <- as.character(clusters)
    clusters[clusters == "0"] <- NA
    
    # assign clusters
    tSNEplotDf$clusters <- clusters
    
    groups <- as.numeric(tSNEplotDf$clusters)
    
    names(groups) <- tSNEplotDf$gene.name
    
    # and remove isolated vertices/nodes from grouping vector
    groups <- groups[!(names(groups) %in% names(isolated))]
    
    # generate grouped edge list for forcenetwork
    NWD3_DF <- igraph_to_networkD3(network_deleted_vertices, group = groups)
    
    
    #################### t-sne plot output ##################
    # plot clusters
    
    tSNEPlot <- ggplot(tSNEplotDf, aes(x = x, y = y, color = clusters)) + 
                        geom_point(alpha = 0.8) #+
                        #scale_fill_viridis(discrete = TRUE) #+ 
      
                        # cluster labels
                        # geom_text_repel(data = labelDF, 
                        #                 aes(label = label),
                        #                 size = 4.5,
                        #                 min.segment.length = 1,) + 
                        
                        # POI labels
                        # geom_text_repel(data = subset(tSNEplotDf, tSNEplotDf$gene.name %in% cols),
                        #                 aes(label = gene.name),
                        #                 min.segment.length = 0.7,
                        #                 size = 2.5)
                      
  } else {
    
    # do walk trap clustering
    groups = cluster_walktrap(network_deleted_vertices) %>% membership()
    
    # output edgelist
    NWD3_DF <- igraph_to_networkD3(network_deleted_vertices, group = groups)
    
    tSNEPlot <- NULL
    
  }
  
  # adjust nodesize based on input genes
  NWD3_DF$nodes$size <- 0
  NWD3_DF$nodes$size[NWD3_DF$nodes$name %in% AvCol] <- 50
  
  
  ForcedNet <- forceNetwork(Links = NWD3_DF$links,
                            Nodes = NWD3_DF$nodes,
                            Source = 'source',
                            Target = 'target',
                            Value = 'value',
                            NodeID = 'name',
                            Group = 'group',
                            Nodesize = "size",
                            zoom = TRUE,
                            fontSize = 15,
                            fontFamily = "serif",
                            linkDistance = 100,
                            charge = -50,
                            linkColour = "#666",
                            opacity = 0.9,
                            legend = TRUE,
                            opacityNoHover = 0.96
  )
  ########################### colour node borders by overall cancer cell dependency ##############
  if(NodeBorders == TRUE){
    
    # add a border to the nodes displaying the overall dependency
    NWD3_DF$nodes <- left_join(NWD3_DF$nodes, 
                               dplyr::select(DependencyScoreDF, gene_name, Percentage),
                               by = c("name" = "gene_name")) 
    
    # create colour gradient based on Percentage cancer cell dependency
    NWD3_DF$nodes$border <- lapply(NWD3_DF$nodes$Percentage, function(x){
      
                                viridis(n = 100, option = "magma")[x + 1]
                                
                              }) %>% unlist()
    
    # add the borders to the network
    ForcedNet$x$nodes$border <- NWD3_DF$nodes$border
    ForcedNet <- htmlwidgets::onRender(ForcedNet, 
                                       'function(el, x) { d3.selectAll("circle").style("stroke", d => d.border); }')
    
  }
  
  
  # adjust border color, to colour input genes by ge
  # NWD3_DF$nodes$borders <- FALSE
  # 
  # NWD3_DF$nodes$borders[NWD3_DF$nodes$name %in% AvCol] <- "#F00" # red
  
  ############################ freeze animation ###############################
  # customJS <- 'function() {
  #               simulation = this;
  #               simulation.stop();
  #               for (var i = 0; i < 300; ++i) simulation.tick();
  #               simulation.nodes().forEach( function(d,i) {
  #                 d.cx = d.x;
  #                 d.cy = d.y;
  #               });
  #               simulation.restart();
  #             }
  #             '
  # htmlwidgets::onRender(ForcedNet, customJS)
  
  
  output <- list(ForcedNet, groups, tSNEPlot)
  
  return(output)
  #return(network_deleted_vertices)
  #return(CoDepencyDF)
  #return(members)
}





# for the walktrap clustering i believe we need a nother output. a list of vectors containing the genes names of the clusters

# plotting function from gProfiler2
# https://github.com/egonw/r-gprofiler2
# adjusted wrap to grid display

gostplot <- function(gostres, capped = TRUE, interactive = TRUE, pal = c("GO:MF" = "#dc3912",
                                                                         "GO:BP" = "#ff9900",
                                                                         "GO:CC" = "#109618",
                                                                         "KEGG" = "#dd4477",
                                                                         "REAC" = "#3366cc",
                                                                         "WP" = "#0099c6",
                                                                         "TF" = "#5574a6",
                                                                         "MIRNA" = "#22aa99",
                                                                         "HPA" = "#6633cc",
                                                                         "CORUM" = "#66aa00",
                                                                         "HP" = "#990099")
){
  # gostres is the GOSt response list (contains results and metadata)
  # This function will plot only the sources that were asked from the query
  
  if( is.null(pal) ){
    pal <- c("GO:MF" = "#dc3912",
             "GO:BP" = "#ff9900",
             "GO:CC" = "#109618",
             "KEGG" = "#dd4477",
             "REAC" = "#3366cc",
             "WP" = "#0099c6",
             "TF" = "#5574a6",
             "MIRNA" = "#22aa99",
             "HPA" = "#6633cc",
             "CORUM" = "#66aa00",
             "HP" = "#990099")
  }
  
  if (!("result" %in% names(gostres))) stop("Name 'result' not found from the input")
  if (!("meta" %in% names(gostres))) stop("Name 'meta' not found from the input")
  
  source_order <- logpval <- term_id <- opacity <- NULL
  term_size <- term_name <- p_value <- term_size_scaled <- NULL
  
  df <- gostres$result
  # Order of data sources comes metadata
  meta <- gostres$meta
  
  # make sure all the essential column names are there
  essential_names <- c("source_order", "term_size", "term_name", "term_id", "source", "significant")
  
  if (!(all(essential_names %in% colnames(df)))) stop(paste("The following columns are missing from the result:", paste0(setdiff(essential_names, colnames(df)), collapse = ", ")))
  
  if (!any(grepl("p_value", colnames(df)))) stop("Column 'p_value(s)' is missing from the result")
  
  # nr_of_terms for every source
  widthscale <- unlist(lapply(meta$query_metadata$sources, function(x) meta$result_metadata[[x]][["number_of_terms"]]))
  names(widthscale) <- meta$query_metadata$sources # all the sources that were queried
  
  # Define the start positions for sources in the plot
  
  # start position for every term
  space <- 1000 # space between different sources
  starts <- c()
  start <- 1
  starts[1] <- start
  
  if(!length(widthscale) < 2) {
    for(idx in 2:length(widthscale)){
      starts[idx] <- starts[idx - 1] + space + widthscale[idx - 1]
    }
  }
  
  names(starts) <- names(widthscale)
  
  # Make sure that all the sources have colors
  
  if (is.null(names(pal))){
    names(pal) = meta$query_metadata$sources[1:length(pal)]
  }
  
  sourcediff = setdiff(meta$query_metadata$sources, names(pal))
  colors = grDevices::colors(distinct = TRUE)[grep('gr(a|e)y|white|snow|khaki|lightyellow', grDevices::colors(distinct = TRUE), invert = TRUE)]
  
  if (length(sourcediff) > 0){
    use_cols = sample(colors, length(sourcediff), replace = FALSE)
    pal[sourcediff] <- use_cols
  }
  
  # If multiquery
  if("p_values" %in% colnames(df)){
    p_values <- query <- significant <- NULL
    # spread the data frame to correct form
    df$query <- list(names(meta$query_metadata$queries))
    df <- tidyr::unnest(data = df, cols = c(p_values, query, significant))
    df <- dplyr::rename(df, p_value = p_values)
  }
  
  # Set sizescale of circles
  logScale <- function(input, input_start = 1, input_end = 50000, output_start = 2, output_end = 10){
    m = (output_end - output_start)/(log(input_end) - log(input_start))
    b = -m*log(input_start) + output_start
    output = m * log(input) + b
    return(output)
  }
  
  # Scale x positions
  xScale <- function(input, input_start = 1, input_end = sum(widthscale) + (length(widthscale) - 1)*space, output_start = 2, output_end = 200){
    m = (output_end - output_start)/(input_end - input_start)
    b = -m*input_start + output_start
    output = m * input + b
    return(output)
  }
  
  # Add values to df needed for plotting
  # add -log10 pval
  df$logpval <- -log10(df$p_value)
  df$opacity <- ifelse(df$significant, 0.8, ifelse(df$p_value == 1, 0, 0.3))
  df$term_size_scaled = logScale(df$term_size)
  # add x axis position
  df <- df %>% dplyr::group_by(source) %>% dplyr::mutate(order = xScale(source_order, input_start = 1, input_end = widthscale[source], output_start = starts[source], output_end = starts[source] + widthscale[source]))
  df$order <- xScale(df$order)
  
  if (capped) {
    df$logpval[df$logpval > 16] <- 17
    ymin <- -1
    ymax <- 18.5
    ticklabels <- c("0", "2", "4", "6", "8", "10", "12", "14", ">16")
    tickvals <- c(0, 2, 4, 6, 8, 10, 12, 14, 16)
  } else {
    ymin <- -1
    ymax <- ceiling(max(df$logpval)) + 5
    ticklabels <- ggplot2::waiver()
    tickvals <- ggplot2::waiver()
  }
  
  if (interactive){
    # Start plotting
    sd <- crosstalk::SharedData$new(df, key = ~term_id)
  } else {
    sd <- df
  }
  
  p <- ggplot2::ggplot(data = sd, ggplot2::aes(x = order, y = logpval, text = paste(term_id, paste0('(', term_size,')'), '<br>', term_name, '<br>', formatC(p_value, format = "e", digits = 3)))) +
    ggplot2::geom_point(ggplot2::aes(color = source, size = term_size_scaled, alpha = opacity),
                        show.legend = FALSE) +
    ggplot2::facet_grid(query ~ ., 
                        #ncol = 1,
                        #rows = length(unique(sd$query)),
                        scales = "free_x",
                        space = "free_x",
                        shrink = FALSE) +
    ggplot2::ylab("-log10(p-adj)") +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position='none',
                   panel.border=ggplot2::element_blank(),
                   strip.text=ggplot2::element_text(size=12, colour = "darkgrey"),
                   strip.background=ggplot2::element_blank(),
                   axis.title.x=ggplot2::element_blank(),
                   axis.text.x=ggplot2::element_text(size = 8, angle=45, hjust = 1),
                   axis.ticks.x=ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_line(color = "grey", size = 0.5),
                   axis.line.x = ggplot2::element_line(color = "grey", size = 0.1),
                   axis.line.y = ggplot2::element_line(size = 0.5, color = "grey"),
                   strip.text.x = ggplot2::element_text(angle = 0, hjust = 0, size = 10),
                   plot.margin = ggplot2::margin(t = 0, r = 5, b = 20, l = 20, unit = "pt"),
                   axis.title.y = ggplot2::element_text(size = 10, margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0))
    ) +
    ggplot2::scale_color_manual(values = pal) +
    ggplot2::scale_alpha(range = c(0, 0.8), limits = c(0, 0.8)) +
    ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(ymin, ymax),
                                labels = ticklabels,
                                breaks = tickvals) +
    ggplot2::scale_x_continuous(expand = c(0, 0), limits = c(0, 210),
                                breaks = (xScale(starts) + xScale(starts + widthscale))/2,
                                labels = names(widthscale))
  
  for (s in names(widthscale)){
    xstart = xScale(starts[s])
    xend = xScale(starts[s] + widthscale[s])
    p <- p + ggplot2::annotate("segment", x = xstart, xend = xend, y = -1, yend = -1,
                               size = 3, colour = pal[s])
  }
  
  if (capped){
    p <- p + ggplot2::annotate(geom = "text", x = 180,
                               y = 16.2, label = "values above this threshold are capped", size = 2, color = "grey") +
      ggplot2::geom_hline(yintercept = 16, linetype = "dashed", size = 0.2, color = "grey")
  }
  
  if (interactive){
    p <- p %>% plotly::ggplotly(tooltip = "text")
    p <- p %>% plotly::highlight(on = "plotly_click", off = "plotly_doubleclick", dynamic = FALSE, persistent = FALSE)
  }
  
  return(p)
}

GoFunction <- function(Memberlist){
  # for member element, separate cluster and gene name vectors.
  GeneNames = names(Memberlist)
  
  Cluster = as.numeric(Memberlist)
  
  # generate a list of list and name element based on clustering
  gProfilerInput <- split(GeneNames, Cluster)
  names(gProfilerInput) <- paste("C", names(gProfilerInput), sep = "")
  
  # remove single element clusters
  gProfilerInput <- gProfilerInput[lengths(gProfilerInput) > 1]
  
  # now we are ready to load the GOST objct
  multiGOSTResult <- gost(query = gProfilerInput,
                          multi_query = TRUE,
                          significant = TRUE,
                          exclude_iea = TRUE,
                          correction_method = "fdr"
                          )
  
  nClusters = length(unique(Cluster))
  OutputList = list(multiGOSTResult, nClusters)
  
  return(OutputList)
}

publish_gosttable <- function(gostres, highlight_terms = NULL, use_colors = TRUE, show_columns = c("source", "term_name", "term_size", "intersection_size"), filename = NULL, ggplot = NULL, dataframe = TRUE){
  # gostres is the GOSt response list (contains results and metadata) or a data frame
  term_id <- p_values <- query <- p_value <- NULL
  
  if (class(gostres) == "list"){
    if (!("result" %in% names(gostres))) stop("Name 'result' not found from the input")
    df <- gostres$result
  } else if (class(gostres) == "data.frame"){
    df <- gostres
  } else {
    stop("The input 'gostres' should be a data frame or a list from the gost() function.")
  }
  
  # make sure all the essential column names are there
  if (!"term_id" %in% colnames(df)) stop("The required column 'term_id' is missing from the input")
  if (!any(grepl("p_value", colnames(df)))) stop("Column 'p_value(s)' is missing from the input")
  
  # selected terms
  if (is.null(highlight_terms)) {
    # show full table if no terms given
    highlight_terms = df
  }
  
  if (is.data.frame(highlight_terms)){
    message("The input 'highlight_terms' is a data.frame. The column 'term_id' will be used.")
    if ("term_id" %in% colnames(highlight_terms)){
      highlight_terms <- highlight_terms$term_id
    }
    else{
      stop("No column named 'term_id'.")
    }
  }
  
  subdf <- base::subset(df, term_id %in% highlight_terms)
  
  if (nrow(subdf) == 0){
    stop("None of the term IDs in the 'highlight_terms' were found from the results.")
  }
  
  highlight_terms <- unique(highlight_terms)
  subdf$id <- match(subdf$term_id, highlight_terms)
  
  # order by id column
  subdf = subdf[order(subdf$id),]
  
  # default column names to show
  show_columns <- unique(append(show_columns, c("id", "term_id", "p_value")))
  gp_colnames <- c("id", "source", "term_id", "term_name", "term_size", "query_size", "intersection_size", "p_value", "intersection_sizes", "query_sizes")
  
  colnames <- gp_colnames[which(gp_colnames %in% show_columns)]
  
  # include non gprofiler columns
  if (length(setdiff(show_columns, gp_colnames)) > 0){
    colnames <- append(colnames, setdiff(show_columns, gp_colnames))
  }
  
  # If multiquery
  if ("p_values" %in% colnames(subdf)){
    if ("meta" %in% names(gostres)){
      meta <- gostres$meta
      subdf$query <- list(names(meta$query_metadata$queries))
    } else {
      qnames = paste("query", seq(1, length(subdf$p_values[[1]])), sep = "_")
      subdf$query <- list(qnames)
    }
    spread_col = c("p_values")
    if ("query_sizes" %in% show_columns){
      spread_col = append(spread_col, "query_sizes")
    }
    if ("intersection_sizes" %in% show_columns){
      spread_col = append(spread_col, "intersection_sizes")
    }
    # spread the data frame to correct form
    subdf <- tidyr::unnest(data = subdf, cols = c(spread_col, query))
    subdf <- dplyr::rename(subdf, p_value = p_values)
    subdf$p_value <- formatC(subdf$p_value, format = "e", digits = 1)
    showdf <- subdf[,stats::na.omit(match(c(colnames, "query"), names(subdf)))]
    showdf <- tidyr::pivot_wider(showdf, names_from = query, values_from = c(p_value, spread_col[spread_col!="p_values"]), names_prefix = ifelse(length(spread_col) == 1,"p_value ", ""))
    
  } else {
    if ("query" %in% names(subdf) & length(unique(subdf$query)) > 1){
      subdf$p_value <- formatC(subdf$p_value, format = "e", digits = 1)
      showdf <- subdf[,stats::na.omit(match(c(colnames, "query"), names(subdf)))]
      spread_col <- c("p_value", "intersection_size", "query_size")
      spread_col <- intersect(colnames(showdf), spread_col)
      spread_col <- intersect(show_columns, spread_col)
      showdf <- tidyr::pivot_wider(showdf, names_from = query, values_from = spread_col, names_prefix = ifelse(length(spread_col) == 1,"p_value ", ""))
      # order columns by query
      if ('meta' %in% names(gostres)){
        input_order <- names(gostres$meta$query_metadata$queries)
        if (length(spread_col) == 1){
          input_order <- paste("p_value", input_order)
        } else{
          input_order <- unlist(lapply(spread_col, function(x) paste(x, input_order, sep = "_")))
        }
        
        showdf <- showdf[c(names(showdf)[stats::na.omit(match(colnames, names(showdf)))], input_order) ]
        
      }
      
    } else {
      subdf$p_value <- formatC(subdf$p_value, format = "e", digits = 1)
      showdf <- subdf[,stats::na.omit(match(colnames, names(subdf)))]
    }
  }
  if(dataframe == TRUE){
    return(showdf)
  }
  
  # find the columns to color
  idx <- which(grepl(pattern = "p_value", x = names(showdf)))
  
  # Prepare table
  if (use_colors){
    order_of_cl = names(showdf)[idx]
    # add empty columns next to pvalue columns that show the color scale (to the right)
    showdf[paste0(1, order_of_cl)] <- NA
    # update the order with extra pvalue color code column
    order_of_cl2 = c(rbind(paste0(1, order_of_cl), order_of_cl))
    showdf = showdf[,c(names(showdf)[1:min(idx)-1], order_of_cl2)]
    colours <- matrix("white", nrow(showdf), ncol(showdf))
    # add colors to empty columns
    temp_df = showdf[, order_of_cl2, drop = F]
    temp_cols <- sapply(temp_df, function(x) ifelse(!is.na(x), mapViridis(-log10(as.numeric(x))), "white"))
    if (nrow(temp_df) == 1){
      temp_cols = data.frame(t(temp_cols), check.names = F, stringsAsFactors = F)
    }
    # switch values
    temp_cols[,seq(1,ncol(temp_cols),2)] = temp_cols[,seq(2,ncol(temp_cols),2)]
    temp_cols[,seq(2,ncol(temp_cols),2)] = "white"
    colours[,which(names(showdf) %in% order_of_cl2)] <- temp_cols
    if (nrow(temp_df) == 1){
      colours = unlist(colours)
    }
    # remove column names from color scale column
    showdf[,startsWith(names(showdf), "1")] = ""
    # rename the column
    names(showdf)[startsWith(names(showdf), "1")] = ""
    
  } else {
    colours <- matrix("white", nrow(showdf), ncol(showdf))
  }
  
  fontcolours <- matrix("black", nrow(showdf), ncol(showdf))
  fontfaces <- matrix("plain", nrow(showdf), ncol(showdf))
  
  th <- gridExtra::ttheme_default(base_size = 10,
                                  padding = grid::unit(c(4, 4), "mm"),
                                  core=list(
                                    padding.h = grid::unit(c(15,15), "mm"),
                                    padding.v = grid::unit(c(15,15), "mm"),
                                    bg_params = list(fill = colours, col="black", lwd = 0.5),
                                    fg_params=list(hjust = 0, x = 0.02, col=fontcolours, fontface=fontfaces)),
                                  colhead=list(bg_params = list(fill = "gray99", lwd = 0.5, col = "black"),
                                               fg_params=list(col="gray39", fontface="bold")),
                                  rowhead=list(fg_params=list(col="black", fontface="bold")))
  
  tb <- gridExtra::tableGrob(showdf, theme = th, rows = NULL)
  h <- grid::unit.c(sum(tb$heights))
  w <- grid::unit.c(sum(tb$widths))
  tg <- gridExtra::arrangeGrob(tb, ncol = 1, widths = w, heights = h, bottom = grid::textGrob("g:Profiler (biit.cs.ut.ee/gprofiler)", x = 0.95, hjust = 1, gp = grid::gpar(fontsize=10, font=8, col = "cornflowerblue")))
  
  if(ggplot){
    p <- ggplot2::ggplot() + ggplot2::annotation_custom(tg) + ggplot2::geom_blank() + ggplot2::theme_void()
  }
  p <- ggplot2::ggplot() + ggplot2::annotation_custom(tg) + ggplot2::geom_blank() + ggplot2::theme_void()
  
  if (is.null(filename)){
    if (ggplot){
      p <- ggplot2::ggplot() + ggplot2::annotation_custom(tg) + ggplot2::geom_blank() + ggplot2::theme_void()
      return(p)
    } else {
      return(tg)
    }
    
  } else{
    imgtype <- strsplit(basename(filename), split="\\.")[[1]][-1]
    
    if (length(imgtype) == 0) {
      filename = paste0(filename, ".pdf")
    }
    
    if (tolower(imgtype) %in% c("png", "pdf", "jpeg", "tiff", "bmp")){
      width = grid::convertWidth(sum(tg$widths), "in", TRUE) + 0.2
      height = grid::convertHeight(sum(tg$heights), "in", TRUE) + 0.2
      p <- ggplot2::ggplot() + ggplot2::annotation_custom(tg) + ggplot2::geom_blank() + ggplot2::theme_void()
      ggplot2::ggsave(filename = filename, plot = p, height = height, width = width)
      message("The image is saved to ", filename)
      return(p)
    } else {
      stop("The given file format is not supported.\nPlease use one of the following extensions: .png, .pdf, .jpeg, .tiff, .bmp")
    }
  }
}

#?imageOutput
# testing
#Test dependencies
# library(tidyr)
# library(dplyr)
# library(data.table)
# library(igraph)
# library(networkD3)
# library(ggraph)
# library(htmlwidgets)
# library(gprofiler2)

# setwd("C:/Users/jeppe/Desktop/InDepthMap")

# gList = c("USP18", "KRAS", "TP53", "SLC9A1")
# 
# CorrMatrix = readRDS("data/CoDepMapCorMatrixFastCorr.rds")
# 
# Correlation = "Positive"
# 
# NumberOfCoDeps = 10
# 
# Threshold = 0.1
# 
# net <- NetworkMapper(gList, CorrMatrix, Correlation, NumberOfCoDeps, Threshold)
# clusterList <- net[[2]]

# 
# GOSTplot <- GoFunction(clusterList)
# 
# resdf <- publish_gosttable(GOSTplot[[1]])
# 
# 
# res <- GOSTplot[[1]]$result
# 
# # If multiquery
# if ("p_values" %in% colnames(subdf)){
#   if ("meta" %in% names(gostres)){
#     meta <- GOSTplot[[1]]$meta
#     res$query <- list(names(meta$query_metadata$queries))
#   } else {
#     qnames = paste("query", seq(1, length(res$p_values[[1]])), sep = "_")
#     res$query <- list(qnames)
#   }
#   
#   spread_col = c("p_values", "query_sizes", "intersection_sizes")
#   
#   # spread the data frame to correct form
#   res <- tidyr::unnest(data = res, cols = c(spread_col, query))
#     
#   res <- dplyr::rename(res, p_value = p_values)
#   res$p_value <- formatC(res$p_value, format = "e", digits = 1)
#   showdf <- res[,stats::na.omit(match(c(colnames, "query"), names(res)))]
#   showdf <- tidyr::pivot_wider(res, names_from = query, values_from = c(p_value, spread_col[spread_col!="p_values"]), names_prefix = ifelse(length(spread_col) == 1,"p_value ", ""))
#   
# }


