

draw_unfiltered_barcode_R <- function(homology_list,
                                      save_path,
                                      edge_type = edge_type,
                                      name = name,
                                      max_x = NULL,
                                      weight_shift = weight_shift,
                                      width = 4,
                                      height_per_panel = 2,
                                      res = 300) {
  
  if (!edge_type %in% c("positive", "negative", "all")) {
    stop("edge_type must be 'positive', 'negative', or 'all'")
  }
  
  all_bars <- list(
    "2" = homology_list[["2"]],
    "1" = homology_list[["1"]],
    "0" = homology_list[["0"]]
  )
  
  colors <- c("#019E4F", "#E94127", "blue")
  
  clean_bars <- function(bars) {
    if (is.null(bars) || length(bars) == 0) {
      return(matrix(numeric(0), ncol = 2))
    }
    
    if (is.vector(bars) && length(bars) == 2) {
      bars <- matrix(bars, ncol = 2, byrow = TRUE)
      return(bars)
    }
    
    bars <- as.matrix(bars)
    
    if (length(bars) == 0 || nrow(bars) == 0) {
      return(matrix(numeric(0), ncol = 2))
    }
    
    if (ncol(bars) != 2) {
      stop("Each homology component must be a 2-column object: [birth, death].")
    }
    
    return(bars)
  }
  
  for (dim_name in names(all_bars)) {
    bars <- clean_bars(all_bars[[dim_name]])
    
    if (nrow(bars) > 0) {
      bars[bars != -1] <- bars[bars != -1] - weight_shift
      
      death_sort <- ifelse(bars[, 2] == -1, Inf, bars[, 2])
      bars <- bars[order(bars[, 1], -death_sort), , drop = FALSE]
    }
    
    all_bars[[dim_name]] <- bars
  }
  
  all_values <- c()
  death_values <- c()
  
  for (bars in all_bars) {
    if (!is.null(bars) && nrow(bars) > 0) {
      birth_values <- bars[, 1]
      death_all <- bars[, 2]
      finite_death <- death_all[death_all != -1]
      
      all_values <- c(all_values, birth_values, finite_death)
      death_values <- c(death_values, finite_death)
    }
  }
  
  if (length(all_values) > 0) {
    min_x_data <- min(all_values[all_values != -100])
    max_x_data <- max(all_values)
  } else {
    min_x_data <- 0
    max_x_data <- 1
  }
  
  if (is.null(max_x)) {
    max_x_plot <- max_x_data
  } else {
    max_x_plot <- max_x
  }
  
  span_death <- max_x_plot - min_x_data
  if (span_death == 0) span_death <- 1
  
  if (edge_type == "positive") {
    xlim <- c(0, max_x_plot+span_death*0.1)
  } else if (edge_type %in% c("negative", "all")) {
    xlim <- c(min_x_data-span_death*0.05,
              max_x_plot + span_death*0.1)
  }
  
  if (any(!is.finite(xlim)) || diff(xlim) <= 0) {
    xlim <- c(min_x_data - 1, max_x_plot + 1)
  }
  
  outfile <- file.path(save_path, paste0(name, "_", edge_type, "_barcode.png"))
  
  png(outfile,
      width = width,
      height = height_per_panel * length(all_bars),
      units = "in",
      res = res)
  
  par(mfrow = c(length(all_bars), 1),
      mar = c(0, 0, 0.5, 0),
      oma = c(4, 8, 0, 1))
  
  for (idx in seq_along(all_bars)) {
    dim_name <- names(all_bars)[idx]
    bars <- all_bars[[idx]]
    col_bar <- colors[idx]
    
    n_bars <- if (!is.null(bars)) nrow(bars) else 0
    max_y <- max(n_bars, 5)
    
    plot(NA, NA,axes=FALSE,
         xlim = xlim,
         ylim = c(0.9, max_y),
         xlab = "",
         ylab = "",
         xaxt = if (idx == length(all_bars)) "s" else "n",
         yaxt = "n",
         cex.axis = 3,
         xaxs = "i")
    if(idx==1){
      mtext(expression(beta[2]), side=2, line=4, las=2,cex = 3)
    }else if(idx==2){
      mtext(expression(beta[1]), side=2, line=4, las=2,cex = 3)
    }else{
      mtext(expression(beta[0]), side=2, line=4, las=2,cex = 3)
    }
    
    
    y_ticks <- pretty(c(1, max_y), n = 3)
    # y_ticks[1] <- 1
    axis(2, at = y_ticks, las = 1,cex.axis = 3)
    
    if (idx == length(all_bars)){
      axis(1, at=seq(xlim[1]+span_death*0.2, xlim[2]-span_death*0.2, length=2), labels=FALSE, tck=-0.02)
      
      mtext(side=1, 
            at=seq(xlim[1]+span_death*0.2, xlim[2]-span_death*0.2, length=2),
            text=sprintf("%0.2f", seq(xlim[1]+span_death*0.2, xlim[2]-span_death*0.2, length=2)),
            line=2,  
            cex=2)
    }
    box(lwd=1)
    
    if (n_bars > 0) {
      for (j in seq_len(n_bars)) {
        birth <- bars[j, 1]
        death <- bars[j, 2]
        y_pos <- j #- 1
        
        if (death == -1) {
          segments(x0 = birth, y0 = y_pos,
                   x1 = max_x_plot + span_death*0.03, y1 = y_pos,
                   col = col_bar, lwd = 1.5)
          
          # arrows(x0 = max_x_plot, y0 = y_pos,
          #        x1 = max_x_plot + 1e-3, y1 = y_pos,
          #        length = 0.08, angle = 20,
          #        col = "black", lwd = 1.5)
          text(max_x_plot + span_death*0.03, y_pos, labels = "▲", cex = 2,srt=270, col = "black")
        } else {
          segments(x0 = birth, y0 = y_pos,
                   x1 = death, y1 = y_pos,
                   col = col_bar, lwd = 1.5)
        }
      }
    }
    
  }
  # mtext("Filtration value", side = 1, outer = TRUE, line = 2.7,cex = 1.2)
  
  # mtext(edge_type,side = 3,line = 0,cex=1.5,outer = TRUE)
  
  dev.off()
}


