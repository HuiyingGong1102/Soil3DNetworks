
split_soil_network_fixed <- function(ii) {
  input_file <- paste("D:/Soil3DNetworks/code/network/links_",ii,".xlsx",sep = "")
  output_dir  <- paste("D:/Soil3DNetworks/code/GLMY/",ii,sep = "") 
  df <- readxl::read_xlsx(input_file)
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  
  required_cols <- c("source", "target", "weight", "edge_type")
  miss_cols <- setdiff(required_cols, colnames(df))
  if (length(miss_cols) > 0) {
    stop(paste(":", paste(miss_cols, collapse = ", ")))
  }
  
  df <- df[, required_cols]
  df$weight <- (df$weight / 5)^10
  
  df$effect <- ifelse(df$edge_type == 1, df$weight,
                      ifelse(df$edge_type == 2, -df$weight, NA))
  
  df <- df[, c("source", "target", "effect", "edge_type")]
  colnames(df) <- c("From","To", "Effect", "edge_type")
  
  get_depth <- function(x) {
    sapply(strsplit(x, "-", fixed = TRUE), function(z) as.integer(z[1]))
  }
  
  source_depth <- get_depth(df$From)
  target_depth <- get_depth(df$To)
  depth_diff <- abs(source_depth - target_depth)
  
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  df$From <- paste0("'", df$From)
  df$To <- paste0("'", df$To)
  
  write.csv(
    df,
    file = file.path(output_dir, paste0("links_", ii, ".csv")),
    row.names = FALSE
  )
  
  for (d in 1:6) {
    sub_df <- df[source_depth == d & target_depth == d, ]
    write.csv(
      sub_df,
      file = file.path(output_dir, paste0("depth_", d, "_within.csv")),
      row.names = FALSE
    )
  }
  
  adjacent_df <- df[depth_diff == 1, ]
  write.csv(
    adjacent_df,
    file = file.path(output_dir, "adjacent.csv"),
    row.names = FALSE
  )
  
  skip1_df <- df[depth_diff == 2, ]
  write.csv(
    skip1_df,
    file = file.path(output_dir, "skip1.csv"),
    row.names = FALSE
  )
  
  skip2_df <- df[depth_diff == 3, ]
  write.csv(
    skip2_df,
    file = file.path(output_dir, "skip2.csv"),
    row.names = FALSE
  )
  
  skip3_df <- df[depth_diff == 4, ]
  write.csv(
    skip3_df,
    file = file.path(output_dir, "skip3.csv"),
    row.names = FALSE
  )
  
  skip4_df <- df[depth_diff == 5, ]
  write.csv(
    skip4_df,
    file = file.path(output_dir, "skip4.csv"),
    row.names = FALSE
  )
  
}







library(jsonlite)
json_data <- fromJSON("D:/Soil3DNetworks/code/cluster/Layer_all/json/k434-1.json")
k <- 434
module<-json_data$max_omega_logi
length(table(module))
table(module)

for (ii in 1:k){
  split_soil_network_fixed(ii)
}



library(jsonlite)

prefixes <- c(
  paste0("depth_", 1:6, "_within"),
  "adjacent", "skip1", "skip2", "skip3", "skip4"
)

# types <- c("all", "positive", "negative")
types <- c("positive", "negative")
var_names <- as.vector(outer(prefixes, types, paste, sep = "_"))

result_list <- setNames(vector("list", length(var_names)), var_names)

for (nm in var_names) {
  result_list[[nm]] <- matrix(numeric(0), ncol = 4)
  colnames(result_list[[nm]]) <- paste0(nm,"_beta", 0:3)
}

read_beta_counts <- function(json_path) {
  tryCatch(
    {
      json <- fromJSON(json_path, simplifyVector = FALSE)
      sapply(1:4, function(k) length(json[[as.character(k)]]))
    },
    error = function(e) {
      rep(0, 4)
    }
  )
}

for (ii in 1:k) {
  file_dir <- file.path("D:/Soil3DNetworks/code/GLMY", as.character(ii))
  
  for (nm in var_names) {
    json_path <- file.path(file_dir, paste0(nm, ".json"))
    beta_counts <- read_beta_counts(json_path)
    result_list[[nm]] <- rbind(result_list[[nm]], beta_counts)
  }
}

list2env(result_list, envir = .GlobalEnv)

result_mat <- do.call(cbind, result_list)
rownames(result_mat) <- paste0("M", 1:k)


