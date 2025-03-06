# parsing_functions.R - Funciones para analizar diferentes formatos de archivos de predicción
#
# Este archivo contiene todas las funciones para analizar (parse) los archivos de
# salida de diferentes predictores de epitopos de células T.

#' Función para obtener miembros de intersección
#'
#' @param x DataFrame con datos numéricos
#' @param ... Columnas para intersección
#' @return DataFrame filtrado con intersección de miembros
get_intersect_members <- function(x, ...) {
  # Selecciona solo columnas numéricas con valores entre 0 y 1
  x <- x[, sapply(x, is.numeric)][, 0 <= colMeans(x[, sapply(x, is.numeric)], na.rm=TRUE) & 
                                      colMeans(x[, sapply(x, is.numeric)], na.rm=TRUE) <= 1]
  n <- names(x)
  x %>% rownames_to_column() -> x
  l <- c(...)
  a <- intersect(names(x), l)
  ar <- vector('list', length(n) + 1)
  ar[[1]] <- x
  i = 2
  
  # Construye condiciones de filtrado
  for (item in n) {
    if (item %in% a) {
      if (class(x[[item]]) == 'integer') {
        ar[[i]] <- paste(item, '>= 1')
        i <- i + 1
      }
    } else {
      if (class(x[[item]]) == 'integer') {
        ar[[i]] <- paste(item, '== 0')
        i <- i + 1
      }
    }
  }
  
  # Aplica filtro y devuelve resultado
  do.call(filter_, ar) %>% column_to_rownames() -> x
  return(x)
}

#' Analiza archivos de predicción NetMHCPAN
#'
#' @param prediction_file Ruta al archivo de predicción
#' @param fasta_file Ruta al archivo FASTA
#' @param score_type Tipo de puntuación ("Rank" o "Score")
#' @return DataFrame con datos analizados
Parse_NetMHCPAN <- function(prediction_file, fasta_file, score_type) {
  # Lee el archivo de predicción
  results <- data.frame(read.table(prediction_file, sep = "\t"))
  
  # Obtiene longitud del archivo FASTA
  faas <- seqinr::read.fasta(fasta_file)
  lengths <- nchar(getSequence(faas, as.string = TRUE)) - 8
  
  # Obtiene nombres del archivo FASTA
  raw_names <- data.frame(names = get.fasta.name(fasta_file))
  metadata_df <- tidyr::separate(raw_names, names, into = c("ID"), sep = " ")
  
  # Crea DataFrame de metadatos
  metadata_df$length <- lengths
  metadata_df$prev <- results[c(-1, -2), 3] %>% unique() %>% unlist()
  
  # Obtiene nombres de alelos para Clase I
  raw_header <- results[1, ]
  alleles_ID <- raw_header[grepl("HLA", raw_header)]
  
  # Agrega nombres a DataFrame de metadatos
  results$ID_2 <- metadata_df$ID[match(results[, 3], metadata_df$prev)]
  results$length <- metadata_df$length[match(results[, 3], metadata_df$prev)]
  
  # Determina índices de columnas según tipo de predicción y puntuación
  if ((ncol(results) - 7) / length(alleles_ID) == 6) {
    if (score_type == "Rank") {
      # Índices para % Ranks Clase I
      col_IDs <- seq(from = 7, by = 6, length.out = length(alleles_ID))
    } else {
      # Índices para concentración nM Clase I
      col_IDs <- seq(from = 9, by = 6, length.out = length(alleles_ID))
    }
  } else { # Predicción sin Score de afinidad de unión (BA)
    if (score_type == "Rank") {
      # Índices para % Ranks Clase I
      col_IDs <- seq(from = 7, by = 4, length.out = length(alleles_ID))
    } else {
      # Índices para % Ranks Clase I (No hay columna BA en este tipo de archivo)
      col_IDs <- seq(from = 7, by = 3, length.out = length(alleles_ID))
    }
  }
  
  # Crea DataFrame analizado
  df_parsed <- results[c(-1, -2), c(2, 1, ncol(results), ncol(results)-1, col_IDs)]
  colnames(df_parsed) <- c("Peptide", "Position", "Length", "ID", alleles_ID)
  df_parsed$Position = as.numeric(df_parsed$Position) + 1
  df_parsed[c(5:ncol(df_parsed))] <- sapply(df_parsed[c(5:ncol(df_parsed))], as.numeric)
  
  # Normaliza nombres de columnas reemplazando caracteres especiales
  colnames(df_parsed) = str_replace_all(colnames(df_parsed), ':', '.') %>% 
                        str_replace_all('-', '.')
  
  return(df_parsed)
}

#' Analiza archivos de predicción NetMHC
#'
#' @param prediction_file Ruta al archivo de predicción
#' @param fasta_file Ruta al archivo FASTA
#' @param score_type Tipo de puntuación ("Rank" o "Score")
#' @return DataFrame con datos analizados
Parse_NetMHC <- function(prediction_file, fasta_file, score_type) {
  # Lee el archivo de predicción
  results <- data.frame(read.table(prediction_file, sep = "\t", row.names = NULL))
  
  # Obtiene longitud del archivo FASTA
  faas <- seqinr::read.fasta(fasta_file)
  lengths <- nchar(getSequence(faas, as.string = TRUE)) - 8
  
  # Obtiene nombres del archivo FASTA
  raw_names <- data.frame(names = get.fasta.name(fasta_file))
  metadata_df <- tidyr::separate(raw_names, names, into = c("ID"), sep = " ")
  
  # Crea DataFrame de metadatos
  metadata_df$length <- lengths
  metadata_df$prev <- results[-1, 3] %>% unique() %>% unlist()
  
  # Obtiene nombres de alelos para Clase I
  raw_header <- colnames(results)
  alleles_ID <- raw_header[grepl("HLA", raw_header)]
  
  # Agrega nombres a DataFrame de metadatos
  results$ID_2 <- metadata_df$ID[match(results[, 3], metadata_df$prev)]
  results$length <- metadata_df$length[match(results[, 3], metadata_df$prev)]
  
  # Determina índices de columnas según tipo de puntuación
  if (score_type == "Rank") {
    # Índices para % Ranks Clase I
    col_IDs <- seq(from = 5, by = 3, length.out = length(alleles_ID))
  } else {
    # Índices para concentración nM Clase I
    col_IDs <- seq(from = 4, by = 3, length.out = length(alleles_ID))
  }
  
  # Crea DataFrame analizado
  df_parsed <- results[-1, c(2, 1, ncol(results), ncol(results)-1, col_IDs)]
  colnames(df_parsed) <- c("Peptide", "Position", "Length", "ID", alleles_ID)
  df_parsed$Position = as.numeric(df_parsed$Position) + 1
  df_parsed[c(5:ncol(df_parsed))] <- sapply(df_parsed[c(5:ncol(df_parsed))], as.numeric)
  
  return(df_parsed)
}

#' Analiza archivos de predicción NetMHCIIPAN
#'
#' @param prediction_file Ruta al archivo de predicción
#' @param fasta_file Ruta al archivo FASTA
#' @param score_type Tipo de puntuación ("Rank" o "Score")
#' @return DataFrame con datos analizados
Parse_NetMHCIIPAN <- function(prediction_file, fasta_file, score_type) {
  # Lee el archivo de predicción
  results <- data.frame(read.table(prediction_file, sep = ",", header = FALSE))
  
  # Obtiene nombres de alelos para Clase II
  alleles_vector = results[1, ] %>% strsplit('\t') %>% unlist()
  alleles_ID = alleles_vector[grepl('D', alleles_vector)]
  
  # Analiza datos
  results = results[c(-1, -2), ] %>% strsplit('\t')
  results = as.data.frame(do.call(rbind, results))
  
  # Obtiene longitud del archivo FASTA
  faas <- seqinr::read.fasta(fasta_file)
  lengths <- nchar(getSequence(faas, as.string = TRUE)) - 8
  
  # Obtiene nombres del archivo FASTA
  raw_names <- data.frame(names = get.fasta.name(fasta_file))
  metadata_df <- tidyr::separate(raw_names, names, into = c("ID"), sep = " ")
  
  # Crea DataFrame de metadatos
  metadata_df$length <- lengths
  metadata_df$prev <- results[, 3] %>% unique() %>% unlist()
  
  # Agrega nombres a DataFrame de metadatos
  results$ID_2 <- metadata_df$ID[match(results[, 3], metadata_df$prev)]
  results$length <- metadata_df$length[match(results[, 3], metadata_df$prev)]
  
  # Determina índices de columnas según tipo de predicción y puntuación
  if ((ncol(results) - 8) / length(alleles_ID) == 6) {
    if (score_type == "Rank") {
      # Índices para % Ranks Clase II
      col_IDs <- seq(from = 7, by = 6, length.out = length(alleles_ID))
    } else {
      # Índices para concentración nM Clase II
      col_IDs <- seq(from = 10, by = 6, length.out = length(alleles_ID))
    }
  } else { # Predicción sin Score de afinidad de unión (BA)
    if (score_type == "Rank") {
      # Índices para % Ranks Clase II
      col_IDs <- seq(from = 7, by = 3, length.out = length(alleles_ID))
    } else {
      # Índices para % Ranks Clase II (No hay columna BA en este tipo de archivo)
      col_IDs <- seq(from = 7, by = 3, length.out = length(alleles_ID))
    }
  }
  
  # Crea DataFrame analizado
  df_parsed <- results[, c(2, 1, ncol(results), ncol(results)-1, col_IDs)]
  colnames(df_parsed) <- c("Peptide", "Position", "Length", "ID", alleles_ID)
  df_parsed$Position = as.numeric(df_parsed$Position)
  df_parsed[c(5:ncol(df_parsed))] <- sapply(df_parsed[c(5:ncol(df_parsed))], as.numeric)
  
  return(df_parsed)
}

#' Analiza archivos de predicción MHCFlurry
#'
#' @param prediction_file Ruta al archivo de predicción
#' @param fasta_file Ruta al archivo FASTA
#' @param score_type Tipo de puntuación ("Rank" o "Score")
#' @return DataFrame con datos analizados
Parse_MHCFlurry <- function(prediction_file, fasta_file, score_type) {
  # Lee el archivo de predicción
  results <- data.frame(read_delim(prediction_file, delim = ",", 
                                  escape_double = FALSE, trim_ws = TRUE))
  
  # Selecciona columnas según tipo de puntuación
  if (score_type == "Rank") {
    # Índices para % Ranks Clase I
    results <- results[, c(3, 2, 1, 8, 9)]
  } else {
    # Índices para concentración nM Clase I
    results <- results[, c(3, 2, 1, 8, 7)]
  }
  
  # Reorganiza resultados
  results <- reshape2::dcast(data = results, formula = peptide + pos + sequence_name ~ best_allele)
  
  # Obtiene longitud del archivo FASTA
  faas <- seqinr::read.fasta(fasta_file)
  lengths <- nchar(getSequence(faas, as.string = TRUE)) - 8
  
  # Obtiene nombres del archivo FASTA
  raw_names <- data.frame(names = get.fasta.name(fasta_file))
  metadata_df <- tidyr::separate(raw_names, names, into = c("ID"), sep = " ")
  
  # Crea DataFrame de metadatos
  metadata_df$length <- lengths
  
  # Agrega nombres a DataFrame de metadatos
  results$length <- metadata_df$length[match(results[, 3], metadata_df$ID)]
  
  # Crea DataFrame analizado
  df_parsed <- results[c(1, 2, ncol(results), 3, 4:(ncol(results)-1))]
  colnames(df_parsed)[1:4] <- c("Peptide", "Position", "Length", "ID")
  df_parsed$Position = as.numeric(df_parsed$Position) + 1
  df_parsed[c(5:ncol(df_parsed))] <- sapply(df_parsed[c(5:ncol(df_parsed))], as.numeric)
  
  return(df_parsed)
}

#' Analiza archivos de predicción IEDB
#'
#' @param prediction_file Ruta al archivo de predicción
#' @param fasta_file Ruta al archivo FASTA
#' @param score_type Tipo de puntuación ("Rank" o "Score")
#' @return DataFrame con datos analizados
Parse_IEDB <- function(prediction_file, fasta_file, score_type) {
  # Lee el archivo de predicción
  results <- read.delim(prediction_file, sep = '\t')
  
  # Selecciona columnas según tipo de puntuación
  if (score_type == "Rank") {
    # Índices para % Ranks Clase I
    results <- results[, c(6, 3, 2, 1, 7)]
  } else {
    # Índices para concentración nM Clase I
    results <- results[, c(6, 3, 2, 1, 8)]
  }
  
  # Reorganiza resultados
  results <- reshape2::dcast(data = results, formula = peptide + start + seq_num ~ allele)
  
  # Obtiene longitud del archivo FASTA
  faas <- seqinr::read.fasta(fasta_file)
  lengths <- nchar(getSequence(faas, as.string = TRUE)) - 8
  
  # Obtiene nombres del archivo FASTA
  raw_names <- data.frame(names = get.fasta.name(fasta_file))
  metadata_df <- tidyr::separate(raw_names, names, into = c("ID"), sep = " ")
  
  # Crea DataFrame de metadatos
  metadata_df$length <- lengths
  
  # Asigna nombres y longitudes
  results$seq_num <- metadata_df[results$seq_num, 1]  
  results$length <- metadata_df$length[match(results$seq_num, metadata_df$ID)]
  
  # Crea DataFrame analizado
  df_parsed <- results[c(1, 2, ncol(results), 3, 4:(ncol(results)-1))]
  colnames(df_parsed)[1:4] <- c("Peptide", "Position", "Length", "ID")
  df_parsed[c(5:ncol(df_parsed))] <- sapply(df_parsed[c(5:ncol(df_parsed))], as.numeric)
  
  return(df_parsed)
}
