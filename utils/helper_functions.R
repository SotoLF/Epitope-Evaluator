# helper_functions.R - Funciones auxiliares para Epitope-Evaluator
#
# Este archivo contiene funciones utilitarias usadas en toda la aplicación
# para estandarizar operaciones comunes.

#' Estandariza nombres de columnas en un DataFrame
#'
#' @param df DataFrame a estandarizar
#' @return DataFrame con nombres de columnas estandarizados
standardize_column_names <- function(df) {
  # Asegura que las primeras cuatro columnas tengan nombres consistentes
  if (ncol(df) >= 4) {
    colnames(df)[1:4] <- c("Peptide", "Pos", "Length", "ID")
  }
  
  # Convierte columnas numéricas a numéricas
  if (ncol(df) > 4) {
    df[, 5:ncol(df)] <- sapply(df[, 5:ncol(df)], as.numeric)
  }
  
  return(df)
}

#' Redondea columnas numéricas a un número especificado de decimales
#'
#' @param df DataFrame con columnas a redondear
#' @param digits Número de decimales (por defecto 2)
#' @return DataFrame con columnas numéricas redondeadas
round_numeric_columns <- function(df, digits = 2) {
  # Encuentra columnas numéricas
  numeric_cols <- sapply(df, is.numeric)
  
  # Redondea columnas numéricas
  if (any(numeric_cols)) {
    df[, numeric_cols] <- round(df[, numeric_cols], digits)
  }
  
  return(df)
}

#' Valida archivos de entrada y configuración
#'
#' @param input_file Ruta al archivo de predicción
#' @param fasta_file Ruta al archivo FASTA
#' @param cond_pred Tipo de predictor
#' @param method Método de puntuación
#' @return NULL si la validación pasa, mensaje de error en caso contrario
validate_input_files <- function(input_file, fasta_file, cond_pred, method) {
  # Verifica si el archivo de salida está vacío
  if (is.null(input_file)) {
    return("Need Output Prediction File.")
  }
  
  # Verifica si el archivo FASTA está vacío
  if (is.null(fasta_file)) {
    return("Need Fasta File.")
  }
  
  # Verifica si el FASTA es correcto
  if (is.na(tryCatch(read.fasta(fasta_file), error = function(e) return(NA)))) {
    return('File submitted is not in Fasta format (See Documentation)')
  }
  
  # Verifica si la tabla de predicción es correcta según el tipo de predictor
  if (cond_pred == "NetMHC") {
    if (is.na(tryCatch(Parse_NetMHC(input_file, fasta_file, method), error = function(e) return(NA)))) {
      return('Prediction File does not correspond to a *NetMHC* Output format (See Documentation)')
    }
  } else if (cond_pred == "NetMHCpan") {
    if (is.na(tryCatch(Parse_NetMHCPAN(input_file, fasta_file, method), error = function(e) return(NA)))) {
      return('Prediction File does not correspond to a *NetMHCPan* Output format (See Documentation)')
    }
  } else if (cond_pred == "NetMHCIIpan") {
    if (is.na(tryCatch(Parse_NetMHCIIPAN(input_file, fasta_file, method), error = function(e) return(NA)))) {
      return('Prediction File does not correspond to a *NetMHCIIPan* Output format (See Documentation)')
    }
  } else if (cond_pred == "MHCFlurry") {
    if (is.na(tryCatch(Parse_MHCFlurry(input_file, fasta_file, method), error = function(e) return(NA)))) {
      return('Prediction File does not correspond to a *MHCFlurry* Output format (See Documentation)')
    }
  } else if (cond_pred == "IEDB Consensus") {
    if (is.na(tryCatch(Parse_IEDB(input_file, fasta_file, method), error = function(e) return(NA)))) {
      return('Prediction File does not correspond to a *IEDB-Consensus* Output format (See Documentation)')
    }
  } else if (cond_pred == "Other") {
    if (is.na(tryCatch(read.delim(input_file, sep = '\t') %>% as.data.frame(), error = function(e) return(NA)))) {
      return('Prediction File does not correspond to the *specified* Output format (See Documentation)')
    }
  }
  
  return(NULL)  # Devuelve NULL si la validación pasa
}

#' Actualiza inputs UI basados en datos analizados
#'
#' @param session Sesión Shiny
#' @param parsed_data DataFrame con datos analizados
#' @return NULL
update_ui_inputs <- function(session, parsed_data) {
  # Actualiza selección de alelos
  alleles <- colnames(parsed_data)[c(-1, -2, -3, -4)]
  updateCheckboxGroupInput(session, "allele", choices = alleles, selected = alleles[1])
  updateCheckboxGroupInput(session, "allele_4", choices = alleles, selected = alleles[1])
  updateCheckboxGroupInput(session, "allele_6", choices = alleles, selected = alleles[1])
  updateCheckboxGroupInput(session, "allele_7", choices = alleles, selected = alleles[1])
  
  # Verifica el tipo de alelo para establecer valores predeterminados apropiados
  allele_example <- colnames(parsed_data)[5]
  if (startsWith(toupper(allele_example), "HLA")) {
    # Valores para alelos HLA (MHC Clase I)
    updateNumericInput(session, "xmin", min = 0, max = 2, value = 0)
    updateNumericInput(session, "xmax", min = 0, max = 2, value = 2)
    updateNumericInput(session, "cutoff", min = 0, max = 2, value = 2)
    updateNumericInput(session, "cutoff_4", min = 0, max = 2, value = 2)
    updateNumericInput(session, "sb_ct", min = 0, max = 2, value = 0.5)
    updateNumericInput(session, "wb_ct", min = 0, max = 10, value = 2)
    updateNumericInput(session, "cutoff_6", min = 0, max = 10, value = 2)
    updateNumericInput(session, "cutoff_7", min = 0, max = 10, value = 2)
  } else {
    # Valores para otros alelos (MHC Clase II)
    updateNumericInput(session, "xmin", min = 0, max = 10, value = 0)
    updateNumericInput(session, "xmax", min = 0, max = 10, value = 10)
    updateNumericInput(session, "cutoff", min = 0, max = 10, value = 10)
    updateNumericInput(session, "cutoff_4", min = 0, max = 10, value = 10)
    updateNumericInput(session, "sb_ct", min = 0, max = 2, value = 2)
    updateNumericInput(session, "wb_ct", min = 0, max = 10, value = 10)
    updateNumericInput(session, "cutoff_6", min = 0, max = 10, value = 10)
    updateNumericInput(session, "cutoff_7", min = 0, max = 10, value = 10)
  }
  
  # Actualiza selección de proteínas
  proteins <- parsed_data %>% dplyr::select(ID) %>% unlist() %>% unique()
  updateNumericInput(session, "n_prots", min = 1, max = length(proteins), value = 2)
  updateSelectInput(session, "protein_4", choices = proteins, selected = proteins[1])
  updateCheckboxGroupInput(session, "protein_7", choices = proteins, selected = proteins[1:min(3, length(proteins))])
  
  # Actualiza número mínimo de alelos
  updateNumericInput(session, "min_al",
                     max = length(alleles),
                     value = length(alleles) - 1)
}

#' Calcula datos de puntuación para análisis de distribución
#'
#' @param parsed_data DataFrame con datos analizados
#' @param allele Vector de nombres de alelos seleccionados
#' @param conditional Método de combinación ("Intersection" o "Union")
#' @return DataFrame con datos de puntuación calculados
calculate_score_data <- function(parsed_data, allele, conditional) {
  # Extrae datos de puntuación para los alelos seleccionados
  ep_data.score <- parsed_data %>% 
    dplyr::select(Peptide, all_of(allele)) %>% 
    unique()
  
  # Calcula puntuación combinada según método seleccionado
  if (length(allele) != 1) {
    if (conditional == "Intersection") {
      ep_data.score$plot_score <- apply(ep_data.score[, allele], 1, max)
    } else if (conditional == "Union") {
      ep_data.score$plot_score <- apply(ep_data.score[, allele], 1, min)
    }
  } else {
    ep_data.score$plot_score <- ep_data.score[, allele]
  }
  
  return(ep_data.score)
}
