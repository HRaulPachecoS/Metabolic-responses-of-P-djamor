#========================================================
# SCRIPT CIENTÍFICO
# Metabolomic responses of Pleurotus djamor
#========================================================

rm(list = ls())
graphics.off()

#========================================================
# 1. LIBRERÍAS
#========================================================

library(readxl)
library(dplyr)
library(tidyr)
library(stringr)

library(ggplot2)
library(ggrepel)
library(viridis)
library(patchwork)
library(cowplot)
library(ragg)

library(FactoMineR)
library(factoextra)

library(vegan)
library(pheatmap)

library(magick)


#========================================================
# 2. IMPORTACIÓN DE DATOS
#========================================================

mt <- read_excel("Data/Raw/ProductionPdjamor.xlsx")

df_raw <- data.frame(mt)

#========================================================
# 3. LIMPIEZA DE METABOLITOS
#========================================================

df <- df_raw %>%

  filter(Qual >= 70) %>%

  mutate(

    SampleID = paste(
      Substrate,
      percentage,
      Sample,
      sep = "_"
    )

  ) %>%

  rename(

    Percentage = percentage,

    Rep = Sample,

    PK_TY = `PK..TY`,

    Pct_Max = `Pct.Max`,

    Pct_Total = `Pct.Total`,

    Area_Pct = `Area.Pct`,

    Compound = `Library.ID`

  )

clean_metabolite_names <- function(x){

  x <- gsub("\\.beta\\.","beta",x)

  x <- gsub("\\.alpha\\.","alpha",x)

  x <- gsub(",.*","",x)

  x <- gsub(
    " (bis|tris|tetrakis|pentakis|hexakis|octakis).*",
    "",
    x
  )

  x <- gsub(
    " (trimethylsilyl|methyl ester|TMS).*",
    "",
    x,
    ignore.case=TRUE
  )

  trimws(x)

}

df <- df %>%

mutate(

Compound = clean_metabolite_names(Compound),

Compound = str_to_title(Compound)

)


#========================================================
# 4. LIMPIEZA DE METABOLITOS
#========================================================

df <- df %>%

  mutate(

    Compound = clean_metabolite_names(Compound),

    Compound = str_to_title(Compound)

  )

#--------------------------------------------------------
# Correcciones manuales de nombres
#--------------------------------------------------------

df <- df %>%

  mutate(

    Compound = case_when(

      grepl("7H-purine", Compound, ignore.case = TRUE) ~ "Guanine",

      Compound == "7-(Trimethylsilyl)-2" ~ "Guanine",

      grepl("Adenosine", Compound, ignore.case = TRUE) ~ "Adenosine",

      Compound == "L-Isoleucine" ~ "L-Isoleucine",

      Compound == "Cis-7" ~ "Adrenic Acid",

      TRUE ~ Compound

    )

  )

#--------------------------------------------------------
# Eliminación de contaminantes
#--------------------------------------------------------

contaminantes <- c(

  "phthalate",

  "siloxane",

  "silane",

  "adipate",

  "ethoxy",

  "silanol",

  "trisiloxane",

  "silanamine",

  "chloro-n",

  "bromazepam",

  "unknown",

  "unidentified"

)

df <- df %>%

  filter(

    nchar(Compound) > 2,

    !grepl("^[0-9]+$", Compound),

    Compound != "",

    !grepl(

      paste(contaminantes, collapse="|"),

      Compound,

      ignore.case = TRUE

    )

  )

#========================================================
# Consolidación de metabolitos
#========================================================

df_final <- df %>%

  group_by(

    SampleID,

    Substrate,

    Percentage,

    Rep,

    Compound

  ) %>%

  summarise(

    Area_Pct = sum(

      Area_Pct,

      na.rm = TRUE

    ),

    .groups = "drop"

  )


#========================================================
# Corrección de nombres de sustrato
#========================================================

df_final$Substrate <-

  gsub(

    "\\s*\\(.*\\)",

    "",

    df_final$Substrate

  )

df_final$Substrate[

  df_final$Substrate == "Coir"

] <- "Coconut fiber"

df_final$Substrate[

  df_final$Substrate == "Coffee husk"

] <- "Coffee husk"

df_final$Substrate[

  df_final$Substrate == "Barley"

] <- "Barley Straw"


#========================================================
# Conversión de variables
#========================================================

df_final <- df_final %>%

  mutate(

    Substrate = factor(Substrate),

    Percentage = as.numeric(as.character(Percentage)),

    Rep = factor(Rep)

  )


#========================================================
# Metadatos
#========================================================

meta <- df_final %>%

  distinct(

    SampleID,

    Substrate,

    Percentage,

    Rep

  ) %>%

  arrange(

    Substrate,

    Percentage,

    Rep

  )

order_samples <- meta$SampleID


#========================================================
# Función para construir matriz metabolómica
#========================================================

make_matrix <- function(data){

  mat <- data %>%

    group_by(
      SampleID,
      Compound
    ) %>%

    summarise(
      Area = mean(Area_Pct, na.rm = TRUE),
      .groups = "drop"
    ) %>%

    pivot_wider(

      names_from = Compound,

      values_from = Area,

      values_fill = 0

    )

  mat <- as.data.frame(mat)

  rownames(mat) <- mat$SampleID

  mat$SampleID <- NULL

  mat[] <- lapply(mat, as.numeric)

  as.matrix(mat)

}


#========================================================
# Matriz metabolómica completa
#========================================================

X <- make_matrix(df_final)

#========================================================
# Filtrado por prevalencia
#========================================================

keep_prev <-

  colSums(

    X > 0

  ) >= 5

X_prev <-

  X[ , keep_prev]

#========================================================
# Filtrado por varianza
#========================================================

vars <-

  apply(

    X_prev,

    2,

    var,

    na.rm = TRUE

  )

threshold <-

  quantile(

    vars,

    0.25

  )

X_filt <-

  X_prev[ ,

    vars > threshold

  ]


#========================================================
# Escalado
#========================================================

X_filt_scaled <-

  as.data.frame(

    scale(X_filt)

  )


#========================================================
# Verificación
#========================================================

cat("\n")

cat("----------------------------------\n")

cat("Número de muestras:",

nrow(X_filt_scaled),

"\n")

cat("Número de metabolitos:",

ncol(X_filt_scaled),

"\n")

cat("----------------------------------\n")


#========================================================
# Matriz final
#========================================================

X_analysis <- X_filt_scaled

library(openxlsx)

write.xlsx(as.data.frame(X_analysis),
           "Matriz_metabolomica.xlsx",
           rowNames = TRUE)

samples_after_cleaning <- df %>%
  distinct(
    SampleID = paste(Substrate, Percentage, Rep, sep = "_")
  )

cat("Muestras después de limpieza:",
    nrow(samples_after_cleaning),
    "\n")

#========================================================
# Preparación de metadatos
#========================================================

meta_perm <- meta %>%

  mutate(

    Percentage = ifelse(

      Substrate == "Barley Straw",

      0,

      Percentage

    )

  )

run_permanova <- function(

    X,

    metadata,

    substrate_name){

  idx <- metadata$Substrate %in%

    c(

      substrate_name,

      "Barley Straw"

    )

  X_sub <- X[idx, ]

  meta_sub <- metadata[idx, ]

  dist_sub <- dist(

    X_sub,

    method = "euclidean"

  )

  set.seed(123)

  adonis <- adonis2(

    dist_sub ~

      Percentage,

    data = meta_sub,

    permutations = 999

  )

  beta <- betadisper(

    dist_sub,

    as.factor(meta_sub$Percentage)

  )

  list(

    permanova = adonis,

    betadisper = beta

  )

}

perm_coco <- run_permanova(

  X = X_analysis,

  metadata = meta_perm,

  substrate_name = "Coconut fiber"

)

perm_coco$permanova

anova(

  perm_coco$betadisper

)


perm_coffee <- run_permanova(

  X = X_analysis,

  metadata = meta_perm,

  substrate_name = "Coffee husk"

)

perm_coffee$permanova

anova(

  perm_coffee$betadisper

)


cat("\n=====================================\n")
cat("PERMANOVA - Coconut fiber\n")
print(perm_coco$permanova)

cat("\nBetadisper\n")
print(anova(perm_coco$betadisper))

cat("\n=====================================\n")
cat("PERMANOVA - Coffee husk\n")
print(perm_coffee$permanova)

cat("\nBetadisper\n")
print(anova(perm_coffee$betadisper))




#========================================================
# Función PCA
#========================================================

run_pca <- function(X, metadata, substrate_name){

  idx <- metadata$Substrate %in%
    c(substrate_name, "Barley Straw")

  X_sub <- X[idx, ]

  meta_sub <- metadata[idx, ]

  pca <- prcomp(
    X_sub,
    center = FALSE,
    scale. = FALSE
  )

  list(
    pca = pca,
    X = X_sub,
    meta = meta_sub
  )

}

#========================================================
# PCA - Coconut fiber
#========================================================

meta_coco <- meta %>%
  filter(
    Substrate %in%
      c("Barley Straw",
        "Coconut fiber")
  )

X_coco <- X_analysis[
  rownames(X_analysis) %in%
    meta_coco$SampleID,
]

X_coco <- X_coco[
  match(
    meta_coco$SampleID,
    rownames(X_coco)
  ),
]

res.pca.coco <- PCA(
  X_coco,
  graph = FALSE
)

#para el caffe

meta_coffee <- meta %>%
  filter(
    Substrate %in%
      c("Barley Straw",
        "Coffee husk")
  )

X_coffee <- X_analysis[
  rownames(X_analysis) %in%
    meta_coffee$SampleID,
]

X_coffee <- X_coffee[
  match(
    meta_coffee$SampleID,
    rownames(X_coffee)
  ),
]

res.pca.coffee <- PCA(
  X_coffee,
  graph = FALSE
)

loadings_coco <-

as.data.frame(
  res.pca.coco$var$coord
)

top20_loadings_coco <-

loadings_coco %>%

mutate(

Compound=rownames(loadings_coco),

importance=

abs(Dim.1)+

abs(Dim.2)

) %>%

arrange(

desc(importance)

) %>%

slice(1:20)

top20_names_coco <-

top20_loadings_coco$Compound

top20_names_coco


loadings_coffee <- as.data.frame(
  res.pca.coffee$var$coord
)

top20_loadings_coffee <- loadings_coffee %>%
  mutate(
    Compound = rownames(loadings_coffee),
    importance = abs(Dim.1) + abs(Dim.2)
  ) %>%
  arrange(
    desc(importance)
  ) %>%
  slice(1:20)

top20_names_coffee <- top20_loadings_coffee$Compound

top20_names_coffee


#matrices para mantel


#========================================================
# MATRICES - COCONUT FIBER
#========================================================

X_coco_all <- make_matrix(

  df_final %>%

    filter(

      Substrate %in%

        c("Barley Straw",
          "Coconut fiber")

    )

)

X_coco_top20 <- make_matrix(

  df_final %>%

    filter(

      Compound %in% top20_names_coco,

      Substrate %in%

        c("Barley Straw",
          "Coconut fiber")

    )

)

#========================================================
# MATRICES - COFFEE HUSK
#========================================================

X_coffee_all <- make_matrix(

  df_final %>%

    filter(

      Substrate %in%

        c("Barley Straw",
          "Coffee husk")

    )

)

X_coffee_top20 <- make_matrix(

  df_final %>%

    filter(

      Compound %in% top20_names_coffee,

      Substrate %in%

        c("Barley Straw",
          "Coffee husk")

    )

)

#reordenar muestras

X_coco_all <-

X_coco_all[

match(

meta_coco$SampleID,

rownames(X_coco_all)

),

]

X_coco_top20 <-

X_coco_top20[

match(

meta_coco$SampleID,

rownames(X_coco_top20)

),

]


X_coffee_all <-

X_coffee_all[

match(

meta_coffee$SampleID,

rownames(X_coffee_all)

),

]

X_coffee_top20 <-

X_coffee_top20[

match(

meta_coffee$SampleID,

rownames(X_coffee_top20)

),

]

#escalar

X_coco_all <- scale(X_coco_all)

X_coco_top20 <- scale(X_coco_top20)

X_coffee_all <- scale(X_coffee_all)

X_coffee_top20 <- scale(X_coffee_top20)


#dimensiones
dim(X_coco_all)
dim(X_coco_top20)
dim(X_coffee_all)
dim(X_coffee_top20)


#========================================================
# MANTEL TEST
#========================================================

library(vegan)

set.seed(123)

mantel_coco <- mantel(
  dist(scale(X_coco_all)),
  dist(scale(X_coco_top20)),
  method = "pearson",
  permutations = 9999
)

set.seed(123)

mantel_coffee <- mantel(
  dist(scale(X_coffee_all)),
  dist(scale(X_coffee_top20)),
  method = "pearson",
  permutations = 9999
)
