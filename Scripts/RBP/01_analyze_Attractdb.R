#'---
#' title: Analyze [ATtRACT](https://attract.cnic.es/) database
#' author: Avsec
#' wb:
#'  input: ["data/public/attract/ATtRACT_db.txt",
#'          "data/public/attract/pwm.txt"]
#'---
#'
#' ## Goals
#'
#' - explore the RBP pwm's
#'
#' ## TODO's
#'
#' - Allow to use RBP motifs in cresim simulation framework
#'
#' --------------------------------------------
opts_chunk$set(echo=TRUE, cache=F, message = FALSE,
               fig.width = 8, fig.height = 6
               )
options(width = 80)
#' ## `README.md`
#+ readme, eval = FALSE, echo = TRUE
Gene_name	( no need to explain :-) right?)
Gene_id	( no need to explain :-) right?)
Mutated	(if the target gene is mutated)
Organism	( no need to explain :-) right?)
Motif	( no need to explain :-) right?)
Len	(lenght of the motif)
Experiment_description(when available)
Database (Database from where the motifs were extracted PDB: Protein data bank, C: Cisbp-RNA, R:RBPDB, S: Spliceaid-F, AEDB:ASD)
Pubmed (pubmed ID)
Experiment (type of experiment; short description)
Family (domain)
Matrix_id (linked to the file PWM.txt)
Score (Qscore refer to the paper)

The field Matrix_id refers to the pwm id that you can find in the pwm.txt file.
The position weight matrices are annotated in fasta format.

#' 
#' ## `ATtRACT_db.txt`
#+ fread
dt <- fread("data/public/attract/ATtRACT_db.txt")

## #+ head, results = 'show'
## head(dt)

#'
#' ### Number of organisms
#+ num_organism, fig.width = 8, fig.height = 8
dtc <- dt[, .N, by = Organism][order(N)][, Organism := as.factor_keep_order(Organism)]

ggplot(dtc, aes(x = Organism, y = N)) +
  geom_bar(stat = "identity") +
  tilt_xlab +
  geom_text(aes(y = N + 100, label = as.character(N)))  +
  coord_flip() +
  ylab("Count")

#'
#' ### Number of database
#+ fig.width = 4, fig.height = 4
dtc <- dt[, .N, by = Database][order(N)][, Database := as.factor_keep_order(Database)]
ggplot(dtc, aes(x = Database, y = N)) +
  geom_bar(stat = "identity") +
  tilt_xlab +
  geom_text(aes(y = N + 100, label = as.character(N)))  +
  coord_flip() +
  ylab("Count")

#'
#' ## Number of Families
#+ fig.height = 10
dtc <- dt[, .N, by = Family][order(N)][, Family := as.factor_keep_order(Family)]
ggplot(dtc, aes(x = Family, y = N)) +
  geom_bar(stat = "identity") +
  tilt_xlab +
  geom_text(aes(y = N + 100, label = as.character(N)))  +
  coord_flip() +
  ylab("Count")

#' 
#' ### Number of unique papers: `r dt[, uniqueN(Pubmed)]`
#'
#' ### Only in human we have the mutated cases
dt[, .N, by = .(Organism, Mutated)]

#'
#' ## `pwm.txt`

pwms <- readLines("data/public/attract/pwm.txt")
## 
#+ results = 'show'
pwms %>% head(10)

#'
#'
#' - Number of pwm's: ***`r sum(grepl("^>", pwms))`***
#' - Number of pwm matrices in the table: ***`r dt[, uniqueN(Matrix_id)]`***
#'
#' ### Pwm id's
dt[, table(Matrix_id, useNA = "always")] %>% hist(breaks = 50, main = "Number of table entries per pwm/matrix")

#'
#' ### Sorted ids match the ones from the table

## grep("^>", pwms, value= TRUE) %>% tstrsplit(split = "\t") %>% .[[1]] %>% gsub(">", "", .) %>% sort %>% tail(10)

ids_pwm <- grep("^>", pwms, value= TRUE) %>% tstrsplit(split = "\t") %>% .[[1]] %>% gsub(">", "", .) %>% sort
ids_dt <- dt[, Matrix_id] %>% unique

stopifnot(setequal(ids_pwm, ids_dt))

#'
#' ## Saccharomyces cerevisiae
ORGANISM <- "Saccharomyces_cerevisiae"
#' 
#' - Number of unique pwm's: ***`r dt[Organism == ORGANISM][, uniqueN(Matrix_id)]`***
#' 
#' ### Experiment description
#+ yeast_exper_descr
dtc <- dt[Organism == ORGANISM][, .(N = .N, N_matrix_id = uniqueN(Matrix_id)), by = Experiment_description][order(N)][, Experiment_description := as.factor_keep_order(Experiment_description)]
ggplot(dtc, aes(x = Experiment_description, y = N)) +
  geom_bar(stat = "identity") +
  tilt_xlab +
  geom_text(aes(y = N + 10, label = as.character(N)))  +
  coord_flip() +
  ylab("Count")

ggplot(dtc, aes(x = Experiment_description, y = N_matrix_id)) +
  geom_bar(stat = "identity") +
  tilt_xlab +
  geom_text(aes(y = N_matrix_id + 3, label = as.character(N_matrix_id)))  +
  coord_flip() +
  ylab("N unique matrix id") +
  ggtitle("N unique matrix id")

#'
#' ### Gene names
#+ echo = TRUE, results = 'show'
dt[Organism == ORGANISM][, table(Gene_name)] %>% sort(decr = T)

dt[Organism == ORGANISM][, .N, by = .(Gene_name, Motif)] %>% print(30)
#'
#' ### Motif length
ggplot(dt[Organism == ORGANISM], aes(x = as.factor(Len))) + geom_bar() + xlab("Motif length")

#'
#' ## Homo sapiens
ORGANISM <- "Homo_sapiens"
#' - Number of unique pwm's: ***`r dt[Organism == ORGANISM][, uniqueN(Matrix_id)]`***
#' 
#' ### Experiment description 
#'

dtc <- dt[Organism == ORGANISM][, .(N = .N, N_matrix_id = uniqueN(Matrix_id)), by = Experiment_description][order(N)][N>50][, Experiment_description := as.factor_keep_order(strtrim(Experiment_description, 30))]

ggplot(dtc, aes(x = Experiment_description, y = N)) +
  geom_bar(stat = "identity") +
  tilt_xlab +
  geom_text(aes(y = N + 20, label = as.character(N)))  +
  coord_flip() +
  ylab("Count")

ggplot(dtc, aes(x = Experiment_description, y = N_matrix_id)) +
  geom_bar(stat = "identity") +
  tilt_xlab +
  geom_text(aes(y = N_matrix_id + 3, label = as.character(N_matrix_id)))  +
  coord_flip() +
  ylab("N unique matrix id") +
  ggtitle("N unique matrix id")

#' ### Gene names
#+ echo = TRUE, results = 'show'
dt[Organism == ORGANISM][, table(Gene_name)] %>% sort(decr = T)

dt[Organism == ORGANISM][, .N, by = .(Gene_name, Motif)] %>% print(30)
#'
#' ### Motif length
ggplot(dt[Organism == ORGANISM], aes(x = as.factor(Len))) + geom_bar() + xlab("Motif length")
