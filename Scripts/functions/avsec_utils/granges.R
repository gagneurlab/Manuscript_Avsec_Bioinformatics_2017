##' dependencies:
##' - useful-fun.R
##'     - is_unique, duplicated_values


##' Get motif positions
##' @param set Set of sequences with unique seqnames
##' @param motif Motif vector
##' @param max_mismatch Maximal number of mismatches allowed
##' @param allow_wildcard TRUE to allow the wild-card search (say using N). This actually maps to negated fixed parameter in vmatchPattern
##' @param omit_boundary_matches Should we discard matches at the boundaries (width < motif_width)
##' @param append_sequence Should we append the motif sequnece
##' @param n_cores Number of cores to use (paralelized across motifs)
##' @param ... Additional arguments to vmatchPattern
##' @return Data.table with columns:
##' - seqnames
##' - start
##' - end
##' - width
##' - motif
##' - sequence (optional) actual DNA sequence of the match
##' - query_sequence_length
##' - motif_i (sequential number of the motif occurrence in the sequence)
get_motif_pos <- function(set, motif, max_mismatch = 1, allow_wildcard = TRUE,
                          omit_boundary_matches = TRUE,
                          append_sequence = TRUE, n_cores = 1, ...) {
  mclapply(motif, get_motif_pos_single,
           set = set,
           max_mismatch = max_mismatch,
           allow_wildcard = allow_wildcard,
           omit_boundary_matches = omit_boundary_matches,
           append_sequence = append_sequence,
           mc.cores = n_cores,
           ...
           ) %>% rbindlist(fill = TRUE)
}
get_motif_pos_single <- function(set, motif, max_mismatch = 1, allow_wildcard = TRUE,
                          omit_boundary_matches = TRUE,
                          append_sequence = TRUE, ...) {
  ## input checking
  if (is.null(names(set))) names(set) <- 1:length(set)
  if (!is_unique(names(set))) stop("Sequence names are not unique")

  motif_length <- nchar(motif)
  
  ## compute the motif locations
  res <- vmatchPattern(motif,set,max.mismatch = max_mismatch, fixed = !allow_wildcard, ... )
  res_unlist <- res %>% unlist

  ## old solution:
  ## 
  ## ## get the actual sequence
  ## if (append_sequence == TRUE) {
  ##   sequence <- lapply(names(set), function(sequence) {
  ##     suppressWarnings(as.character(DNAStringSet(Views(set[[sequence]], res[[sequence]]))))
  ##   }) %>% unlist
  ## }

  ## get the actual sequence, concatenate the string:
  if (append_sequence) {
    gr <- res %>% as("GRanges")
    pos_hash <- data.table(seqnames = names(set), pos = c(0, cumsum(width(set))[-length(set)]))
    big_seq <- do.call(xscat, set)
    setkey(pos_hash, seqnames)
    gr$pos <- pos_hash[as.character(seqnames(gr)), pos]
    gr_shifted <- shift(gr, gr$pos)
    sequence <- Views(big_seq, ranges(gr_shifted)) %>% as.character
  }

  motif_positions <- data.table(seqnames = names(res_unlist), start = start(res_unlist), end = end(res_unlist),
                                width = end(res_unlist) - start(res_unlist) + 1
                                )
  ## sequence length
  motif_positions <- merge(motif_positions,
                           data.table(seqnames = names(set), query_sequence_length = width(set)),
                           by = "seqnames", sort = FALSE)

  ## append motif and motif positions
  motif_positions[, motif := motif]
  if (append_sequence == TRUE)  motif_positions[, sequence := sequence]
  
  ## omit edge cases
  if (omit_boundary_matches) {
    motif_positions <- motif_positions[width == motif_length]
    motif_positions <- motif_positions[start > 0]
    motif_positions <- motif_positions[end <= query_sequence_length]
  }

  
  motif_positions[, motif_i := 1:.N, by = seqnames]
  attr(motif_positions, "id_col") <- attr(set, "id_col")
  attr(motif_positions, "seq_col") <- attr(set, "seq_col")
  return(motif_positions)
}

##' append annotation features from the annotattion granges
##' @param gr GRanges for which to append new columns
##' @param anno_gr GRanges with the annotation (GTF file import)
##' @return GRanges with added columns like: is_promoter, is_five_prime_utr, ...
append_anno_features <- function(gr, anno_gr) {
  features <- anno_gr$type %>% unique %>% as.character

  ## append feature from granges
  append_feature <- function(gr, feature_region, feature_name) {
    ol <- findOverlaps(gr, feature_region)
    feature_bases <- unique(queryHits(ol))
    mcols(gr)[[paste0("is_", feature_name)]] <- 1:length(gr) %in% feature_bases
    return(gr)
  }
  ## feature_region <- anno_gr[anno_gr$type == feature ] %>% reduce

  for (feature in features) {
    feature_region <- anno_gr[anno_gr$type == feature ] %>% reduce
    gr <- append_feature(gr, feature_region, feature_name = feature)
  }

  ## append transcripts 
  
  txdb <- makeTxDbFromGRanges(anno_gr)
  introns <- intronsByTranscript(txdb) %>% unlist %>% reduce
  promoters <- promoters(txdb) %>% reduce

  gr <- append_feature(gr, introns, "intron")
  gr <- append_feature(gr, promoters, "promoter")
  return(gr)
}

##' Append the feature from gr_score to gr as a score along the interval (considering the strand direction)
##'
##' Example: gr is a granges containing the UTR3 regions. gr_score is the conservation score GRanges. We would like to add a column to gr representing the conservation score along the interval (like a view).
##'
##' @param gr GRanges with the region of interest
##' @param gr_score GRanges with a conservation score called `get(score_name)`
##' @param score_name Column name in gr_score representing the score of interest.
##' @param new_colname Name of the appended column to gr
##' @param new_colclass Class of the newly created column (Rle by default)
##' @return GRanges same to gr, but with added column `get(new_colname)`. NA means that the score wasn't present in gr_score.
append_score_view <- function(gr, gr_score, score_name = "score", new_colname = "conservation", new_colclass = NULL) {
  gr <- copy(gr)
  gr_score <- copy(gr_score)

  ## use only common seqlevels
  common_seqlevels <- intersect(seqlevels(gr),seqlevels(gr_score))
  gr <- gr[seqnames(gr) %in% common_seqlevels]
  gr_score <- gr_score[seqnames(gr_score) %in% common_seqlevels]
  seqlevels(gr) <- common_seqlevels
  seqlevels(gr_score) <- common_seqlevels
  ## compute the coverage
  gr_score_coverage <- coverage(gr_score)

  ## coverage can be only 0 or 1
  stopifnot(all((gr_score_coverage %>% lapply(unique) %>% unlist %>% unique %>% sort) %in% c(0,1)))
  
  ## gr_score_coverage %>% lapply(table)  ## sanity check - there only have to be numbers like 0 and 1
  gr_score_score_coverage<- coverage(gr_score, weight=score_name)

  ## add NA's where the coverage is 0
  for (i in 1:length(gr_score_score_coverage)) {
    gr_score_score_coverage[[i]][gr_score_coverage[[i]] == 0] <- NA
  }
  ## gr_score_score_coverage[gr]
  ## gr_score_score_coverage[[1]] %>% range(na.rm= TRUE)

  ## we need the same type as before
  ## gr$order <- 1:length(gr)
  ## gr$conservation <- RleList(lapply(1:length(gr), function(x) NA))
  ## gr_list <- split(gr, seqnames(gr))
  ## v <- Views(gr_score_score_coverage, as(gr_list, "RangesList"))
  ## for (i in 1:length(v)) {
  ##   rle_list <- RleList(v[[i]])
  ##   ## fix the strandedness
  ##   neg_strand <- strand(gr_list[[i]]) == "-"
  ##   rle_list[neg_strand] <- revElements(rle_list[neg_strand])
  
  ##   gr_list[[i]]$conservation <- rle_list
  ## }
  ## names(gr_list) <- NULL
  ## gr <- unlist(gr_list)
  ## gr <- gr[order(gr$order)]
  ## gr$order <- NULL

  ## stopifnot(all(table(width(gr) == elementLengths(mcols(gr)[[new_colname]]))))
  ## return(gr)
  cvg <- gr_score_score_coverage
  cvg.gr <- cvg[as(gr, "RangesList")]
  rle <- unsplit(cvg.gr, rep(seqnames(gr), width(gr)))

  ## rle <- unlist(cvg.gr)
  rle.range <- relist(rle, ranges(gr))
  rle.range[strand(gr)=="-"] <- revElements(rle.range[strand(gr)=="-"])

  mcols(gr)[[new_colname]] <- rle.range

  ## lengths have to match
  stopifnot(all(table(width(gr) == elementNROWS(mcols(gr)[[new_colname]]))))

  if (!is.null(new_colclass)) {
    mcols(gr)[[new_colname]] <- lapply(mcols(gr)[[new_colname]], as, Class = new_colclass)
  }
  
  return(gr)
}


##' Extend (pad) the granges range some bases upstream and some downstream wrt the strand
##' 
##' @param n_upstream Number of bases to extend upstream. 0 is the neutral element
##' @param n_upstream Number of bases to extend downstream. 0 is the neutral element
##' @return granges with extended start and stop with respect to the strand
extend_range <- function(gr, n_upstream = 0, n_downstream = 0) {
  gr <- copy(gr)
  gr <- resize(gr, width(gr) + n_downstream, fix = "start")
  gr <- resize(gr, width(gr) + n_upstream, fix = "end")
  return(gr)
}


##' reverseComplement for character vector
reverseComplement_char <- function(seq) {
  as.character(reverseComplement(DNAStringSet(seq)))
}

##' Create views for DNAStringSet using GRanges
##' 
##' @param seq DNAStringSet
##' @param gr Granges
##' @param reverse_seq Should we reverse_complement the sequence on the minus strand - see tests/test_granges.R for examples
##' @param colname Column name in the returned granges representing the sequence view
##' @return GRanges with a new column `get(colname)`
##' In case of error:
##' Error in fromXStringViewsToStringSet(x, out.of.limits = out.of.limits,  (from granges.R#223) : 
##' x' has "out of limits" views
##' You have to wrap `Views(seq[[chr]],ranges(gr_list[[chr]]))` with `DNAStringSet`
easy_Views <- function(seq, gr, reverse_seq = TRUE, colname = "seq") {
  gr_list <- split(gr, seqnames(gr))

  gr_list <- lapply(intersect(names(seq),names(gr_list)) , function(chr) {
    gr_tmp <- gr_list[[chr]]
    v <- Views(seq[[chr]],ranges(gr_list[[chr]]))

    mcols(gr_tmp)[[colname]] <- tryCatch(as.character(v),
                                         error = function(e) {
                                           warning("Out-of-view sequences exist. Trimming them.")
                                           as.character(trim(v)) #out of view sequences are trimmed
                                         }
                                         )
    return(gr_tmp)
  })

  gr <- Reduce(c, gr_list)
  
  if (isTRUE(reverse_seq)) {
    mcols(gr[strand(gr) == "-"])[[colname]] <- reverseComplement_char(mcols(gr[strand(gr) == "-"])[[colname]])
  }
  
  return(gr)
}

##' Transform the granges with a variable of type "list of atomic vectors" into a long granges, with one row per element (extend the list of atomic vectors)
##'
##' It is an equivalent to dt_expand_vector_list. Strand information is ignored.
##' @param gr GRanges
##' @param variables variable names from dt to be expanded (variables should be of class list of atomic vectors)
##' @param id_name how to name the unique_col_id
## ##' @param reverse_complement Should we reverse back the sequence on the negative strand?
##' @param position_name how to name the position column
##' @param Class Character vector of the corresponding classes. Can be also of length one.
##' @param Separator Character vector of separators. Can be also of length one.
##' @param reverse Should we reverse the variables on the negative strand
##' @param reverse_complement_character Should we reverse and complement the variables on the negative strand that
##' are of class character. 
##' @return Long GRanges with one additional column: "id_name" and changed `get(variable)` columns
##'
##' @details
##' - In case you get NA's, double check your Class attribute.
##'
##' TODO the function is somehow slow with reverse = TRUE
gr_expand_vector_list <- function(gr, variables, id_name = "row_id",
                                  Class = "character", separator = "",
                                  reverse = FALSE,
                                  reverse_complement_character = reverse
                                  ) { 
  position_name <- "position_column"
  if (any(c(id_name, position_name) %in% names(mcols(gr)))) {
    stop("You can't use the columnname: '", position_name, "' or the given id_name in the granges")
  }
  dt <- gr %>% as.data.frame %>% as.data.table

  dtl <- dt_expand_vector_list(dt, variables = variables,
                               id_name = "unique_col_id",
                               position_name = position_name,
                               Class = Class, separator = separator)
  dtl[, start := start + get(position_name) - 1L, by = unique_col_id]

  ## DONE - sanity check - sequence has to have the same length as granges
  if (!all(dtl[, .(start[.N], end[.N]), by = unique_col_id][, V1 == V2])) {
    stop("Sequence length doesn't match with GRanges. Did you use the right separator?")
  }

  dtl[, end := start]
  dtl[, (position_name) := NULL]

  ## reverse the sequence if necessary
  if (isTRUE(reverse)) {
    for (variable in variables) {
      ## for the character vector, reverse complement the sequence
      if (dtl[, class(get(variable))] == "character" & isTRUE(reverse_complement_character)) {
        dtl[strand == "-", (variable) := strsplit(reverseComplement_char(paste(get(variable), collapse = "")), split = "")[[1]],
            by = unique_col_id]
      } else {
        dtl[strand == "-", (variable) := rev(get(variable)), by = unique_col_id]
      }
    }
  }

  setnames(dtl, "unique_col_id", id_name)
  return(as(dtl, "GRanges"))
}

##' Inverse operation of gr_expand_vector_list
gr_shrink_vector_list <- function(gr, variables, id_name,
                                  reverse = TRUE,
                                  reverse_complement_character = reverse,
                                  separator = "") {

  dt <- gr %>% as.data.frame %>% as.data.table

  ## test for uniqueness
  dt[, start := start[1], by = get(id_name)]
  dt[, end := end[.N], by = get(id_name)]

  if (!all(dt[, -c(variables), with = F][, nrow(unique(.SD)), by = get(id_name)]$V1 == 1)) {
    stop("For each ", id_name, ", removing variables, the data table is not unique")
  }

  if (isTRUE(reverse_complement_character) & separator != "") {
    stop("For using the reverse complement character, one has to use the '' separator")
  }

  for (variable in variables) {
    if (is.character(dt[[variable]]) & isTRUE(reverse_complement_character)) {
      dt[, (variable) := if (strand == "-") {
        reverseComplement_char(paste(get(variable), collapse = separator))
      } else {
        paste(get(variable), collapse = separator)
      }, by = .(get(id_name), strand)]
    } else if (isTRUE(reverse)) {
      dt[, (variable) := paste(rev(get(variable)), collapse = separator), by = get(id_name)]
    } else {
      dt[, (variable) := paste(get(variable), collapse = separator), by = get(id_name)]
    }
  }
  
  dt <- dt[, unique(.SD), by = get(id_name)]
  stopifnot(all(dt[, .N, by = get(id_name)][,N] == 1))
  dt[, (id_name) := NULL]
  return(as(dt, "GRanges"))
}


##' Mutate (overwrite) the sequence of a View-like GRanges
##'
##' @param gr_ref Reference GRanges containing the sequence ranges (ala View) and the sequence.
##' This sequence will get mutated on places specified by gr_mutation.
##' @param ref_seq_colname Which column in gr_reference is representing the sequence
##' @param ref_reverse_complement Should we reverse-complement the sequence in gr_ref?
##' @param gr_mut Mutation GRanges containing the mutation sequence ranges (ala View) and the sequence to be replaced
##' @param mut_seq_colname Which column in gr_mutate is representing the sequence
##' @param mut_reverse_complement Should we reverse-complement the sequence in gr_mut?
##' 
##'
##' @return Modified gr_reference, with forced sequence as specified in gr_mutation
mutate_granges <- function(gr_ref,
                           ref_seq_colname = "seq",
                           ref_reverse_complement = TRUE,
                           gr_mut, 
                           mut_seq_colname = "seq",
                           mut_reverse_complement = FALSE
                           ) {
  ## specify protected column names
  grl_ref <- gr_expand_vector_list(gr_ref, variables = ref_seq_colname, Class = "character", separator = "",
                                   reverse = ref_reverse_complement,
                                   reverse_complement_character = ref_reverse_complement
                                   )
  
  grl_mut <- gr_expand_vector_list(gr_mut, variables = mut_seq_colname, Class = "character", separator = "",
                                   reverse = mut_reverse_complement,
                                   reverse_complement_character = mut_reverse_complement
                                   )

  ol <- findOverlaps(grl_ref, grl_mut)

  ## replace the corresponding sequences
  mcols(grl_ref[queryHits(ol)])[[ref_seq_colname]] <- mcols(grl_mut[subjectHits(ol)])[[mut_seq_colname]]

  ## expand gr_ref back
  gr_ref_ret <- gr_shrink_vector_list(grl_ref, ref_seq_colname, id_name = "row_id", reverse = ref_reverse_complement, separator = "")
  return(gr_ref_ret)
}

##' Get the codon context with
##' 
##' @param anno_gr annotation granges (from gtf)
##' @param fasta_seq DNAString i.e. the fasta file
##' @param start_or_stop options are "start" or "stop". Should we use the start or stop codon
##' @param bases_upstream,bases_downstream  by how many bases upstream/downstream should we extend the start/stop codon's
##' @return View with the corresponding sequence
get_codon_context <- function(anno_gr, fasta_seq, start_or_stop = "start", bases_upstream = 5, bases_downstream = 5) {
  stopifnot(start_or_stop %in% c("start", "stop"))
  anno_gr <- anno_gr[anno_gr$type == paste0(start_or_stop,"_codon")]
  anno_gr <- anno_gr[width(anno_gr) == 3]

  anno_gr <- extend_range(anno_gr, n_upstream = bases_upstream,
                          n_downstream = bases_downstream
                          )

  anno_gr <- easy_Views(fasta_seq, anno_gr, reverse_seq = TRUE)
  return(anno_gr)
}


##' Test if the overlaps are unique
is_unique_overlaps <- function(ol) {
  c(queryHits = is_unique(queryHits(ol)),
    subjectHits = is_unique(subjectHits(ol)))
}

##' Generic to duplicated_values but for overlaps
##' returns the rows that are duplicated
duplicated_values_overlaps <- function(ol) {
  ol <- ol[duplicated_values_bool(queryHits(ol)) | duplicated_values_bool(subjectHits(ol))]
  return(ol)
}


##' Wrapper for makeGRangesFromDataFrame
##'
##' This version uses as a default keep.extra.columns = TRUE.
##' Also, all the names specified for defining the name (non-standard names) are additionally appended to the end.
##' That way, when going back to the data.table, we have the same column names.
##' 
##' @param dt data.table or data.frame
##' @param ...
##' 
##' Use case - remove the code excess
##' 
##' gra <- makeGRangesFromDataFrame(dta, keep.extra.columns = TRUE, start.field = "start_gene_full", end.field = "end_gene_full")
##' gra$start_gene_full <- start(gra)
##' gra$end_gene_full <- end(gra)
##'
##' ## is replaced by:
##' gra2 <- dt2granges(dta, start.field = "start_gene_full", end.field = "end_gene_full")
dt2granges <- function(dt, ...) {
  params <- list(...)
  ## return(print(params))
  ## these have to be braught back:
  specified_params <- params[unlist(params) %in% names(dt)]

  dt <- as.data.frame(dt)
  gr <- makeGRangesFromDataFrame(df = dt, keep.extra.columns = TRUE, ...)
  
  add_cols <- unique(unlist(specified_params))

  if (length(add_cols) > 0){
    mcols(gr) <- DataFrame(mcols(gr),                                #old ones
                           dt[, add_cols, drop = FALSE],
                           check.names = FALSE) #extra ones
  }

  ## restore back the important columns
  return(gr)
}


##' Get the granges intersection centers
##'
##' @param gr_list List of GRanges
##' @return GRanges with width = 1
##' @author Žiga Avsec
intersect_center <- function(grlist) {
  if (length(grlist) == 0) return(NULL)
  gr_intersect <- Reduce(GenomicRanges::intersect, grlist) #subsetByOverlaps

  ## take peak middle-point
  gr_intersect <- resize(gr_intersect, width = 1, fix = "center")
  return(gr_intersect)  
}


##' Concatenate a list of genomic ranges together 
##' 
##' @param gr_list - list of GRanges
##' @param idcol name of the new appended column.
##' analogous to idcol in data.table::rbindlist
##' @return GRanges
gr_rbindlist <- function(gr_list, idcol="id") {
  ## all names have to be undefined
  if (length(gr_list) == 0) return(NULL)

  gr_list <- omit_null(gr_list)
  stopifnot(all(sapply(gr_list, function(x) is.null(names(x)))))
  
  grnew <- unlist(GRangesList(gr_list))
  mcols(grnew)[[idcol]] <- names(grnew)
  names(grnew) <- NULL

  stopifnot(sum(sapply(gr_list, length)) == length(grnew))
  return(grnew)
}
