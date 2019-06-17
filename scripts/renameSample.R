renameSample <- function(seq_name) {
    require(xlsx)
    df <- read.xlsx(file.path("/fh/fast/tapscott_s/CompBio/RNA-Seq/hg38.FSHD.biopsyI-IV",
                              "extdata", "WellstoneProj2.Data.Alt.UR.ID.xlsx"),
                    sheetIndex=1, startRow=2, endRow=38, colIndex=1:2,
                    stringsAsFactors=FALSE)
    df[is.na(df$Alt.UR.ID), 1] <- df$record_id[is.na(df$Alt.UR.ID)]
    sample_name <- df$record_id
    sample_name <- sub("--", "-", sample_name)
    sample_name <- sub("01-", "", sample_name)
    alt_id <- df$Alt.UR.ID
    names(alt_id) <- sample_name
    #' add nine control samples
    cntr <- paste0("01-", "00", seq(41,49))
    names(cntr) <- c("2524", "2536", "2538", "2539", "2540",
                     "2550", "2559", "2563", "2578")
    #' 0952-1
    add_on <- "01-0022-1"
    names(add_on) <- "0952-1"
    #' combine
    alt_id <- c(alt_id, cntr, add_on)
    #' santity check
    if (all(seq_name %in% names(alt_id)))
        return(alt_id[seq_name])
    else {
        i <- which(!seq_name %in% names(alt_id))
        message("Cannot find alternative name to match ", paste(seq_name[i], collapse=","))
    }
}
