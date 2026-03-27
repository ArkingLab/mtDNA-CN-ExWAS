#######################################################################
# Genetic Relatedness
#######################################################################

#######################################################################
# Genetic relationship matrix (GRM)
#

snpgdsGRM <- function(gdsobj, sample.id=NULL, snp.id=NULL,
    autosome.only=TRUE, remove.monosnp=TRUE, maf=NaN, missing.rate=NaN,
    method=c("GCTA", "Eigenstrat", "EIGMIX", "Weighted", "Corr", "IndivBeta"),
    num.thread=1L, useMatrix=FALSE, out.fn=NULL, out.prec=c("double", "single"),
    out.compress="LZMA_RA", with.id=TRUE, verbose=TRUE)
{
    # check and initialize ...
    method <- match.arg(method)
    mtxt <- method
    if (method == "Weighted")
    {
        method <- "EIGMIX"
        mtxt <- "Weighted GCTA"
    } else if (method == "Corr")
    {
        mtxt <- "Scaled GCTA (correlation)"
    }

    stopifnot(is.logical(useMatrix), length(useMatrix)==1L)
    stopifnot(is.logical(with.id), length(with.id)==1L)
    ws <- .InitFile2(
        cmd=paste("Genetic Relationship Matrix (GRM, ", mtxt, "):", sep=""),
        gdsobj=gdsobj, sample.id=sample.id, snp.id=snp.id,
        autosome.only=autosome.only, remove.monosnp=remove.monosnp,
        maf=maf, missing.rate=missing.rate, num.thread=num.thread,
        verbose=verbose)
    cat('---HERE---')   
    if (!is.null(out.fn))
    {
        # gds output
        stopifnot(is.character(out.fn), length(out.fn)==1L)
        out.prec <- match.arg(out.prec)
        if (out.prec=="single") out.prec <- "float32"
        # create a gds file
        out.gds <- createfn.gds(out.fn)
        on.exit(closefn.gds(out.gds))
        put.attr.gdsn(out.gds$root, "FileFormat", "SNPRELATE_OUTPUT")
        put.attr.gdsn(out.gds$root, "version",
            paste0("SNPRelate_", packageVersion("SNPRelate")))
        add.gdsn(out.gds, "command", c("snpgdsGRM", paste(":method =", method)))
        add.gdsn(out.gds, "sample.id", ws$sample.id, compress=out.compress,
            closezip=TRUE)
        add.gdsn(out.gds, "snp.id", ws$snp.id, compress=out.compress,
            closezip=TRUE)
        sync.gds(out.gds)
        add.gdsn(out.gds, "grm", storage=out.prec, valdim=c(ws$n.samp, 0L),
            compress=out.compress)
    } else
        out.gds <- NULL

    # call GRM C function
cat('---HERE---')   
 rv <- .Call(gnrGRM, ws$num.thread, method, out.gds, useMatrix, verbose)

    # return
    if (is.null(out.gds))
    {
        if (isTRUE(useMatrix))
            rv <- .newmat(ws$n.samp, rv)
		if (with.id)
		{
			rv <- list(sample.id=ws$sample.id, snp.id=ws$snp.id,
				method=method, grm=rv)
            if (method %in% c("IndivBeta"))
                rv$avg_val <- .Call(gnrGRM_avg_val)
		}
        rv
    } else {
        if (method %in% c("IndivBeta"))
            add.gdsn(out.gds, "avg_val", .Call(gnrGRM_avg_val))
        invisible()
    }
}
