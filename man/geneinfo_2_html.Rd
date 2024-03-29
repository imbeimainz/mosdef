% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/report_helpers.R
\name{geneinfo_2_html}
\alias{geneinfo_2_html}
\title{Information on a gene}
\usage{
geneinfo_2_html(gene_id, res_de = NULL, col_to_use = "SYMBOL")
}
\arguments{
\item{gene_id}{Character specifying the gene identifier for which to retrieve
information}

\item{res_de}{A \code{DESeqResults} object, storing the result of the differential
expression analysis. If not provided, the experiment-related information is not
shown, and only some generic info on the identifier is displayed.
The information about the gene is retrieved by matching on the \code{SYMBOL} column,
which should be provided in \code{res_de}.}

\item{col_to_use}{The column of your res_de object containing the gene symbols.
Default is "SYMBOL"}
}
\value{
HTML content related to a gene identifier, to be displayed in
web applications (or inserted in Rmd documents)
}
\description{
Assembles information, in HTML format, regarding a gene symbol identifier
}
\details{
Creates links to the NCBI and the GeneCards databases
}
\examples{
geneinfo_2_html("ACTB")
geneinfo_2_html("Pf4")
}
