\name{LibraryReport}
\alias{LibraryReport}
\title{Mass Spectral Library Report}
\description{Generate report for the mass spectral library.}
\usage{LibraryReport(spectra = SpecDataEgg,
                          metadata = MetaDataEgg,
                          structures = paste(system.file(package = "SpecLibBlackSkimmerEgg2014"), "/extdata/struct/", sep = ""),
                          pdfFile = "SpecLibBlackSkimmerEgg2014.pdf",
                          pdfTitle = "SpecLibBlackSkimmerEgg2014 Library",
                          xMin = 40)
}

\arguments{
   \item{spectra}{data frame containing the mass spectra.}
    \item{metadata}{data frame containing the metadata.}
    \item{structures}{character string of directory containing the structure image files.}
    \item{pdfFile}{character string. The file name of the PDF report.}
    \item{pdfTitle}{character string. Title string to embed in the PDF file.} 
    \item{xMin}{numeric value. The minimum m/z in the plot.}
}
\details{See the OrgMassSpecR - Mass Spectral Libraries vignette for an
  overview and examples.
}
\value{A PDF report of the mass spectra.  If PDF = FALSE, the report is
  generated using the graphics device.
}

