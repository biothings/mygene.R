### What is this repository for? ###

* R client for accessing MyGene.info annotation and query services
* [Package vignette](https://bytebucket.org/sulab/mygene.r/raw/b5d3312762a3642129c2d01125f8c77b1e053cb2/mygene/inst/doc/mygene.pdf)

### To install ***mygene***: ###

* Install dependencies
```
#!r
source("http://bioconductor.org/biocLite.R")
biocLite(c("IRanges", "httr", "jsonlite", "Hmisc", "BiocStyle"))
```
* You can either install from package archive file
by downloading [***mygene***](https://bitbucket.org/sulab/mygene.r/downloads)
```
#!r

install.packages("mygene_0.99.0.tar.gz", repos = NULL, type = "source") 
```
* or install using devtools
```
#!r

library(devtools) 
install_bitbucket("mygene.R", "sulab", ref="default")
```

### Contact ###

* help@mygene.info or adammark@scripps.edu