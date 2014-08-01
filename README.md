### What is this repository for? ###

* R client for accessing MyGene.info annotation and query services
* [Package vignette](https://bytebucket.org/sulab/mygene.r/raw/b5d3312762a3642129c2d01125f8c77b1e053cb2/mygene/inst/doc/mygene.pdf)

### To install ***mygene***: ###

* Download [***mygene***](https://bitbucket.org/sulab/mygene.r/downloads)


```
#!r

source("http://bioconductor.org/biocLite.R")
biocLite("IRanges")
```


```
#!r

install.packages(c("httr", "jsonlite", "Hmisc"))
```

```
#!r

install.packages("mygene_0.99.0.tar.gz", repos = NULL, type = "source")
```


### Contact ###

* help@mygene.info or adammark@scripps.edu