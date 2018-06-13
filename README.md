### What is this repository for? ###

* R client for accessing MyGene.info annotation and query services
* [Package vignette](http://bioconductor.org/packages/release/bioc/vignettes/mygene/inst/doc/mygene.pdf)

### To install ***mygene***: ###

* Install from Bioconductor by typing in your R console:

```
source("http://bioconductor.org/biocLite.R")
biocLite("mygene")
```

* Or download tarball from [Bioconductor](http://bioconductor.org/packages/release/bioc/html/mygene.html)

### Contact ###

* help@mygene.info

### How to push changes to bioconductor  ###
  
  Ref: https://bioconductor.org/developers/how-to/git/push-to-github-bioc/
  
  * Make sure you have the commit permission to git.bioconductor.org
  
    Setup at https://git.bioconductor.org/BiocCredentials/login/
    
  * Make sure proper git remote upstream is added following the instruction above
  
    ```
    $ git remote -v
    origin  git@github.com:biothings/mygene.R.git (fetch)
    origin  git@github.com:biothings/mygene.R.git (push)
    upstream        git@git.bioconductor.org:packages/mygene.git (fetch)
    upstream        git@git.bioconductor.org:packages/mygene.git (push)
    ```
  
  * Make changes and commit to origin/master
  
  * Merge and push changes to git.bioconductor.org
  
    ```
    git fetch upstream
    git merge upstream/master
    ```
  
    Fix any merge conflicts and push the changes to upstream
    
    ```
    git push merge upstream/master
    ```
