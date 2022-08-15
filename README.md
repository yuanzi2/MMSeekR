# Title: MMSeekR
Segmentation of Bis-seq data using a multi-model sequence aware approach

#### Download the packages from github
*Maybe they are private repositories now, so can't use devtools::install_github to directly downlod them.*
  Use git clone to download them in local position. Open terminal or git bash or Rstudio (Tools-->Shell)

Obtain MMSeekR.data
```
git clone git@github.com:yuanzi2/MMSeekR.data.git
```

Obtain MMSeekR
```
git clone git@github.com:yuanzi2/MMSeekR.git
```

### Enter R environment
Install dependencies ‘MethylSeekR’, ‘valr’, ‘zoocat’
```
if (!requireNamespace("BiocManager", quietly = TRUE)){install.packages("BiocManager")}
BiocManager::install("MethylSeekR")
BiocManager::install("valr")
if(!requireNamespace("devtools", quietly = TRUE)) {install.packages("devtools")}
devtools::install_github("ran-ran/zoocat")  ##Which would you like to update? Please choose 3
```

Enter the path where the MMSeekR.data located
```
devtools::install("MMSeekR.data") ##Which would you like to update? Please choose 3
```

Enter the path where the MMSeekR located
```
devtools::install("MMSeekR")  ##Which would you like to update? Please choose 3
```

#### How to call PMDs
Check the runMultiModel help and you can test the examples.
```
library("MMSeekR")
?runMultiModel
```
