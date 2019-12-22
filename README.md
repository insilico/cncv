
# Consensus Nested Cross-Validation

#### Websites

[insilico Github Organization](https://github.com/insilico)

[insilico McKinney Lab](http://insilico.utulsa.edu/)

#### Related References. 

[2018 EpistasisRank and EpistasisKatz paper in Bioinformatics](https://doi.org/10.1093/bioinformatics/bty965)

### To install:

    >library(devtools)
    
    >install_github("insilico/cncv")  # todo (optional build_vignettes = TRUE)
    >library(cncv)
    >data(package="cncv")
    
    # >vignette(" ") # todo (if you build_vignettes)
    
### Dependencies
To install the `privateEC` collection of R packages:

```
install_github("insilico/privateEC")
```
To install the `PriorKnowledgeEpistasisRank` collection of R packages:

```
install_github("insilico/PriorKnowledgeEpistasisRank")
```

Other packages
```
install.packages(c('randomForest', 'glmnet', 'xgboost', 'caret', 'CORElearn))

```

### Examples

[https://github.com/insilico/cncv/blob/master/Vignettes/cncvExample.md](https://github.com/insilico/cncv/blob/master/Vignettes/cncvExample.md)



### Abstract

Motivation: Appropriate steps must be taken to avoid overfitting when using feature selection to improve the accuracy of machine learning models. Nested cross-validation (nCV) is a common approach that chooses the best classification model and best features to represent a given outer fold based on the maximum inner-fold accuracy. Differential privacy is a related technique to avoid overfitting that uses a privacy preserving noise mechanism to identify features that are stable between training and holdout sets. 
Methods: We develop consensus nested CV (cnCV) that combines the idea of feature stability from differential privacy with nested CV. Feature selection is applied in each inner fold and the consensus of top features across folds is a used as a measure of feature stability or reliability instead of classification used in standard nCV. We use simulated data with main effects, correla-tion, and interactions to compare the classification accuracy and feature selection performance of the new cnCV with standard nCV, Elastic Net optimized by CV, differential privacy, and pri-vate Evaporative Cooling (pEC). We also compare these methods using real RNA-Seq data from a study of major depressive disorder.
Results: The cnCV method has similar training and validation accuracy to nCV, but cnCV has much shorter run times because it does not construct classifiers in the inner folds. The cnCV method chooses a more parsimonious set of features with fewer false positives than nCV. The cnCV method has similar accuracy to pEC and cnCV selects stable features between folds without the need to specify a privacy threshold. We show that cnCV is an effective and efficient approach for combining feature selection with classification. 


#### Contact
[brett.mckinney@gmail.com](brett.mckinney@gmail.com)
