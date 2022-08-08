# Phenome-genome associations prior to genetic diagnosis

This repository is composed of a series of scripts which analyzes time-stamped clinical concepts extracted from electronic medical records (EMR) of a population of which a significant sub-population has a genetic diagnosis and the age at which the respective diagnosis was clinically made.
The clinical concepts are extracted from the free-text of all available EMR of the respective population using a customized form of the cTAKES pipeline. Only clinical concepts cTAKES has tagged as positive (non-negated), non-conditional (e.g., 'patient MAY have [clinical feature]'), and related to the subject (i.e., the patient and not a family member) are used. Clinical concepts are then translated to Human Phenotype Ontology (HPO) terms. Using the time-stamp from the original EMR from which a concept has been extracted, HPO terms are placed into 3 month time bins according to the age of the patient. All higher levels are inferred through 'propagation' and duplicate terms within each patient's time-bin are removed. All individuals' terms prior to their genetic diagnosis are removed.
Consequently, at each age bin we can find clinical concepts that are significantly associated with a particular genetic etiology, group of genetic etiologies, or any genetic prior to clinical genetic diagnosis.


## Scripts:

This [wrapper](https://github.com/shiva-g/The-Cube/blob/master/wrapper.R) script needs to be submitted to run the entire pipeline.

* [Helper file](https://github.com/shiva-g/The-Cube/blob/master/scripts/helper_file.R)  - loads data and cleans it. Creates base and prop hpo files. 
* [3d Array creation](https://github.com/shiva-g/The-Cube/blob/master/scripts/3d_arrays.R) - creates the 3d matrices.
* [Fishers test](https://github.com/shiva-g/The-Cube/blob/master/scripts/hpo_associations.R) - hpo_associations.


## Files: ##

[hpo_is.a_tree](https://github.com/shiva-g/The-Cube/blob/master/files/hpo_is.a_tree.csv) - This file contains the ontological information and definition for every single HPO term. The 'is.a' term is the parent term for each respective HPO term.

[hpo_ancestors](https://github.com/shiva-g/The-Cube/blob/master/files/hpo_ancestors.csv) -  This file contains the higher level terms of each HPO term. This is essential for calculating the MICA (most informative common ancestor) between 2 HPO terms which is one of the first steps to finding the similarity score between patients.

[aed_encounter](https://github.com/shiva-g/The-Cube/blob/master/files/aed_encounter.csv) -  This file contains the first and last encounters for each study ID.

[diagnosis](https://github.com/shiva-g/The-Cube/blob/master/files/diagnosis.csv) -  This file has all the diagnosis for each study ID at all time points along with the HPO terms.

[survival](https://github.com/shiva-g/The-Cube/blob/master/files/survival.csv) -  This file containts information about the genetic diagnosis of each study ID along with gender and diagnosis category.



### Requirements:
  [R](https://www.r-project.org/) with the following packages:
1. [optparse](https://cran.r-project.org/web/packages/optparse/index.html)
2. [yaml](https://cran.r-project.org/web/packages/yaml/index.html)
3. [tidyverse](https://cran.r-project.org/web/packages/tidyverse/index.html)
4. [dplyr](https://cran.r-project.org/web/packages/dplyr/index.html)
5. [ggplot2](https://cran.r-project.org/web/packages/ggplot2/index.html)
6. [ggrepel](https://cran.r-project.org/web/packages/ggrepel/index.html)

### Submitting :

```
Rscript wrapper.R --input input.yml
```
