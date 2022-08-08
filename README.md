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

[hpo_def](https://github.com/galerp/Cube3/blob/main/Files/HPO_def_rl_2020-10-12_dl_2021-08-03.csv) - This file contains every HPO code along with its definition

[hpo_ancestors](https://github.com/galerp/Cube3/blob/main/Files/HPO_ancs_rl_2020-10-12_dl_2021-08-03.csv) -  This file contains the higher level terms of each HPO term. This is essential for propagating HPO terms (e.g., creating the binned prop file from the binned base file.

[gene diagnoses](https://github.com/galerp/Cube3/blob/main/Files/example_gene_data.csv) -  This file contains the genetic diagnoses for those in the cohort as well as the age of genetic diagnosis.

[gene classes](https://github.com/galerp/Cube3/blob/main/Files/gene_classes.csv) -  This file has contains maps the genes into broader gene classes

[binned base](https://github.com/galerp/Cube3/blob/main/Files/example_bin_base.csv) -  This file contains the binned, base HPO terms of the cohor

[binned prop](https://github.com/galerp/Cube3/blob/main/Files/example_bin_prop.csv) -  This file contains the binned, propagated HPO terms of the cohor



### Requirements:
  [R](https://www.r-project.org/) with the following packages:
1. [optparse](https://cran.r-project.org/web/packages/optparse/index.html)
2. [yaml](https://cran.r-project.org/web/packages/yaml/index.html)
3. [tidyverse](https://cran.r-project.org/web/packages/tidyverse/index.html)
4. [hmisc](https://cran.r-project.org/web/packages/hmisc/index.html)
5. [foreach](https://cran.r-project.org/web/packages/foreach/index.html)
6. [doparallel](https://cran.r-project.org/web/packages/doparallel/index.html)

### Submitting :

```
Rscript config.R --input input.yml
```
