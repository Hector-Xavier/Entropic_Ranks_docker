# Entropic Ranks (docker) Readme

# Description:
Performs an Entropic Ranks analysis on a data set, returning a list containing downregulated and upregulated features. May be used supervised, returning the full feature list and printing the suggested cutoff points for later manual trimming, or unsupervised, returning only the information-rich feature list. In the unsupervised mode, the lists of information-rich features may be exported as tab-delimited .txt files automatically.

# Usage within an R script

```console
entropic_ranks(data_under_analysis,population_vector,data_origin=NULL,granularity=1,supervised=FALSE,process_log=FALSE,export_plots=FALSE,create_output_files=FALSE,is_logged=TRUE,logbase=2,huge_feature_list=FALSE)
```

## Arguments in detail:
**data_under_analysis** - Tab-delimited .txt table with rows representing features, columns representing samples and cells containing the values to be compared. Rownames and column names must be unique. (see sample data included in "Test_data_GSE69486")

**population_vector** - Tab-delimited .txt table with a single column, one row per sample and 0 or 1 as the table values. Header and rownames must be included. Denotes the two sample subpopulations to be compared. (see sample data included in "Test_data_GSE69486")

**data_origin** - Tab-delimited .txt table with a single column, one row per sample. Header and rownames must be included. This file is similar in most respects with the "population_vector" .txt file, except that the values must be labels differentiating the origin of each sample (e.g., when performing a meta-analysis of data generated by different labs). If NULL, it defaults to the assumption that data are of the same origin. To be used only if the data are from different experiments or publications.
> default value: null)

**granularity** - The sliding window step, corresponding to the granularity of the partitioning process (e.g. feature-by-feature, or partitioning by 5-feature steps). May be used to speed up analysis and/or verify convergence.
> default value: 1

**supervised** - If TRUE, the full list of differentially behaving features is returned and the tables of suggested cutoff points are printed in the terminal. If FALSE, only the list of information-rich features is returned, the trimming having been performed automatically.
> default value: FALSE

**process_log** - If TRUE, statistics and details of the Entropic Ranks execution will be printed in the terminal and plots of the entropy distributions and clustering qualities will be generated (the latter would only be viewable when running the method in an R GUI, but not as a dockerized application). This would facilitate post-analysis manual assessment of data quality and result reliability, but would be unnecessary if the desired result was the integration in a fully automated pipeline.
> default value: FALSE

**export_plots** - If TRUE, plots of the entropy distributions and clustering qualities will be exported as .png files in a folder system created in the current working directory.
> default value: TRUE

**create_output_files** - If TRUE, the feature lists of information-rich features will be automatically exported in the working directory as tab-delimited .txt files. Ignored if **supervised** is set to TRUE.
> default value: TRUE

**is_logged** - Set to TRUE if the values are log-transformed and the user requires that the Fold Change be readily available (instead of the Log Fold Change). This setting is intended to facilitate integration in automated analysis pipelines and is ignored if **supervised** is set to TRUE.
> default value: TRUE

**logbase** - The base of the log transformation the data have been subjected to as part of their pre-processing and preparation. Ignored if **supervised** is set to TRUE or if **create_output_files** is set to FALSE.
> default value: 2

**huge_feature_list** - Only set to TRUE if the Entropic Ranks fails to run due to huge unsupportable RAM requirements caused by excessively lengthy feature lists (e.g. more than 30000 features) and the researcher is reasonably certain that less than 1000 should be considered differentially expressed and information-rich. If TRUE, entropic_analysis will only investigate the first 20000 features (sorted by rank product) and isolare information-rich features from among them. The issue may consistently come up during analysis intended to identify differentially methylated sites, for example. **This setting should be used with caution and proper mention, since it can introduce bias.**
> default value: FALSE


# Dockerized usage
In the following examples, we will be assuming the docker image is to be named named "entropic_ranks". Moreover, we will only clarify important variables and arguments that have not been referred to so far.

## Docker image creation
(requires "Docker" to be installed)
```console
git clone https://github.com/Hector-Xavier/Entropic_Ranks_docker
```
```console
docker build -t entropic_ranks Entropic_Ranks_docker
```
**-t entropic_ranks** - Sets the name of the resulting docker image to a name of the user's preference, so that repeated usage is facilitated. In our examples, "entropic_ranks" was used.


## Full-parameter usage example
Using default values.
```console
docker run --rm -v "/your/data/here:/data entropic_ranks Rscript Entropic_Ranks.R /data/GSE_data_set.txt /data/vec.txt null 1 FALSE FALSE TRUE TRUE TRUE 2 FALSE
```

**-v "/your/data/here:/data entropic_ranks** - The "-v" flag signifies that Docker will use a folder from the researcher's system as a revolving door between the Docker container's (temporary) existence and the (more permanent) file system of the hosting machine. The "/your/data/here" part of the argument is the path to a folder in the host system where the data to be subjected to analysis are located, and should be set by the researcher. The "/data" part of the argument is a path determined by the container's internal structute and should not be changed.

**entropic_ranks** - The (freely chosen) name of the docker image set during its creation.

**Rscript Entropic_Ranks.R** - This argument calls the Rscript performing the actual analysis and should not be changed by the researcher.

**/data/GSE_data_set.txt** - This argument denotes the .txt file containing the data table to be subjected to analysis by Entropic Ranks. The "/data/" part of the argument should not be changed. The relative path to the data table file should follow it, if it is not directly in the "/your/data/here" path mentioned above.

**/data/vec.txt** - This argument denotes the .txt file containing the population vector used to define the two data population that should be compared against each other. The "/data/" part of the argument should not be changed. The relative path to the population vector file should follow it, if it is not directly in the "/your/data/here" path mentioned above.

**/data/GSE_data_set.txt /data/vec.txt null 1 FALSE FALSE TRUE TRUE TRUE 2 FALSE** - The arguments correspond 1-1 to the arguments of Entropic Ranks.


## Customized usage example
Using a data origin file, supervised execution, without creating output files or plots.
```console
docker run --rm -v /your/data/here:/data entropic_ranks Rscript Entropic_Ranks.R /data/GSE_data_set.txt /data/vec.txt /data/data_origin_file.txt 1 TRUE FALSE FALSE FALSE TRUE 2 FALSE
```

## Fully default usage
The two input files *must* be named "data_table.txt" and "population_vector.txt", all samples *must* have the same origin and all other parameters *will be set to default*.
```console
docker run --rm -v /your/data/here:/data entropic_ranks
```

# Additional useful commands

## Export a docker image to a tar archive file
Useful for archiving purposes and a pre-requisite to transfer the docker image to a machine instead of re-building it (e.g. to work around a lack of internet connectivity).
```console
docker image save -o Entropic_Ranks_image_2020-10-22.tar.gz entropic_ranks
```

**-o Entropic_Ranks_image_2020-10-22.tar.gz  entropic_ranks** - The "-o" flag signifies that Docker should save the image to a file. The "Entropic_Ranks_image_2020-10-22.tar.gz" part of the argument is a freely-chosen filename for the archived docker image. The "entropic_ranks" part of the argument is the name for the docker image chosen at the time of creation.


## Load a docker image from a tar archive file
This allows a previously created tar archive file of a docker image to be loaded and prepared to be used to create new containers. Useful for archiving purposes and a pre-requisite to transfer the docker image to a machine instead of re-building it (e.g. to work around a lack of internet connectivity).
```console
docker image load -i Entropic_Ranks_image_2020-10-22.tar
```

**-i Entropic_Ranks_image_2020-10-22.tar** - The "-i" flag signifies that Docker will attempt to parse and load the contents of a tar archive file. The "Entropic_Ranks_image_2020-10-22.tar" part of the argument is the current name of the docker image tar archive file to be loaded (it can be a full path to the file).


## Getting access to the newly-created files
Traditionally, docker execution tends to create conflicts regarding user rights to files created by a docker and transferred to the host system as Entropic Ranks does with all its non-terminal output. This command will allow a user that should normally have the rights to move and delete files extend them to the newly-created result files.
```console
sudo chown -R $(id -u) /your/data/here
```
