# Entropic Ranks (docker) readme

entropic_ranks

Description: Performs an Entropic Ranks analysis on a data set, returning a list containing downregulated and upregulated features. May be used supervised, returning the full feature list and printing the suggested cutoff points for later manual trimming, or unsupervised, returning only the information-rich feature list. In the unsupervised mode, the lists of information-rich features may be exported as tab-delimited .txt files automatically.

## Usage:
Usage within an R script

```{r}
entropic_ranks(data_under_analysis,population_vector,data_origin=NULL,granularity=1,supervised=FALSE,process_log=FALSE,export_plots=FALSE,create_output_files=FALSE,is_logged=TRUE,logbase=2,huge_feature_list=FALSE)
```

## Arguments:
**data_under_analysis** - Tab-delimited .txt table with rows representing features, columns representing samples and cells containing the values to be compared. Rownames and column names must be unique. (see included test data)

**population_vector** - Tab-delimited .txt table with a single column, one row per sample and 0 or 1 as the table values. Header and rownames must be included. Denotes the two sample subpopulations to be compared. (see included test data)

**data_origin** - Tab-delimited .txt table with a single column, one row per sample. Header and rownames must be included. The data must be labels differentiating the origin of each sample. If NULL, it defaults to assuming that data are of the same origin. To be used only if the data are from different experiments or publications. (default: null)

**granularity** - The sliding window step, corresponding to the granularity of the partitioning process (feature-by-feature, or partitioning by 5-feature steps). (default: 1)

**supervised** - If TRUE, the full list of differentially behaving features is returned and the tables of suggested cutoff points are printed. If FALSE, only the list of information-rich features is returned. (default: FALSE)

**process_log** - If TRUE, statistics of the entropic_ranks execution will be printed and plots of the entropy distributions and clustering qualities will be generated. (default: FALSE)

**export_plots** - If TRUE, png plots of the entropy distributions and clustering qualities will be exported as files in a folder system created in the current working directory. (default: TRUE)

**create_output_files** - If TRUE, the feature lists of information-rich features will be automatically exported in the working directory as tab-delimited .txt files. Ignored if supervised is set to TRUE. (default: TRUE)

**is_logged** - Set to TRUE if the values are log-transformed and you want to export the Fold Change instead of the Log Fold Change in .txt files. Ignored if supervised is set to TRUE. (default: TRUE)

**logbase** - The base of the log transformation. Ignored if supervised is set to TRUE or if create_output_files is set to FALSE. (default: 2)

**huge_feature_list** - Only set to TRUE if the entropic_ranks fails to run due to huge feature lists returned by RankProd (e.g. more than 25000-30000 features) and you can be reasonably sure that less than 1000 are differentially expressed and information-rich. If TRUE, entropic_analysis will only investigate the first 20000 features and isolare information-rich features from among them. The issue may consistently appear due to RAM shortage when analyzing methylation data. (default: FALSE)


# Examples of docker usage:
(assuming the docker image is named "entropic_ranks")


Full-parameter usage (using default values as described above):

```console
docker run --rm -v "/your/data/here:/data entropic_ranks Rscript Entropic_Ranks.R /data/GSE_data_set.txt /data/vec.txt null 1 FALSE FALSE TRUE TRUE TRUE 2 FALSE
```

Full-parameter usage (using a data origin file, supervised without output files and plots):

```console
docker run --rm -v /your/data/here:/data entropic_ranks Rscript Entropic_Ranks.R /data/GSE_data_set.txt /data/vec.txt /data/data_origin_file.txt 1 TRUE FALSE FALSE FALSE TRUE 2 FALSE
```


Default usage (the two input files *must* be named "data_table.txt" and "population_vector.txt" and all samples have the same origin):

```console
docker run --rm -v /your/data/here:/data entropic_ranks
```
