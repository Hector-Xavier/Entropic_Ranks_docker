FROM r-base:latest
USER root
RUN apt-get update && apt-get -y upgrade && apt-get install libcurl4-openssl-dev
RUN apt-get -y install libudunits2-dev libgdal-dev libgeos-dev libproj-dev && apt-get -y install libfontconfig1-dev && apt-get -y install libcairo2-dev
RUN apt-get -y install pandoc

WORKDIR /tmp

RUN wget https://bioconductor.org/packages/3.3/bioc/src/contrib/RankProd_2.44.0.tar.gz 
RUN wget https://cran.r-project.org/src/contrib/entropy_1.2.1.tar.gz
RUN R -e "install.packages(c('/tmp/RankProd_2.44.0.tar.gz','/tmp/entropy_1.2.1.tar.gz'), repos = NULL, type = 'source')"
RUN R -e "install.packages(c('factoextra','plotly'), dependencies = TRUE)"

RUN cd .. && rm -R tmp
RUN mkdir /data
RUN mkdir /home/Entropic_Ranks
ADD Entropic_Ranks.R /home/Entropic_Ranks/Entropic_Ranks.R

WORKDIR /home/Entropic_Ranks

CMD Rscript Entropic_Ranks.R /data/data_table.txt /data/population_vector.txt null 1 FALSE FALSE TRUE TRUE TRUE 2 FALSE
