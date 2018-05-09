# BST 281 QMP Project
- This is a final project for the Spring 2018 course in Genomic Data Manipulation (BST 281) at [Harvard T.H. Chan School of Public Health](https://www.hsph.harvard.edu/) on applying the [Quantitative Microbiome Profiling (QMP)](https://github.com/raeslab/QMP) method to longitudinal fecal samples from a participant with Crohn's Disease (CD) vs. a healthy control from the [NIH Human Microbiome Project (HMP2)](https://www.ibdmdb.org/).
- Group members are [Eunice Yeh](https://github.com/euniceyeh/), Marina Cheng, Anthony Lamattina, and Tian Zhang.
- Most raw data files are publicly available online:
  + [HMP2 Metadata](https://ibdmdb.org/tunnel/cb/document/Public/HMP2/Metadata/hmp2_metadata.csv)
  + [Zipped Metagenome Taxonomic Profiles (RMP)](https://ibdmdb.org/tunnel/cb/document/Public/HMP2/WGS/1812/taxonomic_profiles.tsv.gz)
  + [Select zipped HMP2 metagenome FASTQ raw files](https://ibdmdb.org/tunnel/public/HMP2/WGS/1812/rawfiles)
  + [Select zipped HMP2 Pilot metagenome FASTQ raw files](https://ibdmdb.org/tunnel/public/HMP2_Pilot/WGS/1644/rawfiles)
  + [Zipped Bacterial 16S rRNA gene copy number data files](https://rrndb.umms.med.umich.edu/static/download/rrnDB-5.4_pantaxa_stats_RDP.tsv.zip)
- Raw qPCR 16s DNA quantified in concentration data files from Casey DuLong, Tiffany Poon, and Jason Lloyd-Price are saved in `/data`.
- Analysis-ready datasets are saved within the subdirectory `/data/derived`.
- All relevant data visualization and statistical analysis results can be found under `/outputs`.
- Any python, R, and mathematica program file/script used to derive data or produce outputs by any of the team member is compiled into `/programs`.
- The main objective of this project is to compare our own variant of the novel QMP method proposed by [Vandeputte et al.](https://www.nature.com/articles/nature24460) to a standard relative profiling method widely used in the field by applying these methods on data collected and made publicly available by the HMP2 project, and from additional qPCR performed courtesy of [Curtis Huttenhower](https://huttenhower.sph.harvard.edu/) & colleagues named above from the [Broad Institute](https://www.broadinstitute.org/infectious-disease-microbiome).

Hope you enjoy our project!
