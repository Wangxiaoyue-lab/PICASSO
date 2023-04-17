# Prepare loom files

Install and read <http://velocyto.org/velocyto.py/tutorial/cli.html#run10x-run-on-10x-chromium-samples>

**1. run `velocyto.sh` first:**

*Arguments*
- msk_gtf: To do so you would need to download an appropriate expressed repeat annotation (for example from [UCSC genome browser](https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=611454127_NtvlaW6xBSIRYJEBI0iRDEWisITa&clade=mammal&org=Mouse&db=mm10&hgta_group=allTracks&hgta_track=rmsk&hgta_table=0&hgta_regionType=genome&position=chr12%3A56694976-56714605&hgta_outputType=primaryTable&hgta_outputType=gff&hgta_outFileName=mm10_rmsk.gtf) and make sure to select GTF as output format).
- annotation: The `/public/home/luoliheng/refgenome/refdata-gex-mm10-2020-A/genes/genes.gtf` download from <https://support.10xgenomics.com/single-cell-gene-expression/software/release-notes/build>
- file_path: the output dir of 10x data


One hour for each sample.
**$!$ Make sure you use the correct species and version of annotation files.**
**$!$ Run the script after modifying the parameters. The current script is just an example.**

**2. run `combine.sh`:**

Ckeck if the `argparse` and `loompy` have been installed.

