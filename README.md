# docker-eggnog-mapper
Docker image for eggNOG mapper

[![Docker Repository on Quay](https://quay.io/repository/fhcrc-microbiome/eggnog-mapper/status "Docker Repository on Quay")](https://quay.io/repository/fhcrc-microbiome/eggnog-mapper)

This is a very thin wrapper around the eggNOG mapper tool.
See [https://github.com/jhcepas/eggnog-mapper](https://github.com/jhcepas/eggnog-mapper)
for more details on that tool.

The main addition is the inclusion of a `run_eggnog_mapper.py` script, which handles:

  * Downloading the reference database from an S3 bucket
  * Writing the results to a local path or S3 bucket

There is also a script that will parse the eggNOG output and generate a SQLite database 
with all of the reaction metadata for the detected set of KEGG orthologs fetched from the
KEGG API. The tables and columns are:

  * `ortholog`: [`ortholog`, `name`, `definition`]
  * `reaction`: [`reaction`, `ortholog`, `definition`, `equation`, `enzyme`]
  * `compound`: [`compound`, `name`, `formula`]
  * `pathway`: [`pathway`, `reaction`, `name`, `class`]

NOTE: The `reaction` table maps to all of the other tables, which in the case of the
`compound` table is via the `equation` column, which contains `compound` entry names.