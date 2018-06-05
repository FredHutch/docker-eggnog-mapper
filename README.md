# docker-eggnog-mapper
Docker image for eggNOG mapper

[![Docker Repository on Quay](https://quay.io/repository/fhcrc-microbiome/eggnog-mapper/status "Docker Repository on Quay")](https://quay.io/repository/fhcrc-microbiome/eggnog-mapper)

This is a very thin wrapper around the eggNOG mapper tool.
See [https://github.com/jhcepas/eggnog-mapper](https://github.com/jhcepas/eggnog-mapper)
for more details on that tool.

The main addition is the inclusion of a `run_eggnog_mapper.py` script, which handles:

  * Downloading the reference database from an S3 bucket
  * Writing the results to a local path or S3 bucket
