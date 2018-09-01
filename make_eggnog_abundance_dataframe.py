#!/usr/bin/env python
"""Join together a set of results based on their eggNOG annotations."""

import os
import io
import sys
import uuid
import copy
import time
import gzip
import json
import boto3
import shutil
import logging
import argparse
import traceback
import numpy as np
import pandas as pd
from collections import defaultdict


def exit_and_clean_up(temp_folder):
    """Log the error messages and delete the temporary folder."""
    # Capture the traceback
    logging.info("There was an unexpected failure")
    exc_type, exc_value, exc_traceback = sys.exc_info()
    for line in traceback.format_tb(exc_traceback):
        logging.info(line)

    # Delete any files that were created for this sample
    logging.info("Removing temporary folder: " + temp_folder)
    shutil.rmtree(temp_folder)

    # Exit
    logging.info("Exit type: {}".format(exc_type))
    logging.info("Exit code: {}".format(exc_value))
    sys.exit(exc_value)


def read_json(fp):
    assert fp.endswith((".json", ".json.gz"))
    logging.info("Reading in " + fp)
    if fp.startswith("s3://"):
        # Parse the S3 bucket and key
        bucket_name, key_name = fp[5:].split("/", 1)

        # Connect to the S3 boto3 client
        s3 = boto3.client('s3')

        # Download the object
        retr = s3.get_object(Bucket=bucket_name, Key=key_name)

        if fp.endswith(".gz"):
            # Parse GZIP
            bytestream = io.BytesIO(retr['Body'].read())
            got_text = gzip.GzipFile(
                None, 'rb', fileobj=bytestream).read().decode('utf-8')
        else:
            # Read text
            got_text = retr['Body'].read().decode('utf-8')

        # Parse the JSON
        dat = json.loads(got_text)

    else:
        assert os.path.exists(fp)

        if fp.endswith(".gz"):
            dat = json.load(gzip.open(fp, "rt"))
        else:
            dat = json.load(open(fp, "rt"))

    # Make sure that the sample sheet is a dictionary
    assert isinstance(dat, dict)

    return dat


def parse_gzipped_tsv(fp):
    assert fp.endswith(".tsv.gz")
    logging.info("Reading in " + fp)
    if fp.startswith("s3://"):
        # Parse the S3 bucket and key
        bucket_name, key_name = fp[5:].split("/", 1)

        # Connect to the S3 boto3 client
        s3 = boto3.client('s3')

        # Download the object
        retr = s3.get_object(Bucket=bucket_name, Key=key_name)

        # Parse GZIP
        bytestream = io.BytesIO(retr['Body'].read())
        got_text = gzip.GzipFile(
            None, 'rb', fileobj=bytestream).read().decode('utf-8')

        dat = got_text.split("\n")

    else:
        assert os.path.exists(fp)

        dat = gzip.open(fp, "rt").readlines()

    for line in dat:
        if len(line) == 0 or line[0] == '#':
            continue
        yield line.split("\t")


def read_eggnog_proportion_df(eggnog_annot, sample_sheet, results_key, abundance_key, gene_id_key):
    """Make a single DataFrame with the abundance (depth) from all samples for each eggNOG annotation."""

    # Collect all of the abundance information in this single dict
    dat = {}

    # Iterate over each sample
    for sample_name, sample_path in sample_sheet.items():
        # Get the JSON for this particular sample
        sample_dat = read_json(sample_path)

        # Make sure that the key for the results is in this file
        assert results_key in sample_dat

        # Subset down to the list of results
        sample_dat = sample_dat[results_key]
        assert isinstance(sample_dat, list)

        # Make sure that every element in the list has the indicated keys
        for d in sample_dat:
            assert abundance_key in d
            assert gene_id_key in d

        # Format as a Series
        depth = pd.Series({
            d[gene_id_key]: d[abundance_key]
            for d in sample_dat
        })

        # Sum up the depths by eggNOG annotations
        eggnog_depth = pd.Series({
            eggnog_id: np.sum([depth.get(gene_id, 0) for gene_id in list(gene_id_list)])
            for eggnog_id, gene_id_list in eggnog_annot.items()
        })

        # Remove the eggNOG annots with zero depth
        eggnog_depth = eggnog_depth.loc[eggnog_depth > 0]
        logging.info("Read in {:,} eggNOG annotations for {}".format(
            eggnog_depth.shape[0],
            sample_name
        ))

        # Save the proportion of the total sample assigned to each annotation
        dat[str(sample_name)] = eggnog_depth / depth.sum()

    logging.info("Formatting as a DataFrame")
    dat = pd.DataFrame(dat).fillna(0)

    logging.info("Read in data for {:,} eggNOG annotations across {:,} samples".format(
        dat.shape[0],
        dat.shape[1]
    ))

    return dat


def return_results(df, log_fp, output_prefix, output_folder, temp_folder):
    """Write out all of the results to a file and copy to a final output directory"""

    # Make sure the output folder ends with a '/'
    if output_folder.endswith("/") is False:
        output_folder = output_folder + "/"

    if output_folder.startswith("s3://"):
        s3 = boto3.resource('s3')

    for suffix, obj in [
        (".feather", df),
        (".logs.txt", log_fp)
    ]:
        if obj is None:
            "Skipping {}{}, no data available".format(output_prefix, suffix)
            continue

        fp = os.path.join(temp_folder, output_prefix + suffix)
        if suffix.endswith(".feather"):
            obj.reset_index().to_feather(fp)
        elif suffix.endswith(".json.gz"):
            json.dump(obj, gzip.open(fp, "wt"))
        elif suffix.endswith(".txt"):
            with open(fp, "wt") as f:
                f.write(obj)
        else:
            raise Exception(
                "Object cannot be written, no method for " + suffix)

        if output_folder.startswith("s3://"):
            bucket, prefix = output_folder[5:].split("/", 1)

            # Make the full name of the destination key
            file_prefix = prefix + output_prefix + suffix

            # Copy the file
            logging.info("Copying {} to {}/{}".format(
                fp,
                bucket,
                file_prefix
            ))
            s3.Bucket(bucket).upload_file(fp, file_prefix)

        else:
            # Copy as a local file
            logging.info("Copying {} to {}".format(
                fp, output_folder
            ))
            shutil.copy(fp, output_folder)


def read_eggnog_annot(eggnog_tsv_fp, eggnog_annot_field):
    """Read in the eggNOG TSV and group queries by annotation."""
    assert eggnog_annot_field in ["eggNOG", "KO", "GO"], eggnog_annot_field

    # Keys are the annotation, values are sets of gene IDs (queries)
    eggnog_annot = defaultdict(set)

    # Parse the eggNOG output
    for line in parse_gzipped_tsv(eggnog_tsv_fp):
        
        # Name for the query
        try:
            gene_id = line[0]
        except:
            logging.info("Problem")
            logging.info(line)
            assert 1 == 0

        # Parse the eggNOG annotation field
        if eggnog_annot_field == "KO":
            annots = line[6].split(",")
        elif eggnog_annot_field == "GO":
            annots = line[5].split(",")
        else:
            assert eggnog_annot_field == "eggNOG"
            annots = [line[1]]

        # Add this gene to the set of annotations
        for a in annots:
            if len(a) > 0:
                eggnog_annot[a].add(gene_id)

    return eggnog_annot


def calculate_proportions_by_eggnog_annot(df, eggnog_annot):
    """Calculate the proportion of each sample contained within each of the annotations."""

    # Get the intersection of genes with annotations that also are in the abundance DataFrame
    df_gene_id_set = set(df.index.values)
    eggnog_annot = {
        eggnog_id: gene_list & df_gene_id_set
        for eggnog_id, gene_list in eggnog_annot.items()
    }
    # Remove the annotations with 0 genes in this set of samples
    eggnog_annot_id_list = []
    eggnog_annot_gene_list_of_lists = []
    for annot_id, gene_list in eggnog_annot:
        if len(gene_list) > 0:
            eggnog_annot_id_list.append(annot_id)
            eggnog_annot_gene_list_of_lists.append(gene_list)

    logging.info("Grouping {:,} genes (out of a total set of {:,} genes) into {:,} eggNOG annotations".format(
        sum(map(len, eggnog_annot_gene_list_of_lists)),
        df.shape[0],
        len(eggnog_annot_gene_list_of_lists)
    ))

    eggnog_df = pd.DataFrame([
        df.loc[gene_list].sum()
        for gene_list in eggnog_annot_gene_list_of_lists
    ], index=eggnog_annot_id_list)

    # Calculate the proportion of each sample going into each group
    logging.info("Calculating proportional abundance per sample")
    eggnog_df = eggnog_df / df.sum()

    return eggnog_df


def make_eggnog_abundance_dataframe(
    eggnog_tsv_fp=None,
    eggnog_annot_field="KO",
    sample_sheet=None,
    output_prefix=None,
    output_folder=None,
    temp_folder="/scratch",
    results_key="results",
    abundance_key="depth",
    gene_id_key="id",
):
    # Make a new temp folder
    temp_folder = os.path.join(temp_folder, str(uuid.uuid4())[:8])
    os.mkdir(temp_folder)

    # Set up logging
    log_fp = os.path.join(temp_folder, "log.txt")
    logFormatter = logging.Formatter(
        '%(asctime)s %(levelname)-8s [eggNOG_abund_DF] %(message)s'
    )
    rootLogger = logging.getLogger()
    rootLogger.setLevel(logging.INFO)

    # Write to file
    fileHandler = logging.FileHandler(log_fp)
    fileHandler.setFormatter(logFormatter)
    rootLogger.addHandler(fileHandler)
    # Also write to STDOUT
    consoleHandler = logging.StreamHandler()
    consoleHandler.setFormatter(logFormatter)
    rootLogger.addHandler(consoleHandler)

    # READING IN DATA

    # Read in the eggNOG annotations
    logging.info("Reading in the eggNOG annotations from {}, grouping by {}".format(
        eggnog_tsv_fp,
        eggnog_annot_field
    ))
    try:
        eggnog_annot = read_eggnog_annot(
            eggnog_tsv_fp,
            eggnog_annot_field
        )
    except:
        exit_and_clean_up(temp_folder)

    # Read in the sample_sheet
    logging.info("Reading in the sample sheet from " + sample_sheet)
    try:
        sample_sheet = read_json(sample_sheet)
    except:
        exit_and_clean_up(temp_folder)

    # Make the abundance DataFrame
    logging.info("Making the abundance DataFrame")
    try:
        df = read_eggnog_proportion_df(
            eggnog_annot,
            sample_sheet,
            results_key,
            abundance_key,
            gene_id_key
        )
    except:
        exit_and_clean_up(temp_folder)

    # Read in the logs
    logs = "\n".join(open(log_fp, "rt").readlines())

    # Return the results
    logging.info("Returning results to " + output_folder)
    try:
        return_results(
            df,
            logs,
            output_prefix,
            output_folder,
            temp_folder
        )
    except:
        exit_and_clean_up(temp_folder)

    # Delete any files that were created for this sample
    logging.info("Removing temporary folder: " + temp_folder)
    shutil.rmtree(temp_folder)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="""Join together a set of results based on their eggNOG annotations"""
    )

    parser.add_argument("--eggnog-annot-field",
                        type=str,
                        default="KO",
                        help="""Annotation to normalize by, accepts eggNOG, KO, or GO.""")
    parser.add_argument("--eggnog-tsv-fp",
                        type=str,
                        required=True,
                        help="""Compressed TSV with eggNOG output.""")
    parser.add_argument("--sample-sheet",
                        type=str,
                        required=True,
                        help="""Location for sample sheet (.json[.gz]).""")
    parser.add_argument("--output-prefix",
                        type=str,
                        required=True,
                        help="""Prefix for output files.""")
    parser.add_argument("--output-folder",
                        type=str,
                        required=True,
                        help="""Folder to place results.
                                (Supported: s3://, or local path).""")
    parser.add_argument("--temp-folder",
                        type=str,
                        default="/scratch",
                        help="Folder for temporary files.")
    parser.add_argument("--results-key",
                        type=str,
                        default="results",
                        help="Key identifying the list of gene abundances for each sample JSON.")
    parser.add_argument("--abundance-key",
                        type=str,
                        default="depth",
                        help="Key identifying the abundance value for each element in the results list.")
    parser.add_argument("--gene-id-key",
                        type=str,
                        default="id",
                        help="Key identifying the gene ID for each element in the results list.")

    args = parser.parse_args(sys.argv[1:])

    # Sample sheet is in JSON format
    assert args.sample_sheet.endswith((".json", ".json.gz"))

    # Normalization factor is absent, 'median', or 'sum'
    assert args.eggnog_annot_field in ["eggNOG", "KO", "GO"]

    # Make sure the temporary folder exists
    assert os.path.exists(args.temp_folder), args.temp_folder

    make_eggnog_abundance_dataframe(
        **args.__dict__
    )
