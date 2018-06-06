#!/usr/bin/env python
"""Run the eggNOG mapper tool."""

import os
import sys
import uuid
import shutil
import logging
import argparse
import traceback
import subprocess


def exit_and_clean_up(temp_folder):
    """Log the error messages and delete the temporary folder."""
    # Capture the traceback
    logging.info("There was an unexpected failure")
    exc_type, exc_value, exc_traceback = sys.exc_info()
    for line in traceback.format_tb(exc_traceback):
        logging.info(line.encode("utf-8"))

    # Delete any files that were created for this sample
    logging.info("Removing temporary folder: " + temp_folder)
    shutil.rmtree(temp_folder)

    # Exit
    logging.info("Exit type: {}".format(exc_type))
    logging.info("Exit code: {}".format(exc_value))
    sys.exit(exc_value)


def run_cmds(commands, retry=0, catchExcept=False, stdout=None):
    """Run commands and write out the log, combining STDOUT & STDERR."""
    logging.info("Commands:")
    logging.info(' '.join(commands))
    if stdout is None:
        p = subprocess.Popen(commands,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT)
        stdout, stderr = p.communicate()
    else:
        with open(stdout, "wt") as fo:
            p = subprocess.Popen(commands,
                                 stderr=subprocess.PIPE,
                                 stdout=fo)
            stdout, stderr = p.communicate()
        stdout = False
    exitcode = p.wait()
    if stdout:
        logging.info("Standard output of subprocess:")
        for line in stdout.decode("latin-1").split('\n'):
            logging.info(line.encode("utf-8"))
    if stderr:
        logging.info("Standard error of subprocess:")
        for line in stderr.decode("latin-1").split('\n'):
            logging.info(line.encode("utf-8"))

    # Check the exit code
    if exitcode != 0 and retry > 0:
        msg = "Exit code {}, retrying {} more times".format(exitcode, retry)
        logging.info(msg)
        run_cmds(commands, retry=retry - 1)
    elif exitcode != 0 and catchExcept:
        msg = "Exit code was {}, but we will continue anyway"
        logging.info(msg.format(exitcode))
    else:
        assert exitcode == 0, "Exit code {}".format(exitcode)


def get_file_from_url(file_url, temp_folder):
    """Get a file from a URL and place it in the temporary folder."""
    logging.info("Fetching " + file_url)

    filename = file_url.split('/')[-1]
    local_path = os.path.join(temp_folder, filename)

    logging.info("Filename: " + filename)
    logging.info("Local path: " + local_path)

    # Get files from AWS S3
    if file_url.startswith('s3://'):
        logging.info("Getting reads from S3")
        run_cmds([
            'aws', 's3', 'cp', '--quiet', '--sse',
            'AES256', file_url, temp_folder
        ])

    # Get files from an FTP server
    elif file_url.startswith('ftp://'):
        logging.info("Getting reads from FTP")
        run_cmds(['wget', '-P', temp_folder, file_url])

    elif file_url.startswith(('s3://', 'ftp://')):
        logging.info("Treating as local path")
        msg = "Input file does not exist ({})".format(file_url)
        assert os.path.exists(file_url), msg
        logging.info("Making symbolic link in temporary folder")
        os.symlink(file_url, local_path)
        return local_path

    return local_path


def safe_copy_file(path_from, path_to):
    """Copy a file, either locally or to S3."""
    if path_to.startswith("s3://"):
        try:
            run_cmds(["aws", "s3", "cp", path_from, path_to])
        except:
            exit_and_clean_up(temp_folder)
    else:
        try:
            os.rename(path_from, path_to)
        except:
            exit_and_clean_up(temp_folder)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""
    Run eggNOG mapper on a set of protein sequences in FASTA format.
    """)

    parser.add_argument("--input",
                        type=str,
                        required=True,
                        help="""Location for input file.
                                (Supported: s3://, ftp://, or local path).""")
    parser.add_argument("--db",
                        type=str,
                        required=True,
                        help="""Folder containing eggNOG mapper database files.""")
    parser.add_argument("--output",
                        type=str,
                        required=True,
                        help="""Output in TSV format.""")
    parser.add_argument("--cpu",
                        type=int,
                        required=True,
                        help="""Number of CPUs to use.""")
    parser.add_argument("--temp-folder",
                        type=str,
                        default='/share',
                        help="Folder used for temporary files.")

    args = parser.parse_args()

    # Check that the temporary folder exists
    assert os.path.exists(args.temp_folder)

    # Set a random string, which will be appended to all temporary files
    random_string = str(uuid.uuid4())[:8]

    # Make a temporary folder within the --temp-folder with the random string
    temp_folder = os.path.join(args.temp_folder, str(random_string))
    # Make sure it doesn't already exist
    msg = "Collision, {} already exists".format(temp_folder)
    assert os.path.exists(temp_folder) is False, msg
    # Make the directory
    os.mkdir(temp_folder)

    # Set up logging
    log_fp = '{}/log.txt'.format(temp_folder)
    logFormatter = logging.Formatter(
        '%(asctime)s %(levelname)-8s [run_eggnog_mapper.py] %(message)s')
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

    # Get the reference database
    local_db_folder = os.path.join(temp_folder, "db") + "/"
    logging.info("Downloading the reference database from {}, writing to {}".format(
        args.db, local_db_folder
    ))
    try:
        run_cmds([
            "aws", "s3", "sync", "--quiet", args.db, local_db_folder
        ])
    except:
        exit_and_clean_up(temp_folder)

    logging.info("Processing file: " + args.input)

    local_input_file = get_file_from_url(args.input, temp_folder)

    # Run eggNOG mapper
    local_output_prefix = os.path.join(temp_folder, "output")
    try:
        run_cmds([
            "emapper.py", 
            "-i", local_input_file,
            "--output", local_output_prefix,
            "-m", "diamond",
            "--cpu", str(args.cpu),
            "--data_dir", local_db_folder,
            "--scratch_dir", temp_folder
        ])
    except:
        exit_and_clean_up(temp_folder)

    # Move the results to the final output location
    local_output_file = local_output_prefix + '.emapper.annotations'
    try:
        assert os.path.exists(local_output_file)
    except:
        exit_and_clean_up(temp_folder)

    logging.info("Copying output to " + args.output)
    safe_copy_file(local_output_file, args.output)

    logging.info("Copying logs to {}.log".format(args.output.replace(".tsv", "")))
    safe_copy_file(log_fp, args.output.replace(".tsv", "") + ".log")

    # Delete any files that were created in this process
    logging.info("Deleting temporary folder: {}".format(temp_folder))
    shutil.rmtree(temp_folder)

    # Stop logging
    logging.info("Done")
    logging.shutdown()
