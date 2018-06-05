FROM quay.io/biocontainers/eggnog-mapper:1.0.3--py27_0
MAINTAINER Samuel Minot, PhD sminot@fredhutch.org

# Install BCW
RUN pip install bucket_command_wrapper==0.3.0 awscli

# Add the wrapper script
ADD run_eggnog_mapper.py /usr/local/bin/
