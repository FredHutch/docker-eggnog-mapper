FROM quay.io/biocontainers/eggnog-mapper:1.0.3--py27_0
MAINTAINER Samuel Minot, PhD sminot@fredhutch.org

# Install BCW
RUN pip install --upgrade bucket_command_wrapper==0.3.0 awscli pandas numpy feather-format

# Add the wrapper scripts
ADD run_eggnog_mapper.py /usr/local/bin/
ADD make_eggnog_abundance_dataframe.py /usr/local/bin/
