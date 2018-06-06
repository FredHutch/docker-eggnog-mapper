#!/usr/bin/env python3

import os
import json
import gzip
import sqlite3
import requests
import argparse
from collections import defaultdict


def open_tsv(fp, skip=0):
    if fp.endswith(".gz"):
        f = gzip.open(fp, "rt")
    else:
        f = open(fp, "rt")

    for ix, line in enumerate(f):
        if ix >= skip:
            yield(line.rstrip("\n").split("\t"))

    f.close()


def fetch_kegg_api(data_type, kegg_id):
    r = requests.get("http://rest.kegg.jp/get/{}:{}".format(data_type, kegg_id))

    data = defaultdict(list)

    line_label = None
    for line in r.text.split("\n"):
        if len(line[:12].rstrip(" ")) > 0:
            line_label = line[:12].strip(" ")
        line_value = line[12:]

        data[line_label].append(line_value)

    return data

def get_kegg_reaction_metadata(input_tsv=None, output_db=None):
    """Get reaction metadata for the KEGG entries from eggNOG output, write to SQLite."""

    # Make sure the tables exist for orthology, reaction, pathway, and compound

    conn = sqlite3.connect(output_db)

    c = conn.cursor()
    c.execute(
        "create table if not exists ortholog(ortholog TEXT, name TEXT, definition TEXT);"
    )
    c.execute(
        "create unique index if not exists ortholog_ix on ortholog (ortholog);"
    )
    c.execute(
        "create table if not exists reaction(reaction TEXT, ortholog TEXT, definition TEXT, equation TEXT, enzyme TEXT);"
    )
    c.execute(
        "create unique index if not exists reaction_ix on reaction (reaction);"
    )
    c.execute(
        "create table if not exists compound(compound TEXT, name TEXT, formula TEXT);"
    )
    c.execute(
        "create unique index if not exists compound_ix on compound (compound);"
    )
    c.execute(
        "create table if not exists pathway(pathway TEXT, name TEXT, class TEXT);"
    )
    c.execute(
        "create unique index if not exists pathway_ix on pathway (pathway);"
    )
    conn.commit()

    # Get the set of KEGG IDs in the input TSV
    kegg_ids = set([])

    f = open_tsv(input_tsv, skip=3)

    header = next(f)

    for line in f:
        if len(line) < len(header):
            continue
        entry = dict(zip(header, line))
        assert "KEGG_KOs" in entry, entry
        kos = entry["KEGG_KOs"]
        if len(kos) > 0:
            for ko in kos.split(","):
                kegg_ids.add(ko)

    print("There are {:,} KEGG IDs in the input TSV".format(len(kegg_ids)))

    for kegg_id in list(kegg_ids):
        print(kegg_id)
        
        kegg_data = fetch_kegg_api("ko", kegg_id)

        sql_cmd = "insert or replace into ortholog (ortholog, name, definition) values ('{}', '{}', '{}');".format(
            sql_safe_string(kegg_id), 
            sql_safe_string(kegg_data["NAME"][0]),
            sql_safe_string(kegg_data["DEFINITION"][0]),
        )
        try:
            c.execute(sql_cmd)
        except:
            print(sql_cmd)
            assert 1 == 0

        # Check to see if there is any reaction data for this KEGG
        for reaction_id in [
            reaction_id
            for db_link in kegg_data.get("DBLINKS", [])
            for reaction_id in db_link[4:].split(" ")
            if db_link.startswith("RN: ")
        ]:
            reaction_data = fetch_kegg_api("rn", reaction_id)
            
            sql_cmd = "insert or replace into reaction(reaction, ortholog, definition, equation, enzyme) values ('{}', '{}', '{}', '{}', '{}');".format(
                sql_safe_string(reaction_id),
                sql_safe_string(kegg_id),
                sql_safe_string(reaction_data.get("DEFINITION", [''])[0]),
                sql_safe_string(reaction_data.get("EQUATION", [''])[0]),
                sql_safe_string(reaction_data.get("ENZYME", [''])[0]),
            )
            
            try:
                c.execute(sql_cmd)
            except:
                print(sql_cmd)
                assert 1 == 0


            for pathway_id in [
                p.split(" ", 1)[0]
                for p in reaction_data["PATHWAY"]
            ]:
                print(pathway_id)
                pathway_data = fetch_kegg_api("path", pathway_id)
                
                sql_cmd = "insert or replace into pathway(pathway, name, class) values ('{}', '{}', '{}');".format(
                    sql_safe_string(pathway_id),
                    sql_safe_string(pathway_data["NAME"][0]),
                    sql_safe_string(pathway_data.get("CLASS", [''])[0]),
                )
                try:
                    c.execute(sql_cmd)
                except:
                    print(sql_cmd)
                    assert 1 == 0


            for compound_id in [
                c
                for e in reaction_data["EQUATION"]
                for c in e.split(" ")
                if c.startswith("C")
            ]:
                print(compound_id)
                compound_data = fetch_kegg_api("cpd", compound_id)

                for compound_name in compound_data["NAME"]:
                    sql_cmd = "insert or replace into compound(compound, name, formula) values ('{}', '{}', '{}');".format(
                        sql_safe_string(compound_id),
                        sql_safe_string(compound_name),
                        sql_safe_string(compound_data.get("FORMULA", [''])[0]),
                    )
                try:
                    c.execute(sql_cmd)
                except:
                    print(sql_cmd)
                    assert 1 == 0

    conn.commit()
    conn.close()

def sql_safe_string(s):
    for c in ["[", "]", ":", ";", ",", "'", '"']:
        s = s.replace(c, "")
    return s


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""
    Get reaction metadata for the KEGG entries from eggNOG output, write to SQLite.
    """)

    parser.add_argument("--input-tsv",
                        type=str,
                        required=True,
                        help="""Location for local input path.""")
    parser.add_argument("--output-db",
                        type=str,
                        required=True,
                        help="""SQLite database output path.""")

    args = parser.parse_args()

    get_kegg_reaction_metadata(**args.__dict__)
