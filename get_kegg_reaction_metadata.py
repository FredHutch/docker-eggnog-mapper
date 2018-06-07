#!/usr/bin/env python3

import argparse
import os
import json
import gzip
import sqlite3
import requests
from collections import defaultdict
import logging
from multiprocessing import Pool
from functools import partial
import re

# Set up logging
logFormatter = logging.Formatter(
    '%(asctime)s %(levelname)-8s [kegg_reaction_metadata] %(message)s'
)
rootLogger = logging.getLogger()
rootLogger.setLevel(logging.INFO)

# Also write to STDOUT
consoleHandler = logging.StreamHandler()
consoleHandler.setFormatter(logFormatter)
rootLogger.addHandler(consoleHandler)


def chunks(l, n):
    for i in range(0, len(l), n):
        yield l[i:i + n]


def sql_safe_string(s):
    for c in ["[", "]", ":", ";", ",", "'", '"']:
        s = s.replace(c, "")
    return s


def open_tsv(fp, skip=0):
    if fp.endswith(".gz"):
        f = gzip.open(fp, "rt")
    else:
        f = open(fp, "rt")

    for ix, line in enumerate(f):
        if ix >= skip:
            yield(line.rstrip("\n").split("\t"))

    f.close()


def fetch_kegg_api(kegg_id, data_type):
    r = requests.get("http://rest.kegg.jp/get/{}:{}".format(data_type, kegg_id))

    data = defaultdict(list)

    line_label = None
    for line in r.text.split("\n"):
        if len(line[:12].rstrip(" ")) > 0:
            line_label = line[:12].strip(" ")
        line_value = line[12:]

        data[line_label].append(line_value)

    return data


def get_kegg_reaction_metadata(input_tsv=None, output_db=None, threads=1, chunk_size=100):
    """Get reaction metadata for the KEGG entries from eggNOG output, write to SQLite."""

    # Make sure the tables exist for orthology, reaction, pathway, and compound
    conn = sqlite3.connect(output_db)

    c = conn.cursor()
    c.execute(
        "create table if not exists ortholog(ortholog_id TEXT, name TEXT, definition TEXT);"
    )
    c.execute(
        "create unique index if not exists ortholog_ix on ortholog (ortholog_id);"
    )
    c.execute(
        """create table if not exists ortholog_has_reaction
        (ortholog_id TEXT, reaction_id TEXT, UNIQUE(ortholog_id, reaction_id))
        ;"""
    )
    c.execute(
        """create table if not exists reaction
        (reaction_id TEXT PRIMARY KEY, definition TEXT, equation TEXT, enzyme TEXT,
        direction INT);"""
    )
    c.execute(
        "create unique index if not exists reaction_ix on reaction (reaction_id);"
    )
    c.execute(
        """create table if not exists reaction_compound
        (reaction_id TEXT, compound_id TEXT, stoichiometry INT, side TEXT,
        UNIQUE(compound_id, reaction_id))
        ;"""
    )
    c.execute(
        "create table if not exists compound(compound_id PRIMAY KEY TEXT, formula TEXT);"
    )
    c.execute(
        """CREATE TABLE IF NOT EXISTS compound_name
            (compound_id TEXT, name TEXT, UNIQUE(compound_id, name));
        """
    )
    c.execute(
        "create unique index if not exists compound_ix on compound (compound_id);"
    )
    c.execute(
        """create table if not exists query_ortholog
        (query_id TEXT, ortholog_id TEXT, UNIQUE(ortholog_id, query_id))
        ;"""
    )

    # Pathway, modules and links tables
    c.execute(
        """create table if not exists pathway
        (pathway_id TEXT PRIMARY_KEY,
        name TEXT,
        class TEXT,
        description TEXT);"""
    )
    c.execute(
        "create unique index if not exists pathway_ix on pathway (pathway_id);"
    )
    c.execute(
        """CREATE TABLE IF NOT EXISTS module
        (
            module_id TEXT PRIMARY_KEY,
            name TEXT,
            class TEXT
        ); """
    )
    c.execute("create unique index if not exists module_ix on module (module_id);")
    c.execute(
        """create table if not exists module_pathway
        (module_id TEXT, pathway_id TEXT, UNIQUE(module_id, pathway_id))
        ;"""
    )
    c.execute(
        """create table if not exists module_compound
        (module_id TEXT, compound_id TEXT, UNIQUE(module_id, compound_id))
        ;"""
    )
    c.execute(
        """create table if not exists module_reaction
        (module_id TEXT, reaction_id TEXT, UNIQUE(module_id, reaction_id))
        ;"""
    )
    c.execute(
        """create table if not exists pathway_compound
        (pathway_id TEXT, compound_id TEXT, UNIQUE(pathway_id, compound_id))
        ;"""
    )
    c.execute(
        """create table if not exists pathway_reaction
        (pathway_id TEXT, reaction_id TEXT, UNIQUE(pathway_id, reaction_id))
        ;"""
    )
    conn.commit()

    # Get the set of KEGG IDs in the input TSV
    kegg_ids = set([])
    # and mapping of query_name to keggs
    query_keggs = defaultdict(set)

    f = open_tsv(input_tsv, skip=3)

    header = next(f)

    # Get the set of KEGG ids represented here
    for line in f:
        if len(line) < len(header):
            continue
        entry = dict(zip(header, line))
        assert "KEGG_KOs" in entry, entry
        assert "#query_name" in entry, entry
        kos = entry.get("KEGG_KOs")
        query_name = entry.get("#query_name")
        l_ko = {ko.strip() for ko in kos.split(',') if not ko.strip() == ''}
        kegg_ids.update(l_ko)
        query_keggs[query_name].update(l_ko)

    logging.info("There are {:,} KEGG IDs in the input TSV".format(len(kegg_ids)))
    c.executemany(
        """INSERT OR REPLACE INTO query_ortholog
        (query_id, ortholog_id)
        VALUES (?, ?)
        """,
        [
            (query, q_kegg)
            for query, q_keggs in query_keggs.items()
            for q_kegg in q_keggs
        ]
    )
    conn.commit()

    # Figure out if we have any existing KEGGS in our db
    c.execute('select ortholog_id from ortholog;')
    existing_orthologs = {i[0] for i in c.fetchall()}
    logging.info("There are {:,} KEGG IDs already in the database".format(
        len(existing_orthologs)))

    kegg_ids_needed = list(kegg_ids - existing_orthologs)
    logging.info("There are {:,} KEGG IDs we need to download".format(
        len(kegg_ids_needed)))

    # Begin downloading from KEGG
    kegg_pool = Pool(threads)

    reaction_ids = set()
    logging.info("Starting download of KEGG data")
    # The raw data by kegg_id
    for kegg_id_chunk_i, kegg_id_chunk in enumerate(chunks(kegg_ids_needed, chunk_size)):
        logging.info("Working on KEGG ID chunk {:,} of {:,}".format(
            kegg_id_chunk_i + 1,
            int(len(kegg_ids_needed) / chunk_size) + 1
        ))
        raw_kegg_data = kegg_pool.map(
            partial(
                fetch_kegg_api,
                data_type='ko'),
            list(kegg_id_chunk)
            )
        logging.info("Completed downloading base data for chunk {}".format(
            kegg_id_chunk_i + 1,
        ))
        logging.info("Inserting orthologs")
        try:
            c.executemany(
                """INSERT OR REPLACE INTO ORTHOLOG
                    (ortholog_id, name, definition)
                    VALUES (?, ?, ?);""",
                [(
                    sql_safe_string(kegg_id),
                    sql_safe_string(kegg_data.get("NAME", [""])[0]),
                    sql_safe_string(kegg_data.get("DEFINITION", [""])[0]),
                ) for (kegg_id, kegg_data) in zip(kegg_id_chunk, raw_kegg_data)]
            )
            conn.commit()
        except Exception as e:
            logging.error("Insertion of orthologs failed with exception {}".format(e))
            assert 1 == 0
        # Update our reaction id set
        ortholog_reactions = {
            kegg_id: set(dblink[4:].split(' '))
            for (kegg_id, kegg_data) in zip(kegg_id_chunk, raw_kegg_data)
            for dblink in kegg_data.get("DBLINKS", []) if dblink.startswith('RN: ')
        }
        logging.info("Inserting ortholog-reaction-links")
        for ortholog, o_rxns in ortholog_reactions.items():
            try:
                c.executemany(
                    """INSERT OR REPLACE INTO ortholog_has_reaction
                        (ortholog_id, reaction_id)
                        VALUES (?, ?)
                    """,
                    [
                        (sql_safe_string(ortholog), sql_safe_string(rxn_id))
                        for rxn_id in o_rxns
                    ]
                )
            except Exception as e:
                logging.error("Insertion of ortholog-reactions failed with exception {}".format(e))
                assert 1 == 0
        for rxns in ortholog_reactions.values():
            reaction_ids.update(rxns)

    # Now work on reactions
    c.execute("select reaction_id from reaction;")
    existing_reaction_ids = {i[0] for i in c.fetchall()}
    c.execute("SELECT DISTINCT reaction_id FROM ortholog_has_reaction")
    linked_reactions = {i[0] for i in c.fetchall()}
    reaction_ids_needed = list(reaction_ids.union(linked_reactions) - existing_reaction_ids)
    logging.info("There are {:,} reactions we need to retrieve from KEGG".format(
        len(reaction_ids_needed))
        )

    re_compound = re.compile(r'((?P<stoich>\d+)\s+|)(?P<compound_id>(?P<type>G|C)\d+)')
    pathways = set()
    compound_ids = set()
    for rxn_id_chunk_i, rxn_id_chunk in enumerate(chunks(reaction_ids_needed, chunk_size)):
        logging.info("Working on reaction chunk {:,} of {:,}".format(
            rxn_id_chunk_i + 1,
            int(len(reaction_ids_needed) / chunk_size) + 1
        ))
        raw_rxn_data = kegg_pool.map(
            partial(
                fetch_kegg_api,
                data_type='rn'),
            list(rxn_id_chunk)
            )
        logging.info("Inserting {} reactions".format(
            len(raw_rxn_data)
        ))
        try:
            c.executemany(
                """INSERT OR REPLACE INTO reaction
                    (reaction_id, definition, equation, enzyme)
                    VALUES (?, ?, ?, ?);
                """,
                [(
                    sql_safe_string(rxn_id),
                    sql_safe_string(rxn_data.get("DEFINITION", [""])[0]),
                    sql_safe_string(rxn_data.get("EQUATION", [""])[0]),
                    sql_safe_string(rxn_data.get("ENZYME", [""])[0]),
                ) for (rxn_id, rxn_data) in zip(rxn_id_chunk, raw_rxn_data)
                ]
            )
            conn.commit()
        except Exception as e:
            logging.error("Insertion of reactions failed with exception {}".format(e))
            assert 1 == 0

        logging.info("Parsing reaction equations")
        for rxn_id, rrd in zip(rxn_id_chunk, raw_rxn_data):
            if 'EQUATION' not in rrd:
                logging.warning("{} has no equation".format(
                    rxn_id
                ))
            # implicit else
            e = rrd['EQUATION'][0]
            # Split the equation by =
            e_L, e_R = e.split("=")
            # Determine directionality
            if e_L[-1] == "<" and e_R[0] == ">":
                rrd['direction'] = 2
            elif e_L[-1] == "<":
                rrd['direction'] = 1
            elif e_R[0] == ">":
                rrd['direction'] = 0
            else:
                rrd['direction'] = -1
            c.execute(
                "UPDATE reaction SET direction = ? WHERE reaction_id ==?",
                (rrd['direction'], rxn_id)
                )
            # And compounds
            L_compounds = [re_compound.search(c.strip()) for c in e_L[:-1].split(' + ')]
            R_compounds = [re_compound.search(c.strip()) for c in e_R[1:].split(' + ')]
            logging.info("Inserting {:,} reaction-compound links".format(
                len(L_compounds)+len(R_compounds)
            ))
            try:
                c.executemany(
                    """INSERT OR REPLACE INTO reaction_compound
                        (reaction_id, compound_id, side, stoichiometry)
                        VALUES (?, ?, ?, ?)
                    """,
                    [
                        (
                            sql_safe_string(rxn_id),
                            sql_safe_string(m.group('compound_id')),
                            sql_safe_string("L"),
                            sql_safe_string(m.group('stoich'))
                        ) if m.group('stoich') != None
                        else
                        (
                            sql_safe_string(rxn_id),
                            sql_safe_string(m.group('compound_id')),
                            sql_safe_string("L"),
                            sql_safe_string(str(1))
                        )
                        for m in L_compounds
                    ]+[
                        (
                            sql_safe_string(rxn_id),
                            sql_safe_string(m.group('compound_id')),
                            sql_safe_string("R"),
                            sql_safe_string(m.group('stoich'))
                        ) if m.group('stoich') != None
                        else
                        (
                            sql_safe_string(rxn_id),
                            sql_safe_string(m.group('compound_id')),
                            sql_safe_string("R"),
                            sql_safe_string(str(1))
                        )
                        for m in R_compounds
                    ]
                )
                conn.commit()
            except Exception as e:
                logging.error("Insertion of reaction-compound failed with exception {}".format(e))
                assert 1 == 0
            compound_ids.update({m.group('compound_id') for m in L_compounds+R_compounds})

            rxn_pathway = [
                (sql_safe_string(rxn_id), sql_safe_string(p.split(" ", 1)[0]))
                for rxn_id, rxn in zip(rxn_id_chunk, raw_rxn_data)
                for p in rxn['PATHWAY']
            ]
            logging.info("Now working on {:,} reaction-pathway links".format(
                len(rxn_pathway)
            ))
            pathways.update({rp[1] for rp in rxn_pathway})
            try:
                c.executemany(
                    """INSERT OR REPLACE INTO pathway_reaction
                        (reaction_id,pathway_id)
                        VALUES (?,?)""",
                    rxn_pathway
                )
                conn.commit()
            except Exception as e:
                logging.error("Insertion of reaction-pathway failed with exception {}".format(e))
                assert 1 == 0

    # Pathways
    c.execute("SELECT pathway_id from pathway;")
    existing_pathways = {i[0] for i in c.fetchall()}
    c.execute("SELECT DISTINCT pathway_id from pathway_reaction;")
    pathway_links = {i[0] for i in c.fetchall()}
    pathways_needed = list(pathways.union(pathway_links) - existing_pathways)
    logging.info("We need to download {:,} pathways".format(len(pathways_needed)))
    for id_chunk_i, id_chunk in enumerate(chunks(pathways_needed, chunk_size)):
        logging.info("Working on pathway chunk {:,} of {:,}".format(
            id_chunk_i + 1,
            int(len(pathways_needed) / chunk_size) + 1
        ))
        raw_data = kegg_pool.map(
            partial(
                fetch_kegg_api,
                data_type='path'),
            list(id_chunk)
            )
        logging.info("Inserting {:,} pathways".format(len(raw_data)))
        c.executemany(
            """INSERT OR REPLACE INTO pathway
            (pathway_id, name, class, description)
            VALUES (?, ?, ?, ?)
            """,
            [(
                sql_safe_string(path_id),
                sql_safe_string(path_data.get('NAME', [''])[0]),
                sql_safe_string(path_data.get('CLASS', [''])[0]),
                sql_safe_string(path_data.get('DESCRIPTION', [''])[0]),
            )
                for path_id, path_data in
                zip(id_chunk, raw_data)]
        )
        c.executemany(
            """INSERT OR REPLACE INTO pathway_reaction
            (pathway_id, reaction_id)
            VALUES (?, ?)
            """,
            [(
                sql_safe_string(path_id),
                sql_safe_string(rxn)
            )
                for path_id, path_data in
                zip(id_chunk, raw_data)
                for rxn_block in path_data['REACTION']
                for rxn in rxn_block.split(' ')[0].split(',')
                ]
        )
        c.executemany(
            """INSERT OR REPLACE INTO pathway_compound
            (pathway_id, compound_id)
            VALUES (?, ?)
            """,
            [(
                sql_safe_string(path_id),
                sql_safe_string(comp.split(" ")[0].strip())
            )
                for path_id, path_data in
                zip(id_chunk, raw_data)
                for comp in path_data['COMPOUND']
                ]
        )
        c.executemany(
            """INSERT OR REPLACE INTO module_pathway
            (pathway_id, module_id)
            VALUES (?, ?)
            """,
            [(
                sql_safe_string(path_id),
                sql_safe_string(mod.split(" ")[0].strip())
            )
                for path_id, path_data in
                zip(id_chunk, raw_data)
                for mod in path_data['MODULE']
                ]
        )

        conn.commit()

    # Modules
    c.execute("SELECT module_id from module;")
    existing_modules = {i[0] for i in c.fetchall()}
    c.execute("SELECT DISTINCT module_id from module_pathway;")
    module_links = {i[0] for i in c.fetchall()}
    modules_needed = list(module_links - existing_modules)
    logging.info("We need to download {:,} modules".format(len(modules_needed)))
    for id_chunk_i, id_chunk in enumerate(chunks(modules_needed, chunk_size)):
        logging.info("Working on module chunk {:,} of {:,}".format(
            id_chunk_i + 1,
            int(len(modules_needed) / chunk_size) + 1
        ))
        raw_data = kegg_pool.map(
            partial(
                fetch_kegg_api,
                data_type='md'),
            list(id_chunk)
            )
        logging.info("Inserting {:,} modules".format(len(raw_data)))
        c.executemany(
            """INSERT OR REPLACE INTO module
            (module_id, name, class)
            VALUES (?,?,?)
            """,
            [(
                sql_safe_string(mod_id),
                sql_safe_string(mod_data.get('NAME', [''])[0]),
                sql_safe_string(mod_data.get('CLASS', [''])[0]),
            )
                for mod_id, mod_data in
                zip(id_chunk, raw_data)]
        )
        c.executemany(
            """INSERT OR REPLACE INTO module_compound
            (module_id, compound_id)
            VALUES (?,?)
            """,
            [(
                sql_safe_string(mod_id),
                sql_safe_string(comp.split(" ")[0].strip())
            )
                for mod_id, mod_data in
                zip(id_chunk, raw_data)
                for comp in mod_data['COMPOUND']
                ]
        )
        c.executemany(
            """INSERT OR REPLACE INTO module_reaction
            (module_id, reaction_id)
            VALUES (?,?)
            """,
            [(
                sql_safe_string(mod_id),
                sql_safe_string(rxn)
            )
                for mod_id, mod_data in
                zip(id_chunk, raw_data)
                for rxn_block in mod_data['REACTION']
                for rxn in rxn_block.split(' ')[0].split(',')
                ]
        )
        conn.commit()

    # Revisit reactions to capture everything in our modules / pathways
    # (but no further spinning out from here)
    c.execute("SELECT DISTINCT reaction_id FROM pathway_reaction;")
    rxn_link_path = {i[0] for i in c.fetchall()}
    c.execute("SELECT DISTINCT reaction_id FROM module_reaction;")
    rxn_link_mod = {i[0] for i in c.fetchall()}
    c.execute("SELECT reaction_id from reaction;")
    existing_rxn = {i[0] for i in c.fetchall()}
    rxn_needed = list(rxn_link_mod.union(rxn_link_path) - existing_rxn)
    logging.info("Downloading {:,} more reactions from pathways and modules".format(
        len(rxn_needed)))
    for rxn_id_chunk_i, rxn_id_chunk in enumerate(chunks(rxn_needed, chunk_size)):
        logging.info("Working on reaction chunk {:,} of {:,}".format(
            rxn_id_chunk_i + 1,
            int(len(rxn_needed) / chunk_size) + 1
        ))
        raw_rxn_data = kegg_pool.map(
            partial(
                fetch_kegg_api,
                data_type='rn'),
            list(rxn_id_chunk)
            )
        logging.info("Inserting {} reactions".format(
            len(raw_rxn_data)
        ))
        try:
            c.executemany(
                """INSERT OR REPLACE INTO reaction
                    (reaction_id, definition, equation, enzyme)
                    VALUES (?, ?, ?, ?);
                """,
                [(
                    sql_safe_string(rxn_id),
                    sql_safe_string(rxn_data.get("DEFINITION", [""])[0]),
                    sql_safe_string(rxn_data.get("EQUATION", [""])[0]),
                    sql_safe_string(rxn_data.get("ENZYME", [""])[0]),
                ) for (rxn_id, rxn_data) in zip(rxn_id_chunk, raw_rxn_data)
                ]
            )
            conn.commit()
        except Exception as e:
            logging.error("Insertion of reactions failed with exception {}".format(e))
            assert 1 == 0

        logging.info("Parsing reaction equations")
        for rxn_id, rrd in zip(rxn_id_chunk, raw_rxn_data):
            if 'EQUATION' not in rrd:
                logging.warning("{} has no equation".format(
                    rxn_id
                ))
                continue
            # Implicit else
            e = rrd['EQUATION'][0]
            # Split the equation by =
            e_L, e_R = e.split("=")
            # Determine directionality
            if e_L[-1] == "<" and e_R[0] == ">":
                rrd['direction'] = 2
            elif e_L[-1] == "<":
                rrd['direction'] = 1
            elif e_R[0] == ">":
                rrd['direction'] = 0
            else:
                rrd['direction'] = -1
            c.execute(
                "UPDATE reaction SET direction = ? WHERE reaction_id ==?",
                (rrd['direction'], rxn_id)
                )
            # And compounds
            L_compounds = [re_compound.search(c.strip()) for c in e_L[:-1].split(' + ')]
            R_compounds = [re_compound.search(c.strip()) for c in e_R[1:].split(' + ')]
            logging.info("Inserting {:,} reaction-compound links".format(
                len(L_compounds)+len(R_compounds)
            ))
            try:
                c.executemany(
                    """INSERT OR REPLACE INTO reaction_compound
                        (reaction_id, compound_id, side, stoichiometry)
                        VALUES (?, ?, ?, ?)
                    """,
                    [
                        (
                            sql_safe_string(rxn_id),
                            sql_safe_string(m.group('compound_id')),
                            sql_safe_string("L"),
                            sql_safe_string(m.group('stoich'))
                        ) if m.group('stoich') != None
                        else
                        (
                            sql_safe_string(rxn_id),
                            sql_safe_string(m.group('compound_id')),
                            sql_safe_string("L"),
                            sql_safe_string(str(1))
                        )
                        for m in L_compounds
                    ]+[
                        (
                            sql_safe_string(rxn_id),
                            sql_safe_string(m.group('compound_id')),
                            sql_safe_string("R"),
                            sql_safe_string(m.group('stoich'))
                        ) if m.group('stoich') != None
                        else
                        (
                            sql_safe_string(rxn_id),
                            sql_safe_string(m.group('compound_id')),
                            sql_safe_string("R"),
                            sql_safe_string(str(1))
                        )
                        for m in R_compounds
                    ]
                )
                conn.commit()
            except Exception as e:
                logging.error("Insertion of reaction-compound failed with exception {}".format(e))
                assert 1 == 0

    # Now compounds
    c.execute("select compound_id from compound;")
    existing_compound_ids = {i[0] for i in c.fetchall()}
    c.execute("SELECT DISTINCT compound_id from reaction_compound;")
    linked_compound_ids = {i[0] for i in c.fetchall()}
    compound_ids_needed = list(linked_compound_ids - existing_compound_ids)
    compound_ids_needed_c = [c for c in compound_ids_needed if c[0] == 'C']
    # Also account for glycans
    compound_ids_needed_g = [c for c in compound_ids_needed if c[0] == 'G']
    logging.info("We need to download {:,} compounds".format(len(compound_ids_needed)))

    for comp_id_chunk_i, comp_id_chunk in enumerate(chunks(compound_ids_needed_c, chunk_size)):
        logging.info("Working on compound chunk {:,} of {:,}".format(
            comp_id_chunk_i + 1,
            int(len(compound_ids_needed) / chunk_size) + 1
        ))
        raw_comp_data = kegg_pool.map(
            partial(
                fetch_kegg_api,
                data_type='cpd'),
            list(comp_id_chunk)
            )
        c.executemany(
            """INSERT OR REPLACE INTO compound
                (compound_id, formula)
                VALUES (?, ?)
            """,
            [(comp_id, comp_data.get('FORMULA', [""])[0])
             for comp_id, comp_data in zip(comp_id_chunk, raw_comp_data)]
        )
        c.executemany(
            """INSERT OR REPLACE INTO compound_name
            (compound_id, name)
            VALUES (?,?)
            """,
            [(comp_id, name.replace(";", ""))
                for comp_id, comp_data in zip(comp_id_chunk, raw_comp_data)
                for name in comp_data['NAME']]
        )
        conn.commit()
    # Glycans
    for comp_id_chunk_i, comp_id_chunk in enumerate(chunks(compound_ids_needed_g, chunk_size)):
        logging.info("Working on glycan chunk {:,} of {:,}".format(
            comp_id_chunk_i + 1,
            int(len(compound_ids_needed) / chunk_size) + 1
        ))
        raw_comp_data = kegg_pool.map(
            partial(
                fetch_kegg_api,
                data_type='gl'),
            list(comp_id_chunk)
            )
        c.executemany(
            """INSERT OR REPLACE INTO compound
                (compound_id, formula)
                VALUES (?, ?)
            """,
            [(comp_id, comp_data.get('COMPOSITION', [""])[0])
             for comp_id, comp_data in zip(comp_id_chunk, raw_comp_data)]
        )
        c.executemany(
            """INSERT OR REPLACE INTO compound_name
            (compound_id, name)
            VALUES (?,?)
            """,
            [(comp_id, name.replace(";", ""))
                for comp_id, comp_data in zip(comp_id_chunk, raw_comp_data)
                for name in comp_data['NAME']]
        )
        conn.commit()

    conn.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""
    Get reaction metadata for the KEGG entries from eggNOG output, write to SQLite.
    """)

    parser.add_argument("--input-tsv",
                        type=str,
                        required=True,
                        help="""Location for local input path.""")
    parser.add_argument("--threads",
                        type=int,
                        default=1,
                        help="""Number of concurrent calls to KEGG""")
    parser.add_argument("--chunk-size",
                        type=int,
                        default=100,
                        help="""Number of KEGG ids to work on at a time""")
    parser.add_argument("--output-db",
                        type=str,
                        required=True,
                        help="""SQLite database output path.""")

    args = parser.parse_args()

    get_kegg_reaction_metadata(**args.__dict__)
