#!/usr/bin/env python
from operator import itemgetter
import itertools
import os
import pandas as pd
from pymatgen.analysis.substrate_analyzer import SubstrateAnalyzer
from pymatgen import MPRester
from tqdm import tqdm
from maptool.online_extract import check_apikey


def make_connect(mpid=None,struct=None):
    mpr = check_apikey()

    # Get list of material IDs for common substrates
    mids = mpr.get_all_substrates()
    
    substrates = mpr.query({"material_id": {"$in": mids}}, ["structure", "material_id"])
    
    # Pull an example film structure. You will need to construct a
    # `pymatgen.core.structure.Structure`
    if mpid is not None:
       film = mpr.get_structure_by_material_id(mpid)
    if struct is not None:
       film = struct
    return film,substrates

def groupby_itemkey(iterable, item):
    """groupby keyed on (and pre-sorted by) itemgetter(item)."""
    itemkey = itemgetter(item)
    return itertools.groupby(sorted(iterable, key=itemkey), itemkey)

def get_subs(film,substrates):
    all_matches = []
    sa = SubstrateAnalyzer()
    for s in tqdm(substrates):
        substrate = s["structure"]
    
        # Calculate all matches and group by substrate orientation
        matches_by_orient = groupby_itemkey(
            sa.calculate(film, substrate, lowest=True),
            "sub_miller")
    
        # Find the lowest area match for each substrate orientation
        lowest_matches = [min(g, key=itemgetter("match_area"))
                          for k, g in matches_by_orient]
    
        for match in lowest_matches:
            db_entry = {
                "sub_id": s["material_id"],
                "orient": " ".join(map(str, match["sub_miller"])),
                "sub_form": substrate.composition.reduced_formula,
                "film_orient": " ".join(map(str, match["film_miller"])),
                "area": match["match_area"],
            }
            all_matches.append(db_entry)
    
    df = pd.DataFrame(all_matches)
    df.set_index("sub_id", inplace=True)
    df.sort_values("area")
    return df
   
