import numpy as np
import pandas as pd
import tfs
from cpymad.madx import Madx
from pymadng import MAD

# from turn_by_turn import read_tbt
# from turn_by_turn.structures import TbtData

madx = Madx()
madx.call("ex-fodo-acd.madx")

madx_trk_results: pd.DataFrame = madx.table.trackone.dframe()
madx_trk_results.index = madx_trk_results.index.str.replace("#s", "$start")
madx_trk_results.index = madx_trk_results.index.str.replace("#e", "$end")
madx_trk_results.index.name = "name"
madx_trk_results = madx_trk_results[madx_trk_results["turn"] != 0]

madx_tws_results = tfs.read("twiss_acd_x.tfs", index="NAME")

mad = MAD()
mad.loadfile("ex-fodo-acd.mad")

mad_trk_results = tfs.read("track_acd_n.tfs", index="name")
mad_trk_results.index = mad_trk_results.index.str.lower()

mad_tws_results = tfs.read("twiss_acd_n.tfs", index="name")

# Convert the column names from MAD-NG to MADX
mad_tws_results = mad_tws_results.rename(
    columns={
        "beta11": "BETX",
        "beta22": "BETY",
        "mu1": "MUX",
        "mu2": "MUY",
        "alfa11": "ALFX",
        "alfa22": "ALFY",
    }
)

print("Average difference of MAD-NG vs MADX comparison for twiss results:")
for col in ["BETX", "BETY", "MUX", "MUY", "ALFX", "ALFY"]:
    mean_diff = np.mean(
        np.abs(mad_tws_results[col].values - madx_tws_results[col].values)
    )
    print(f"{col}: {mean_diff}")

print("\nAverage difference of MAD-NG vs MADX comparison for track results:")
for col in ["x", "px", "y", "py"]:
    mad_trk_results[col] = mad_trk_results[col]
    madx_trk_results[col] = madx_trk_results[col]
    print(f"{col}: {np.mean(np.abs(mad_trk_results[col] - madx_trk_results[col]))}")
