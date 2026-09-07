# Script 04b - Determine FIM per reach
#
# Reads every segment gpkg in the folder and writes ONE plan CSV
# (raster_write_plan.csv).
#
# WHAT EACH GPKG IS
#   Each hydraulic_results_*.gpkg covers one mainstem. Its filename names two
#   reaches (wb-<id_a>_wb-<id_b>) and a flow range (<lo>-cfs_to_<hi>-cfs). The
#   gpkg holds one row per reach, and each reach row has several flow_<n>
#   columns giving the modeled flow (cfs) for that run. Every planned output is
#   one raster per (reach, flow): wsel_<flowpath>_<flow>.tif.
#
# RULE 1 - keep a reach only if it belongs to this segment's span.
#   A segment carries an extra upstream "boundary" reach beyond the two named
#   in its filename; that spillover reach should come from its OWN segment, not
#   this one, so we drop it here.
#
#   HOW IT'S DECIDED : we do NOT look up true span membership. We
#   take the two filename ids as a numeric range [min(id_a,id_b), max(...)] and
#   keep a reach only if its integer id falls inside that range. Reaches
#   outside the range are marked SKIP_out_of_span.
#
#   This is only equivalent to real membership if wb-#### ids increase steadily
#   along the reach. If they don't, a boundary reach whose id happens to land
#   inside the range will be kept by mistake, and a real in-span reach whose id
#   lands outside will be dropped. Verify against the data before trusting it.
#
# RULE 2 - drop duplicates across segments, keeping the lower-flow copy.
#   Adjacent segments share a boundary flow, so the same (reach, flow) output
#   can be produced by two in-span segments. Among the kept (in-span) rows we
#   sort by the segment's upper flow bound (seg_q_hi) ascending and keep the
#   first occurrence of each output_name; later ones are marked SKIP_duplicate.
#   Result: the copy from the segment with the LOWER flow range wins.
#
# OUTPUT
#   Every row lands in the CSV with a "decision" column of WRITE,
#   SKIP_out_of_span, or SKIP_duplicate, so nothing is hidden - you can see why
#   each candidate was kept or dropped. A per-reach summary prints at the end.

# Created by: Andy Carter, PE
# Created - 2026.09.07

# ************************************************************
import os
import re
import glob
import geopandas as gpd
import pandas as pd

#import concurrent.futures
import argparse
import time
import datetime
import warnings
# ************************************************************


# ----------------
def fn_str_to_bool(value):
    if isinstance(value, bool):
        return value
    if value.lower() in {'true', 't', '1'}:
        return True
    elif value.lower() in {'false', 'f', '0'}:
        return False
    else:
        raise argparse.ArgumentTypeError(f"Boolean value expected. Got '{value}'.")
# ----------------


# ......................
def fn_reach_id(str_fp):
    m = re.search(r"(\d+)", str_fp)
    return int(m.group(1)) if m else None
# ......................


# ---------------------------------------------
def fn_determine_fims_per_reach(str_input_dir,
                                b_print_output):

    
    # supress all warnings
    warnings.filterwarnings("ignore", category=UserWarning )
    
    if b_print_output:
        print(" ")
        print("+=================================================================+")
        print("|           DETERMINE VALID FIMS PER EACH STREAM REACH            |")
        print("|                Created by Andy Carter, PE of                    |")
        print("|             Center for Water and the Environment                |")
        print("|                 University of Texas at Austin                   |")
        print("+-----------------------------------------------------------------+")
    
        
        print("  ---(i) INPUT HYDRAULIC RESULTS GEOPACKAGES (GPKG): " + str_input_dir)
        print("  ---[r] PRINT OUTPUT: " + str(b_print_output))
        print("===================================================================")
        
    str_gpkg_folder = str_input_dir
    
    str_layer = "01_flowpath_flooded_cells_ar"
    str_pattern = "*.gpkg"
    
    # pull the two span ids and the flow range out of the filename
    re_name = re.compile(r"wb-(\d+)_wb-(\d+).*?_(\d+)-cfs_to_(\d+)-cfs", re.IGNORECASE)

    rows = []
    for str_path in sorted(glob.glob(os.path.join(str_gpkg_folder, str_pattern))):
        str_name = os.path.basename(str_path)
        m = re_name.search(str_name)
        if not m:
            print("  !! filename not parsed, skipping:", str_name)
            continue
        id_a, id_b, q_lo, q_hi = int(m.group(1)), int(m.group(2)), int(m.group(3)), int(m.group(4))
        span_lo, span_hi = min(id_a, id_b), max(id_a, id_b)
    
        try:
            gdf = gpd.read_file(str_path, layer=str_layer)
        except Exception as e:
            print("  !! could not read", str_name, "-", e)
            continue
        flow_cols = [c for c in gdf.columns if c.startswith("flow_")]
    
        for _, row in gdf.iterrows():
            str_fp = str(row["flowpath"])
            rid = fn_reach_id(str_fp)
            in_span = (rid is not None) and (span_lo <= rid <= span_hi)
            for c in flow_cols:
                if pd.isna(row[c]):
                    continue
                int_q = int(round(float(row[c])))
                run_idx = int(c.split("_")[1])
                rows.append({
                    "flowpath": str_fp,
                    "flow_cfs": int_q,
                    "output_name": "wsel_%s_%d.tif" % (str_fp, int_q),
                    "source_gpkg": str_name,
                    "seg_q_hi": q_hi,
                    "run_index": run_idx,
                    "in_span": in_span,
                })
    
    df = pd.DataFrame(rows)
    
    # decision: out-of-span -> skip; else keep lowest-segment copy per output_name
    df["decision"] = "WRITE"
    df.loc[~df["in_span"], "decision"] = "SKIP_out_of_span"
    
    writeable = df[df["decision"] == "WRITE"].sort_values("seg_q_hi")
    seen = set()
    for idx, r in writeable.iterrows():
        if r["output_name"] in seen:
            df.at[idx, "decision"] = "SKIP_duplicate"
        else:
            seen.add(r["output_name"])
    
    str_out = os.path.join(str_gpkg_folder, "raster_write_plan.csv")
    df.sort_values(["flowpath", "flow_cfs", "seg_q_hi"]).to_csv(str_out, index=False)
    
    # summary
    if b_print_output:
        print("\n  === planned WRITE count per reach ===")
        w = df[df["decision"] == "WRITE"]
        for fp in sorted(df["flowpath"].unique()):
            flows = sorted(w[w["flowpath"] == fp]["flow_cfs"].tolist())
            print("    %s: %2d files  max=%d  -> %s" %
                  (fp, len(flows), max(flows) if flows else 0, flows))
        
        print("\n  === totals: WRITE=%d  SKIP_out_of_span=%d  SKIP_duplicate=%d" % (
            (df["decision"] == "WRITE").sum(),
            (df["decision"] == "SKIP_out_of_span").sum(),
            (df["decision"] == "SKIP_duplicate").sum()))
        
        print("+-----------------------------------------------------------------+")
        print("plan written:", str_out)
        print("+-----------------------------------------------------------------+")
# ---------------------------------------------


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if __name__ == '__main__':

    flt_start_run = time.time()
    
    parser = argparse.ArgumentParser(description='======= DETERMINE FIMS PER STREAM REACH =======')
    
    parser.add_argument('-i',
                        dest = "str_input_dir",
                        help=r'REQUIRED: directory of hydraulic result geopackages: E:\ras2fim_2d_output_20260906\04_wsel_and_hr',
                        required=True,
                        metavar='DIR',
                        type=str)
    

    parser.add_argument('-r',
                        dest = "b_print_output",
                        help=r'OPTIONAL: Print output messages Default: True',
                        required=False,
                        default=True,
                        metavar='T/F',
                        type=fn_str_to_bool)
    

    args = vars(parser.parse_args())
    
    str_input_dir = args['str_input_dir']
    b_print_output = args['b_print_output']
    
    fn_determine_fims_per_reach(str_input_dir,
                                b_print_output)

    flt_end_run = time.time()
    flt_time_pass = (flt_end_run - flt_start_run) // 1
    time_pass = datetime.timedelta(seconds=flt_time_pass)
    
    print('Compute Time: ' + str(time_pass))
 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~