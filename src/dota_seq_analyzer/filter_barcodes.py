#!/usr/bin/env python3
import pandas as pd

def filter_barcodes_in_df(df, min_16s_reads: int = 5, max_contam: float = 0.1, min_barcodes: int = 10):
    """
    Filter barcode summaries while preserving the existing scientific rules.

    The input DataFrame is updated in place for compatibility with existing callers.
    """
    # 2026-09-04: Validate the numeric taxonomy threshold before filtering.
    # Reason: negative values have no meaningful minimum-cell interpretation.
    if min_barcodes < 0:
        raise ValueError("min_barcodes must be non-negative")

    original_num_barcodes = len(df)
    print("Pre-filtering # of barcodes:", original_num_barcodes)

    print(f"Stage 1: Filter out low-signal records: specifically, barcodes with low 16s reads (<={min_16s_reads}),")
    print(f" high contamination (>={max_contam}), or taxonomy only classified to bacteria-level...")

    # 2026-09-04: Apply the existing Stage 1 predicates as one vectorized mask.
    # Reason: row-wise .loc access and repeated drop() made filtering scale poorly.
    taxonomy = df["Predicted taxonomy"].astype("string")
    total_16s_reads = pd.to_numeric(df["Total # of 16s reads"])
    contamination = pd.to_numeric(df["Contamination"])
    unclassified_taxonomy = "P - None | C - None | O - None | F - None | G - None | S - None"
    stage1_mask = (
        (total_16s_reads > min_16s_reads)
        & (contamination < max_contam)
        & ~taxonomy.str.contains(unclassified_taxonomy, regex=False, na=False)
    )
    filtered_df = df.loc[stage1_mask].copy()

    # 2026-09-04: Preserve the original first-occurrence Klebsiella rewrite.
    # Reason: the old code replaced everything from the first matching genus onward,
    # but only for rows surviving Stage 1.
    klebsiella_mask = filtered_df["Predicted taxonomy"].astype("string").str.contains(
        "G - Klebsiella", regex=False, na=False
    )
    filtered_taxonomy = filtered_df.loc[klebsiella_mask, "Predicted taxonomy"].astype("string")
    filtered_df.loc[klebsiella_mask, "Predicted taxonomy"] = (
        filtered_taxonomy.str.split("G - Klebsiella", n=1).str[0]
        + "G - Klebsiella pneumoniae complex | S - None"
    )

    # 2026-09-04: Return the filtered frame instead of clearing/rebuilding the input.
    # Reason: callers now assign the result, avoiding dtype/index reconstruction risks.
    df = filtered_df
    stage_1_barcodes = len(df)
    print("  # of barcodes filtered out:", (original_num_barcodes - stage_1_barcodes))
    print("  # of barcodes remaining:", stage_1_barcodes)

    # 2026-09-04: Preserve the existing early exit when Stage 2 is disabled.
    # Reason: min_barcodes=0 means Stage 1 only, with no taxonomy sorting/filtering.
    if df.empty:
        print("Barcode filtering complete - final # of barcodes: 0")
        return df
    if min_barcodes == 0:
        print("Stage 2: Taxonomic minimum-cell filtering disabled.")
        print("Barcode filtering complete - final # of barcodes:", stage_1_barcodes)
        return df

    print(f"Stage 2: Filter out low-signal cells: specifically, taxonomic classifications that have <{min_barcodes} barcodes...")

    # 2026-09-04: Preserve the old taxonomy sort, then filter complete taxon groups.
    # Reason: the old implementation sorted before its group-wise deletion loop.
    filtered_df = df.sort_values(by=["Predicted taxonomy"]).copy()
    taxon_counts = filtered_df.groupby("Predicted taxonomy")["Predicted taxonomy"].transform("size")
    df = filtered_df.loc[taxon_counts >= min_barcodes].copy()

    stage_2_barcodes = len(df)
    print("  # of barcodes filtered out:", (stage_1_barcodes - stage_2_barcodes))
    print("  # of barcodes remaining:", stage_2_barcodes)
    print("Barcode filtering complete - final # of barcodes:", stage_2_barcodes)
    return df
