import pandas as pd

def filter_barcodes_in_df(df, min_16s_reads: int = 5, max_contam: float = 0.1, min_barcodes: int = 10):
    """
    Input: pandas dataframe indexed by barcode, and with the following columns: 
    predicted_taxonomy, confidence, contamination, total 16s reads, technical noise count, ARGs
    """
    original_num_barcodes = len(df)
    print("Pre-filtering # of barcodes:", original_num_barcodes)    


    # Stage 1 Filtering

    print(f"Stage 1: Filter out low-signal records: specifically, barcodes with low 16s reads (<={min_16s_reads}),")
    print(f" high contamination (>={max_contam}), or taxonomy only classified to bacteria-level...")

    for barcode in df.index:
        
        predicted_taxonomy = df.loc[barcode].to_list()[0]
        confidence, contamination, total_16s_reads = map(float, df.loc[barcode].to_list()[1:4])

        # filter out barcode rows
        if any([total_16s_reads <= min_16s_reads, contamination >= max_contam, \
        "P - None | C - None | O - None | F - None | G - None | S - None" in predicted_taxonomy]):            
            df.drop(index=barcode, inplace=True)

        # Normalize Klebsiella assignments
        # 2026-08-10: Only normalize rows that survived Stage 1 filtering.
        # Reason: assigning with .loc after drop() can recreate a filtered barcode.
        if barcode in df.index and "G - Klebsiella" in predicted_taxonomy:
            idx = predicted_taxonomy.find("G - Klebsiella")
            df.loc[barcode, "Predicted taxonomy"] = predicted_taxonomy[:idx] + "G - Klebsiella pneumoniae complex | S - None"


    stage_1_barcodes = len(df)
    print("  # of barcodes filtered out:", (original_num_barcodes - stage_1_barcodes))
    print("  # of barcodes remaining:", stage_1_barcodes)



    # Stage 2 Filtering

    # 2026-08-10: Stop cleanly when Stage 1 removes every barcode.
    # Reason: df.iloc[0] below raises IndexError for an empty filtered table.
    if df.empty:
        print("Barcode filtering complete - final # of barcodes: 0")
        return

    print(f"Stage 2: Filter out low-signal cells: specifically, taxonomic classifications that have <{min_barcodes} barcodes...")

    df.sort_values(by = ["Predicted taxonomy"], inplace = True)
    bcs_with_this_tax = []
    current_tax = df.iloc[0]["Predicted taxonomy"]

    for barcode in df.index:
        tax = df.loc[barcode, "Predicted taxonomy"]
        if tax == current_tax:
            if len(bcs_with_this_tax) < min_barcodes:
                bcs_with_this_tax.append(barcode)
        else:
            if len(bcs_with_this_tax) < min_barcodes:
                for bc in bcs_with_this_tax:
                    df.drop(index=bc, inplace=True)
            current_tax = tax
            bcs_with_this_tax = [barcode]

    # account for barcodes from the last taxonomy in the dataframe
    if len(bcs_with_this_tax) < min_barcodes:
        for bc in bcs_with_this_tax:
            df.drop(index=bc, inplace=True)

    stage_2_barcodes = len(df)
    print("  # of barcodes filtered out:", (stage_1_barcodes - stage_2_barcodes))
    print("  # of barcodes remaining:", stage_2_barcodes)

    print("Barcode filtering complete - final # of barcodes:", stage_2_barcodes)
