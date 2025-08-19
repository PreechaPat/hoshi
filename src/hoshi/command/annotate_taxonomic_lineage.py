import pandas as pd
from collections import OrderedDict
from hoshi.lib.taxdb import get_taxdb, get_filtered_lineage


def annotate_with_taxonomic_lineage(
    df,
    taxdb=None,
    taxid_col="tax_id",
    abundance_col="abundance",
    desired_ranks=None,
    skip_invalid=True,
):
    if desired_ranks is None:
        desired_ranks = [
            "species",
            "genus",
            "family",
            "order",
            "class",
            "phylum",
            "kingdom",
        ]

    if taxdb is None:
        taxdb = get_taxdb()

    lineage_records = []
    for _, row in df.iterrows():
        taxid = row[taxid_col]
        lineage_data = OrderedDict((rank, "NA") for rank in desired_ranks)

        try:
            taxid_int = int(taxid)
            lineage = get_filtered_lineage(taxid_int, taxdb, desired_ranks)
            for k, v in lineage.items():
                lineage_data[k] = v
        except (ValueError, TypeError):
            if skip_invalid:
                continue

        merged = {**row.to_dict(), **lineage_data}
        lineage_records.append(merged)

    return pd.DataFrame(lineage_records)


def main():
    import argparse
    from hoshi.lib.ingress import read_input_table

    parser = argparse.ArgumentParser(
        description="Annotate tax_id table with lineage columns."
    )
    parser.add_argument("input_file", help="Path to input table with 'tax_id' column.")
    parser.add_argument(
        "-o",
        "--output_file",
        help="Path to save annotated output (TSV). If not provided, prints to stdout.",
    )
    parser.add_argument(
        "--taxid_col",
        default="tax_id",
        help="Name of the tax ID column (default: 'tax_id').",
    )
    parser.add_argument(
        "--abundance_col",
        default="abundance",
        help="Name of the abundance column (default: 'abundance').",
    )

    args = parser.parse_args()

    df_input = read_input_table(
        args.input_file, required_columns=[args.taxid_col, args.abundance_col]
    )

    df_annotated = annotate_with_taxonomic_lineage(
        df_input,
        taxid_col=args.taxid_col,
        abundance_col=args.abundance_col,
        skip_invalid=False,
    )

    if args.output_file:
        df_annotated.to_csv(args.output_file, sep="\t", index=False)
    else:
        print(df_annotated.to_csv(sep="\t", index=False))


if __name__ == "__main__":
    main()
