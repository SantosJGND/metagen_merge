import os

import pandas as pd

from modules.ncbi_connect import entrez_fetch_taxid_from_org_description_curate


def process_televir(televir_reports):
    telebac_found = televir_reports[
        [
            "Sample",
            "Taxid",
            "Description",
            "accID",
            "Cov (%)",
            "Depth",
            "DepthC",
            "Mapped reads",
            "Windows Covered",
            "Warning",
            "Control",
        ]
    ].drop_duplicates()
    telebac_found["Windows Covered"] = telebac_found["Windows Covered"].apply(
        lambda x: x.replace("/", "-")
    )
    telebac_found["Taxid"] = telebac_found["Taxid"].astype(str)

    return telebac_found


def read_televir(televir_filepath):
    """
    read tsv file from televir export
    """
    reports = pd.read_csv(televir_filepath, sep="\t")
    return reports


def read_panel(report, panel="Microorganisms"):
    """
    read excel extract spreadsheet name Microorganisms

    """
    panel = pd.read_excel(report, sheet_name=panel)
    return panel


def get_illumina_found(rpip_panel, upip_panel):
    """
    merge rpip and upip panels, retrieve taxids
    """
    illumina_found_full = pd.concat([rpip_panel, upip_panel], axis=0, ignore_index=True)

    ## discard description= "None"
    illumina_found = illumina_found_full[
        illumina_found_full["Microorganism Name"] != "None"
    ].drop_duplicates(subset=["Microorganism Name"])
    illumina_found = illumina_found[
        ["Accession", "Microorganism Name", "Coverage", "ANI", "Median Depth", "RPKM"]
    ]

    illumina_found["Taxid"] = illumina_found["Microorganism Name"].apply(
        entrez_fetch_taxid_from_org_description_curate
    )

    illumina_found.rename(
        columns={"Accession": "Sample", "Microorganism Name": "Description"},
        inplace=True,
    )
    illumina_found = illumina_found.explode("Taxid").drop_duplicates()

    return illumina_found


def merge_panels(illumina_found, telebac_found):
    """
    merge illumina_found and televir_reports
    """

    all_samples = illumina_found.Sample.unique()

    final_set = pd.DataFrame()

    for sample in all_samples:
        illumina_sample = illumina_found[
            illumina_found["Sample"] == sample
        ].reset_index(drop=True)
        telebac_sample = telebac_found[telebac_found["Sample"] == sample].reset_index(
            drop=True
        )

        sample_taxids = (
            pd.concat([illumina_sample["Taxid"], telebac_sample["Taxid"]], axis=0)
            .drop_duplicates()
            .reset_index(drop=True)
        )
        new_set = pd.merge(
            illumina_sample, telebac_sample, on=["Sample", "Taxid"], how="outer"
        )

        new_set = new_set.drop_duplicates(
            subset=["Sample", "Description_x", "accID"]
        ).reset_index(drop=True)

        final_set = pd.concat([final_set, new_set], axis=0, ignore_index=True)

    return final_set
