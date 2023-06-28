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
        [
            "Accession",
            "Microorganism Name",
            "Class Type",
            "Coverage",
            "ANI",
            "Median Depth",
            "RPKM",
        ]
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

        new_set = pd.merge(
            illumina_sample, telebac_sample, on=["Sample", "Taxid"], how="outer"
        )

        new_set = new_set.drop_duplicates(
            subset=["Sample", "Description_x", "accID"]
        ).reset_index(drop=True)

        def support_for_taxid(taxid: str):
            in_illumina = taxid in illumina_sample["Taxid"].values
            in_telebac = taxid in telebac_sample["Taxid"].values

            if in_illumina and in_telebac:
                return "Both"
            elif in_illumina:
                return "Panels"
            elif in_telebac:
                return "TELEVir"
            else:
                return None

        new_set["Support"] = new_set["Taxid"].apply(support_for_taxid)
        new_set = new_set[
            [
                "Sample",
                "Support",
                "Taxid",
                "Description_y",
                "accID",
                "Cov (%)",
                "Depth",
                "DepthC",
                "Mapped reads",
                "Windows Covered",
                "Warning",
                "Control",
                "Description_x",
                "Class Type",
                "Coverage",
                "ANI",
                "Median Depth",
                "RPKM",
            ]
        ]

        new_set.rename(
            columns={
                "Description_x": "Description",
                "Description_y": "TELEVir Description",
            },
            inplace=True,
        )
        # sort by Mapped reads, then windows covered float
        new_set["windows_covered_float"] = new_set["Windows Covered"].apply(
            lambda x: float(x.split("-")[0]) / float(x.split("-")[1])
        )
        new_set = new_set.sort_values(
            by=["Mapped reads", "windows_covered_float"], ascending=False
        ).reset_index(drop=True)
        # drop windows covered float
        new_set = new_set.drop(columns=["windows_covered_float"])

        ## drop duplicate taxids, keep the one with the highest Mapped reads
        new_set = new_set.drop_duplicates(subset=["Taxid"], keep="first")

        final_set = pd.concat([final_set, new_set], axis=0, ignore_index=True)

    return final_set
