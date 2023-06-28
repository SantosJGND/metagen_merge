import csv
import os
import random

from pywebio import start_server
from pywebio.input import *
from pywebio.output import *

from modules.process import (
    get_illumina_found,
    merge_panels,
    process_televir,
    read_panel,
    read_televir,
)

output_dir = "assets/"


def app():
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    ## get session id
    session_id = random.randint(0, 100000)
    # session_id = session_id.split("-")[0]

    put_markdown("## metagenomics file merger")

    put_markdown("### 1. upload files")
    file1 = file_upload("upload televir reports", accept="*")
    file2 = file_upload("upload rpip report", accept="*")
    file3 = file_upload("upload upip report", accept="*")

    put_markdown("### 2. merge files")
    ## read files using pandas

    open(output_dir + f"{session_id}_televir.tsv", "wb").write(file1["content"])
    open(output_dir + f"{session_id}_rpip.xlsx", "wb").write(file2["content"])
    open(output_dir + f"{session_id}_upip.xlsx", "wb").write(file3["content"])

    ## display waiting message
    put_markdown("### processing files...")

    try:
        televir_reports = read_televir(output_dir + f"{session_id}_televir.tsv")
        rpip_panel = read_panel(
            output_dir + f"{session_id}_rpip.xlsx", panel="Microorganisms (RPIP)"
        )
        upip_panel = read_panel(
            output_dir + f"{session_id}_upip.xlsx", panel="Microorganisms (UPIP)"
        )

        illumina_found = get_illumina_found(rpip_panel, upip_panel)
        telebac_found = process_televir(televir_reports)

        merged_panel = merge_panels(illumina_found, telebac_found)
        merged_panel.to_csv(
            output_dir + f"{session_id}_merged_panel.tsv", sep="\t", index=False
        )
        ## file to bytes
        merged_panel_bytes = open(
            output_dir + f"{session_id}_merged_panel.tsv", "rb"
        ).read()

        put_markdown("### 3. download merged file")
        put_file(f"{session_id}_merged_panel.tsv", merged_panel_bytes)

    except Exception as e:
        raise (e)

    finally:
        # remove files
        os.remove(output_dir + f"{session_id}_televir.tsv")
        os.remove(output_dir + f"{session_id}_rpip.xlsx")
        os.remove(output_dir + f"{session_id}_upip.xlsx")
        os.remove(output_dir + f"{session_id}_merged_panel.tsv")


if __name__ == "__main__":
    start_server(app, host="127.0.0.1", port=8000, debug=False)
