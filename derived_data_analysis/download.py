"""
Script to download unmerged outputs from hyperloop
"""

import argparse
import os

def download(infile_name, suffix, max_files_to_merge):
    """
    Main function for download
    It requires access to alien
    """

    files_to_merge = []
    with open(infile_name, "r") as infile:
        for idir, directory in enumerate(infile):
            if os.path.isdir(f"input_{idir:03d}"): # we clean previously downloaded files
                os.rmdir(f"input_{idir:03d}")
            os.mkdir(f"input_{idir:03d}")
            if os.path.isfile("files_to_download_part.txt"):
                os.system("rm files_to_download_part.txt")
            os.system("touch files_to_download_part.txt")
            os.system(f"alien_find alien://{directory.strip()} 0*/AO2D.root > files_to_download_part.txt")
            with open("files_to_download_part.txt", "r") as file_part:
                for ifile, line_part in enumerate(file_part):
                    os.system(f"alien_cp alien://{line_part.strip()} file:input_{idir:03d}/AO2D_{ifile:03d}.root")
                    os.system(f"alien_cp alien://{line_part.strip().replace('AO2D', 'AnalysisResults')} "
                              f"file:input_{idir:03d}/AnalysisResults_{ifile:03d}.root")
                    files_to_merge.append(f"input_{idir:03d}/AO2D_{ifile:03d}.root")

    nbatches = 0
    for ifile, file_to_merge in enumerate(files_to_merge):
        if ifile % max_files_to_merge == 0:
            nbatches += 1
            if os.path.isfile(f"files_to_merge_{nbatches-1}.txt"):
                os.system(f"rm files_to_merge_{nbatches-1}.txt")
            os.system(f"touch files_to_merge_{nbatches-1}.txt")
        os.system(f"echo {file_to_merge} >> files_to_merge_{nbatches-1}.txt")

    if nbatches > 1:
        if os.path.isfile("files_to_merge.txt"):
            os.system("rm files_to_merge.txt")
        os.system("touch files_to_merge.txt")
        for batch in range(nbatches):
            os.system(f"o2-aod-merger --input files_to_merge_{batch}.txt --output "
                      f"AO2D{suffix}_batch{batch}.root --max-size 1000000000")
            os.system(f"echo AO2D{suffix}_batch{batch}.root >> files_to_merge.txt")
        os.system("o2-aod-merger --input files_to_merge.txt --output "
                  f"AO2D{suffix}.root --max-size 1000000000")
    else:
        os.system("o2-aod-merger --input files_to_merge_0.txt --output "
                  f"AO2D{suffix}.root --max-size 1000000000")

    os.system(f"hadd -n {max_files_to_merge} AnalysisResults{suffix}.root input_*/AnalysisResults_*.root")

    os.system("rm files_to_download_part.txt")
    os.system("rm files_to_merge*.txt")
    os.system("rm *batch*.root")
    os.system("rm -r input_*")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Arguments")
    parser.add_argument("--infile", "-i", metavar="text", default="files_to_download_B0_LHC24.txt",
                        help="text file with input directories (from hyperloop)")
    parser.add_argument("--suffix", "-s", metavar="text", default="_LHC24_pass1_skimmed",
                        help="output suffix")
    parser.add_argument("--max_files_to_merge", "-m", type=int, default=200,
                        help="maximum number of files to merge")
    args = parser.parse_args()

    download(args.infile, args.suffix, args.max_files_to_merge)
