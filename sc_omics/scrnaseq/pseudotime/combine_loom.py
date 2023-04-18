import argparse
import loompy

parser = argparse.ArgumentParser(description="Combine multiple loom files into one")
parser.add_argument("files", nargs="+", help="Input loom files to combine")
parser.add_argument("-o", "--output", required=True, help="Output file name")

args = parser.parse_args()

loompy.combine(args.files, args.output, key="Accession")