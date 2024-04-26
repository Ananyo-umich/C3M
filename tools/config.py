# Configure file to parse the arguments for C3M simulation
# Following are the list of arguments
# Input file
# Output file
# log file
# Start time, end time, time steps

import argparse, glob

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("--i", default="", help="input file")

parser.add_argument("--o", default="out.nc", help="output file")

parser.add_argument("--ex", default="", help="output file")

parser.add_argument("--tstart", default="0", help="start time")

parser.add_argument("--tend", default="10", help="end time")

parser.add_argument("--tsteps", default="10", help="Time steps")
args = vars(parser.parse_args())
