import argparse
import sys

parser = argparse.ArgumentParser()
parser.add_argument('--conf', type=str, default='conf.ini', help='Path to config file')
    
# 2. Subparsers setup
subparsers = parser.add_subparsers(dest="command", required=False)

# 'init' subcommand
subparsers.add_parser("init", help="Initialize a new project")

# 'rebin' subcommand
rebin_parser = subparsers.add_parser("rebin", help="Rebin measurement files")
rebin_parser.add_argument("measurement_file", type=str, help="Path to the measurement file")
rebin_parser.add_argument("bin_number", type=int, help="The rebinning number")

args = parser.parse_args()

if args.command and '--conf' in sys.argv:
    parser.error("Argument '--conf' cannot be used with subcommands (init, rebin).")
# return parser.parse_args()

print(args.command)

