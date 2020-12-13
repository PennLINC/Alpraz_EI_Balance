# A python script for printing the inputs and config of an analysis on flywheel
# Arguments:
#     The Session ID (e.g. 5c8ab80bf546b6002f6d78e9)
#     The analysis label (must not be exact but will fail if not close enough) (e.g. fmriprep-hpc)

import flywheel
import sys

# get arguments from command line
if len(sys.argv) is not 3:
    print("Incorrect arguments given! (see file header)")
    sys.exit(1)
args = sys.argv

# get flywheel client
client = flywheel.Client()

# get the target session
session = client.get(args[1])

if not session:
    print("session not found:", args[1])

# get the target analysis
target_analysis = [x for x in session.analyses if args[2] in x.label]

if not target_analysis:
    print("analysis not found:", args[2])
    sys.exit(1)

print('analysis inputs:', [x.name for x in target_analysis[0].inputs])
print('analysis config:', [x.job.config['config'] for x in target_analysis][0])
print("Done!")

