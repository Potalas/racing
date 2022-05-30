#!/opt/homebrew/bin/python3

###############################################################################
# This script is the command that is executed every run.
# Check the examples in examples/
#
# This script is run in the execution directory (execDir, --exec-dir).
#
# PARAMETERS:
# argv[1] is the candidate configuration ID
# argv[2] is the instance ID
# argv[3] is the seed
# argv[4] is the instance name
# The rest (argv[5:]) are parameters to the run
#
# RETURN VALUE:
# This script should print one numerical value: the cost that must be minimized.
# Exit with 0 if no error, with 1 in case of error
###############################################################################

import datetime
import os.path
import re
import subprocess
import sys
import pandas as pd

# algo executable path
exe = "./ga"
# config file path for the algo
config_file = "config.csv"
# pre-defined buckets
buckets = {1:[1,1],2:[2,2],3:[3,3],4:[4,4],5:[5,5],6:[6,10]}

if __name__=='__main__':
    if len(sys.argv) < 5:
        print("\nUsage: ./target-runner.py <configuration_id> <instance_id> <seed> <instance_path_name> <list of parameters>\n")
        sys.exit(1)

# Get the parameters as command line arguments.
configuration_id = sys.argv[1]
instance_id = sys.argv[2]
seed = sys.argv[3]

# read the config file and create a config dataframe
dfconfig = pd.read_csv(config_file, sep=";")

# update the seed in the config dataframe
dfconfig.loc[dfconfig["name"] == "seed", "value"] = seed

# update the instance (problem for IOH) in the config dataframe
dfconfig.loc[dfconfig["name"] == "problem", "value"] = instance_id

# convert the arguments from irace into a string
cl_probs = sys.argv[5:]
values = ""
for i in range(len(cl_probs)):
    values += re.sub(".*=", "", cl_probs[i])
    if i < len(cl_probs) - 1:
        values += ","

# convert the bucket list into a string
sizes = ""
for i in buckets.keys():
    if (sizes != ""):
        sizes += ","
    sizes += str(buckets[i][0])
    sizes += ":"
    sizes += str(buckets[i][1])

# update the seed in the config dataframe
dfconfig.loc[dfconfig["name"] == "BucketBitMutation", "argument"] = sizes + "|" + values

# rewrite the config file
dfconfig.to_csv(config_file, sep=";", index = False)

# run command
cmd = [exe, config_file]

# Define the stdout and stderr files.
out_file = "c" + str(configuration_id) + "-" + str(instance_id) + str(seed) + ".stdout"
err_file = "c" + str(configuration_id) + "-" + str(instance_id) + str(seed) + ".stderr"

def target_runner_error(msg):
    now = datetime.datetime.now()
    print(str(now) + " error: " + msg)
    sys.exit(1)

def check_executable(fpath):
    fpath = os.path.expanduser(fpath)
    if not os.path.isfile(fpath):
        target_runner_error(str(fpath) + " not found")
    if not os.access(fpath, os.X_OK):
        target_runner_error(str(fpath) + " is not executable")

# This is an example of reading a number from the output.
def parse_output(out):
    match = re.search(r'Best ([-+0-9.eE]+)', out.strip())
    if match:
        return match.group(1);
    else:
        return "No match"
        
check_executable (exe)

outf = open(out_file, "w")
errf = open(err_file, "w")
return_code = subprocess.call(cmd, stdout = outf, stderr = errf)

outf.close()
errf.close()

if return_code != 0:
    target_runner_error("command returned code " + str(return_code))

if not os.path.isfile(out_file):
    target_runner_error("output file " + out_file  + " not found.")

cost = parse_output (open(out_file).read())
print(open(out_file).read().strip())

os.remove(out_file)
os.remove(err_file)
sys.exit(0)

