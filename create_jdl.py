#!/usr/bin/env python

from os.path import isdir, join
import sys
if not isdir("jdl_creator"):
    raise Exception("Could not find directory \"jdl_creator\".")
sys.path.append("jdl_creator")
from os import listdir
from classes.JDLCreator import JDLCreator
import argparse


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("process", type=str, help="Physics process to be processed")
    parser.add_argument("output", type=str, help="Output directory for JDL file")
    return parser.parse_args()


def main(args):
    jobs = JDLCreator("docker")

    jobs.executable = "job.sh"
    jobs.wall_time = 1 * 60 * 60
    jobs.memory = 2048
    jobs.accounting_group = "cms.higgs"
    jobs.image = "stwunsch/slc6-condocker:smhtt"

    # Build list of arguments
    arguments = []
    if not isdir("data"):
        raise Exception("Could not find directory \"data\".")
    count = 0
    for f in listdir("data"):
        if args.process in f:
            for line in open(join("data", f)).readlines():
                arguments.append("{} {} {}".format(count, args.process, line).rstrip())
                count += 1
    jobs.arguments = arguments

    # The job requires lots of CPU resources
    jobs.requirements = '(Target.ProvidesCPU == True) && (Target.ProvidesEKPResources == True) && (Target.CloudSite != "blade")'
    jobs.job_folder = args.output
    jobs.WriteJDL()


if __name__ == "__main__":
    args = parse_arguments()
    main(args)
