#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  untitled.py
#
#  Copyright 2018 Unknown <mknigge@mknigge>
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#
#

# IMPORTS
import sys
import os
import re
import uuid
import time
import subprocess


# GLOBALS
file = open("cell.counts.proportion.txt", "r")
count = 0 # global counter
job_limit = 100


# FUNCTIONS
def global_list(phenotype):
    global phenotype_list
    phenotype_list = phenotype


def global_count():
    global count
    count += 1


def make_bash(database, phenotype):

    print("Executing job for {} cohort, phenotype {}".format(database, phenotype))

    bash = """#!/bin/bash
#SBATCH --job-name={0}.{1}
#SBATCH --output=log/{0}.{1}.out
#SBATCH --err=log/{0}.{1}.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --get-user-env=L
#SBATCH --mem=120gb

module load R

Rscript glmnet.perform.elastic.net.R -p {0} -d {1}""".format(phenotype, database)

    name = "{0}.{1}".format(database, phenotype)
    output = open("{}.sh".format(name), "w")
    output.write(bash)
    output.close()
    os.system("sbatch {}.sh".format(name))
    os.remove("{}.sh".format(name))

def main(args):

    phenotype_list = list()

    print("Building job list.")

    for line in file:
        # get phenotype and database
        tmp = line.strip().split()
        # prepare
        phenotype = "{0}|{1}".format(tmp[0], tmp[1])
        # bind to list
        phenotype_list.append(phenotype)
    global_list(phenotype_list)

    print("Executing jobs, as starting point.")

# TODO code ends here

    for x in range(0, job_limit):
        # database
        database = phenotype_list[x].split("|")[1]
        # phenotype
        phenotype  = phenotype_list[x].split("|")[0]
        # create job
        make_bash(database = database, phenotype=phenotype)
        # keep track
        global_count()

    print("Done, start listening. When jobs are finished, execute difference. Always maintain {} active jobs.".format(job_limit))

    while True:
        # get amount of running jobs
        jobs = int(subprocess.check_output("squeue -u umcg-mknigge | wc -l", shell=True).strip())-2
        # if jobs are done, the amount will not be the same
        if jobs != job_limit:
            print("{} jobs are running, executing {} jobs to maintain {} active threshold.".format(jobs, job_limit-jobs, job_limit))
            # create new jobs for the difference
            for i in range(0, job_limit-jobs):
                # database
                database = phenotype_list[count].split("|")[1]
                # phenotype
                phenotype  = phenotype_list[count].split("|")[0]
                # create job
                make_bash(database = database, phenotype=phenotype)
                # keep track
                global_count()



    return 0

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
