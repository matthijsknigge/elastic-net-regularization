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


def main(args):
	
	
	file = open("cytokines.txt", "r")
	for line in file:

		for i in range(1, 11):


			bash = """#!/bin/bash
#SBATCH --job-name={0}.fold{1}
#SBATCH --output=log/{0}.fold{1}.out
#SBATCH --err=log/{0}.fold{1}.err
#SBATCH --time=6:00:00
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --get-user-env=L
#SBATCH --mem=20gb

module load R

Rscript ../cytokines.R -c {0} -f {1}
			   """.format(line.strip(), i)
			bash_file = open("bash/{0}.fold{1}.sh".format(line.strip(), i), "w")
			bash_file.write(bash)
			bash_file.close()

	
	return 0

if __name__ == '__main__':
	import sys
	sys.exit(main(sys.argv))
