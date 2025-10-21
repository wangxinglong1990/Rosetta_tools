import RosettaTool
import argparse

RosettaTool.seq()
parser = argparse.ArgumentParser()
parser.add_argument('-sites',dest='mut_sites',type=str)
parser.add_argument('-cores',dest='cpu_number',type=int)
args = parser.parse_args()
RosettaTool.select_site(mut_sites='%s'%args.mut_sites, pdb_name='1', cpu_number='%s'%args.cpu_number)
