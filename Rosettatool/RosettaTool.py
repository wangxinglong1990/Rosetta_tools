import os
def select_site(mut_sites,pdb_name,cpu_number):
    mutlist = []
    filepath = os.getcwd()
    if not os.path.exists('run'):
        os.mkdir('run')
    os.system('cp 1.pdb ./run/1.pdb')
    residues_list = ['A','G','D','N','F','P','Q','C','E','H','I','K','L','M','R','S','T','V','W','Y']
    with open ('%s/run/cartddg_flag'%filepath,'w') as ddg_flag:
        ddg_flag.write ('-ddg:iterations 3\n-ddg::cartesian\n-ddg::dump_pdbs False\n-bbnbrs 1\n-fa_max_dis 9.0\n-score::weights ref2015_cart\n-relax::cartesian\n-relax:min_type lbfgs_armijo_nonmonotone\n-ex1\n-ex2\n-use_input_sc\n-flip_HNQ\n-optimization::default_max_cycles 200\n-crystal_refine\n-score:extra_improper_file LIG.tors')
        ddg_flag.close()
    with open (r'seq.txt') as seq_file:
        seqlist=seq_file.readlines()
        for i in seqlist:
            mut_sites=mut_sites.split()
            for k in mut_sites:
                k = k.replace(',',' ')
                k = k.split()
                for l in range (int(k[0]),int(k[1])+1):
                    for residue in residues_list:
                        mut = str(i[int(l)-1])+str(int(l))+str(residue)
                        mutlist.append(mut)
                        f = open ('%s/run/%s.mutfile'%(filepath,mut),'w')
                        f.write('total 1\n1\n%s'%(mut))
                        f.close()
    m = len(mutlist)//int(cpu_number)
    f=open('run.sh','w')
    f.write('cd %s\n'%filepath)
    f.close()
    for j in range(int(cpu_number)):
        f1=open('%s/run%s.sh'%(filepath,j),'w')
        f1.close()
        f2=open('run.sh','a+')
        f2.write('nohup bash run%s.sh &\n'%j)
        f2.close()
        if len(mutlist)-(j+1)*m+1 > m:
            for k in range (j*m+1,(j+1)*m+1):
                f3=open ('%s/run%s.sh'%(filepath,j),'a+')
                f3.write('\ncd %s/run\ncartesian_ddg.mpi.linuxgccrelease -s %s.pdb @cartddg_flag -ddg:mut_file %s.mutfile'%(filepath,pdb_name,mutlist[(k-1)]))
                f3.close()
        else:
            for k in range (j*m+1,len(mutlist)+1):
                f3=open ('%s/run%s.sh'%(filepath,j),'a+')
                f3.write('\ncd %s/run\ncartesian_ddg.mpi.linuxgccrelease -s %s.pdb @cartddg_flag -ddg:mut_file %s.mutfile'%(filepath,pdb_name,mutlist[(k-1)]))
                f3.close()
                
def Pscan(mut_sites,pdb_name,cpu_number):
    mutlist = []
    residue = 'P'
    filepath = os.getcwd()
    if not os.path.exists('run'):
        os.mkdir('run')
    os.system('cp 1.pdb ./run/1.pdb')
    residues_list = ['A','G','D','N','F','P','Q','C','E','H','I','K','L','M','R','S','T','V','W','Y']
    with open ('%s/run/cartddg_flag'%filepath,'w') as ddg_flag:
        ddg_flag.write ('-ddg:iterations 3\n-ddg::cartesian\n-ddg::dump_pdbs False\n-bbnbrs 1\n-fa_max_dis 9.0\n-score::weights ref2015_cart\n-relax::cartesian\n-relax:min_type lbfgs_armijo_nonmonotone\n-ex1\n-ex2\n-use_input_sc\n-flip_HNQ\n-optimization::default_max_cycles 200\n-crystal_refine\n-score:extra_improper_file LIG.tors')
        ddg_flag.close()
    with open (r'seq.txt') as seq_file:
        seqlist=seq_file.readlines()
        for i in seqlist:
            mut_sites=mut_sites.split()
            for k in mut_sites:
                k = k.replace(',',' ')
                k = k.split()
                for l in range (int(k[0]),int(k[1])+1):
                    mut = str(i[int(l)-1])+str(int(l))+str(residue)
                    mutlist.append(mut)
                    f = open ('%s/run/%s.mutfile'%(filepath,mut),'w')
                    f.write('total 1\n1\n%s'%(mut))
                    f.close()
    m = len(mutlist)//int(cpu_number)
    f=open('run.sh','w')
    f.write('cd %s\n'%filepath)
    f.close()
    for j in range(int(cpu_number)):
        f1=open('%s/run%s.sh'%(filepath,j),'w')
        f1.close()
        f2=open('run.sh','a+')
        f2.write('nohup bash run%s.sh &\n'%j)
        f2.close()
        if len(mutlist)-(j+1)*m+1 > m:
            for k in range (j*m+1,(j+1)*m+1):
                f3=open ('%s/run%s.sh'%(filepath,j),'a+')
                f3.write('\ncd %s/run\ncartesian_ddg.mpi.linuxgccrelease -s %s.pdb @cartddg_flag -ddg:mut_file %s.mutfile'%(filepath,pdb_name,mutlist[(k-1)]))
                f3.close()
        else:
            for k in range (j*m+1,len(mutlist)+1):
                f3=open ('%s/run%s.sh'%(filepath,j),'a+')
                f3.write('\ncd %s/run\ncartesian_ddg.mpi.linuxgccrelease -s %s.pdb @cartddg_flag -ddg:mut_file %s.mutfile'%(filepath,pdb_name,mutlist[(k-1)]))
                f3.close()
                
def seq():
    f = open('prep.sh','w')
    f.write('source activate python27\npython $Rose\get_fasta_from_pdb.py 1.pdb A 1.fa\nsource activate')
    f.close()
    os.system('bash prep.sh')
    f = open ('1.fa','r')
    a = f.readlines()
    f.close()
    c = 0
    for i in a:
        c += 1
        if c == 2:
            f = open ('seq.txt','w')
            f.write(i)
            f.close()

def SC_pos(charge):
    f = open('SCrun.sh','w')
    f.write('supercharge.mpi.linuxgccrelease '
            '-s 1.pdb -use_input_sc -ignore_unrecognized_res -jd2:no_output -dont_mutate_glyprocys true '
            '-dont_mutate_correct_charge true -dont_mutate_hbonded_sidechains true -include_asp '
            '-include_glu -refweight_asp -0.6 -refweight_glu -0.8 -surface_residue_cutoff 16 '
            '-target_net_charge_active -target_net_charge -%s -nstruct 100 > log'%charge)
    f.close()

def SC_neg(charge):
    f = open('SCrun.sh','w')
    f.write('supercharge.mpi.linuxgccrelease '
            '-s 1.pdb -use_input_sc -ignore_unrecognized_res -jd2:no_output -dont_mutate_glyprocys true '
            '-dont_mutate_correct_charge true -dont_mutate_hbonded_sidechains true -include_arg '
            '-include_lys -refweight_arg -1.98 -refweight_lys -1.65 -surface_residue_cutoff 16 '
            '-target_net_charge_active -target_net_charge -%s -nstruct 100 > log'%charge)
    f.close()
    
def relax1(cores,nstruct):
    f = open('relax.sh','w')
    f.write('mpirun -np %s relax.mpi.linuxgccrelease -s 1.pdb -relax:constrain_relax_to_start_coords '
            '-ramp_constraints false -relax:coord_constrain_sidechains '
            '-nstruct %s -ex1 -ex2 -use_input_sc -flip_HNQ -no_optH false'%(cores,nstruct))
    f.close()

def relax2(cores,nstruct):
    f = open('relax.sh','w')
    f.write('mpirun -np %s relax.mpi.linuxgccrelease -s 1.pdb '
            '-nstruct %s -ex1 -ex2 -use_input_sc -flip_HNQ -no_optH false '
            '-relax:cartesian -score:weights ref2015_cart -crystal_refine'%(cores,nstruct))
    f.close()
    
if __name__ == '__main__':
    seq()
    select_site(mut_sites='', pdb_name='1', cpu_number='8')
    Pscan(mut_sites='', pdb_name='1', cpu_number='8')
    ##eg: select_site(mut_sites='45,48 77,79 116,119 124,126 160,161 163,164 264,269', pdb_name='1', cpu_number='16')
    SC_pos(charge='')
    SC_neg(charge='')
    relax1(cores='8',nstruct=50)
    relax2(cores='8',nstruct=50)
