from distutils.command.clean import clean
from tqdm import tqdm
import argparse
import sys
import os
import subprocess

tmp_dir="~/s/tmp/"
code_dir="~/s/proj4/scr/"
file_map=[]

#list of all memttime
class SizeCollect:
    size_a=0
    size_b=0

class MemTime:
    mem_kb=0
    time_s=0

def str_mem_time(base, tool):
    sh("mkdir -p mem_time")
    return "/usr/bin/time  -f \"%M\t%e\" --output-file=mem_time/{base}_{tool}.txt ".format(base, tool)   

def out_sh(cmd):
    return subprocess.check_output(cmd)
    
def sh(cmd):
    subprocess.call(cmd, shell=True)  
    
def generate_kmer_sets_in_spss(list_fq, kmer_size, ab, run=True, cleanup=False):
    #fout = open("demofile2.txt", "w")
    if(run==True):
        sh("rm -rf list_essd")
        count=0
        with open(list_fq, 'r') as list_file:
            for fq in tqdm(list_file, desc="ESS-Compress running on individual files"):
                filename_fq = fq.rstrip()
                file_map.append(filename_fq)
                basename_fq = os.path.basename(filename_fq)
                subprocess.call("mkdir -p mem_time", shell=True)
                cmd = str_mem_time(basename_fq, "essc_indi") + "essCompress -i {} -k {} -a {}".format(fq.rstrip(), kmer_size, ab)
                subprocess.call(cmd, shell=True)  
                
                sh("mv {}.essc s{}.essc".format(filename_fq, count))
                sh(str_mem_time("s{}.essc", "essd_indi".format(count)) + "essDecompress s{}.essc".format(count))
                count += 1
        
        sh("mkdir -p sets ; mv *.essd sets/")        
        sh("mkdir -p tmp_essc ; mv *.essc tmp_essc/")
        
        print("Output in "+ os.getcwd()+"/sets")
        
        sh("ls {}/sets/*.essd > list_essd".format( os.getcwd()))
    else:
        print("Skipping method "+ __name__)
    
    if(cleanup):
        sh("rm -rf tmp_essc")        
        
    return os.getcwd()+"/sets"


def spss_union_kmerlist(list_essd, kmer_size, run=True, cleanup=False):
#requires list_essd
    if(run==True):
        basename_list = os.path.basename(list_essd)
        sh(str_mem_time(basename_list, "essc_union") + "essCompress -i {} -k {}".format(list_essd, kmer_size)) 
        sh(str_mem_time(basename_list, "essd_union") + "essDecompress {}.essc".format(list_essd))         
        sh("mv {l}.essc union.essc; mv {l}.essd union.essd".format(l=list_essd))
        print("Output in union.essc and union.essd")
        sh(str_mem_time(basename_list, "kmerlist") + "{}/kmerlist -k {} -u {}".format(code_dir, kmer_size, "union.essd"))
        print("Output in ess_kmer_id.txt")
    else:
        print("Skipping method "+ __name__)
    return

def rename_input_and_keep_mapping(list_input, ext):
    #accepted extension: "fa, fq, fq.gz, fa.gz, fasta, fasta.gz"
    #outputs a file called mapping - list_filename.colors 
    #for all file in list
    count=0
    with open(list_input, 'r') as list_file:
        for fq in tqdm(list_file, desc="Renaming input files"):
            filename_f= fq.rstrip()
            sh("mv {} r{}.{}".format(filename_f, count, ext))
            count += 1
                        
def essd_to_kmcascii(essd_file, k, ab):
    #expected s0.essd, s1.essd, s2.essd... 
    bz=out_sh("basename {} .essd".format(essd_file))
    sh("mkdir -p kmc_tmp_dir/")
    sh(str_mem_time(bz, "kmc_only") +  "kmc -k{k} -m{m} -ci{ab} -fa {fasta} {bz}.kmc kmc_tmp_dir/".format(k,24,ab,essd_file, bz))  
    sh(str_mem_time(bz, "kmc_only") +  "kmc -k{k} -m{m} -ci{ab} -fa {fasta} {bz}.kmc kmc_tmp_dir/".format(k,24,ab,essd_file, bz))  
    sh("kmc_dump -ci{ab} {bz} {bz}.kmers".format(ab=ab, bz=bz))
    return

def kmer_sets_to_colormatrix(list_essd, kmer_size, run=True, cleanup=False):
    if(run==True):
        basename_list = os.path.basename(list_essd) #s0.essd
        with open(list_essd, 'r') as list_file:
            for fq in tqdm(list_file, desc="Renaming input files"):
                filename_f= fq.rstrip()
                essd_to_kmcascii(filename_f,kmer_size, 1)
                sh("mv {} r{}.{}".format(filename_f, count, ext))
                count += 1
    else:
        print("Skipping method "+ __name__)
    return
    return

def spss_to_kmc():
    return

def kmc_to_joincounts():
    return

def join_jc_ess_kmer_id():
    return

def bit_compression():
    return

def get_stats():
    return

def main(args):
    os.chdir(args.outdir)
    print("Current working directory: {0}.".format(os.getcwd()))
    generate_kmer_sets_in_spss(args.list, args.kmersize, args.abundance, run=True, cleanup=False)
    spss_union_kmerlist("list_essd", args.kmersize, run=True, cleanup=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Mega pipeline for sets of k-mer set compression.')
    parser.add_argument('-l', '--list', help='input list file of fastq.gz (include full pathname)')
    parser.add_argument('-k', '--kmersize', help='k value')
    parser.add_argument('-a', '--abundance', help='ab value')
    # parser.add_argument('-o', '--outdir', help='output working directory')
    args = parser.parse_args()
    main(args)
    
