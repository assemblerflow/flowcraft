#!/usr/bin/env python3
from os.path import basename
import argparse
import subprocess

def main():
    def msg(name=None):
        return '''sample_fastq.py -P1 [R1.fastq.gz] -P2 [R2.fastq.gz] -gs [Genome size (MB)] -tdep [Target Depth]'''


    parser = argparse.ArgumentParser(description="This script uses seqtk to sample a file down to target depth", usage=msg())
    parser.add_argument('-P1', nargs='?', type=str, help='Paired End fastq.gz 1', required=True)
    parser.add_argument('-P2', nargs='?', type=str, help='Paired End fastq.gz 2', required=True)
    parser.add_argument('-gs', nargs='?', type=float, help='genome size (MB)', required=True)
    parser.add_argument('-tdep', nargs='?', type=int, help="Target desired depth", required=True)
    args = parser.parse_args()

    GenomeSize = args.gs
    TargetDepth = args.tdep
    P1 = args.P1
    P2 = args.P2
    BN1=".".join(basename(P1).split('.')[:-2])
    BN2=".".join(basename(P2).split('.')[:-2])

    R1_fqchk = subprocess.Popen(['seqtk', 'fqchk', P1], stdout=subprocess.PIPE)
    R1_stdout, R1_stderr = R1_fqchk.communicate()
    B_P1=int(R1_stdout.splitlines()[2].split()[1])
    print("Bases P1:"+str(B_P1))

    R2_fqchk = subprocess.Popen(['seqtk', 'fqchk', P2], stdout=subprocess.PIPE)
    R2_stdout, R2_stderr = R2_fqchk.communicate()
    B_P2= int(R2_stdout.splitlines()[2].split()[1])
    print("Bases P2:"+str(B_P2))
    print("")

    EstCov = (B_P1 + B_P2)/ (GenomeSize * 1E6)
    print ("Estimated coverage: "+str(EstCov))
    Ratio = TargetDepth/EstCov

    print("Subsample target ratio:"+str(Ratio))
    if Ratio < 1:
        #print ("Writing R1.fq.gz")
        ps = subprocess.Popen(('seqtk', 'sample','-s100', P1, str(Ratio)), stdout=subprocess.PIPE)
        with open('{}_ss.fq.gz'.format(BN1),'w') as outfile:
            output = subprocess.Popen(('pigz', '--fast', '-c'), stdin=ps.stdout, stdout=outfile )
        ps.wait()

        #print ("Writing R2.fq.gz")
        ps = subprocess.Popen(('seqtk', 'sample', '-s100', P2, str(Ratio)), stdout=subprocess.PIPE)
        with open('{}_ss.fq.gz'.format(BN2),'w') as outfile:
            output = subprocess.Popen(('pigz', '--fast', '-c'), stdin=ps.stdout, stdout=outfile )
        ps.wait()
        #print("All done. Have a nice day!")
        # Get real path of the symlink
        for fq in [P1,P2]:
            rp = os.path.realpath(fq)
            print("removing temporary fastq file path: {}".format(rp))
            # remove only when the file is in the work directory
            if re.match(".*/work/.{2}/.{30}/.*", rp):
                os.remove(rp)
if __name__ == "__main__":
        main()
