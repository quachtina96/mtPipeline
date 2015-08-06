import vcf
import mtVariantCaller
# defines mathematical operations
# sum
path="/gpfs/home/quacht/ID18/ID18_Father/bam/myMtoolbox_out/" 
tablefile = "/gfps/home/quacht/ID18/ID18_Father/results/ID18_Father-mtDNAassembly-table.txt"
sam_handle = path+'ID18_Father.rg.ra.marked.norg.sam"'
mt_table_handle = tablefile
sample_name = "ID18_Father"

#open the marked sam file with no read group for reading
sam = open(sam_handle, 'r')
mtable = open(tablefile, 'r').readlines()
tail=5
#mut_events = mtvcf_main_analysis(mt_table, sam_file, str(sample_name), tail=5)
mtDNA = []
for i in mtable:
    print i
    mtDNA.append((i[1]).strip())
    print mtDNA
mtDNAseq = "".join(mtDNA)


