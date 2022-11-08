"""

"""

samtools_test = "/Users/winx/Documents/Bioinformatics/stecode/tests/data/test1_eaesubtype.tab"
stecvir_test = "/Users/winx/Documents/Bioinformatics/stecode/tests/data/test1_stecvir.tab"
recAstxeaehly_test = "/Users/winx/Documents/Bioinformatics/stecode/tests/data/test1_2recAstxeaehly.txt"

import parsers.samtools_parse as sp
import parsers.recAstxeaehly_parse as rp
import parsers.stecvir_parse as vp

def merge_all_NNs(stfile, recAfile, virfile, longread):
    #NN1 = sp.samtools_input(stfile)
    NN2 = rp.recA_input(recAfile, longread)
    #NN3456 = vp.vfdb_input(virfile)

    return NN2
print(merge_all_NNs(samtools_test, recAstxeaehly_test, stecvir_test, False))
