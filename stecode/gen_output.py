"""

"""

samtools_test = "/Users/winx/Documents/Bioinformatics/stecode/tests/data/test1_eaesubtype.tab"
stecvir_test = "/Users/winx/Documents/Bioinformatics/stecode/tests/data/test1_stecvir.tab"
recAstxeaehly_test = "/Users/winx/Documents/Bioinformatics/stecode/tests/data/test1_2recAstxeaehly.txt"

import parsers.samtools_parse as sp

def merge_all_NNs():
    sp.samtools_input(samtools_test)