"""Write _one_<section>.tex: report.tex with only the named section input, so
each section compiles on its own and two sections can compile at once.
Usage: python one_section.py 03_genetic"""
import re
import sys

name = sys.argv[1]
text = open("report.tex").read()
text = re.sub(r"\\input\{sections/[^}]*\}\n", "", text)
text = text.replace("\\maketitle", "\\maketitle\n\\input{sections/%s}" % name, 1)
open("_one_%s.tex" % name, "w").write(text)
