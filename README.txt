
No install required. Can be executed directly predecessing 'python3  ' or './' in the call.


Requirements: 
- Python3 (or higher) 
- numpy


usage: pwanalygen.py [-h] [-w WORKDIR] [-i MODE_INTERACTIVE] [--pval PVAL]
                     [--pval-cpatt PVAL_CPATT] [--pval-let PVAL_LET]
                     [--pval-num PVAL_NUM] [--pval-spec PVAL_SPEC]
                     inpwfile outpwfile

pwanalygen.py is a pw word list sec tool, that includes a sophisticated, data-
science based word list analyzer and a compatible word list generator, which
also builds the new frequency efficient word list following the analyzed pw
construction patterns as a proof-of-concept implementation.

positional arguments:
  inpwfile
  outpwfile

optional arguments:
  -h, --help            show this help message and exit
  -w WORKDIR, --workdir WORKDIR
  -i MODE_INTERACTIVE, --interactive-mode MODE_INTERACTIVE
  --pval PVAL
  --pval-cpatt PVAL_CPATT
  --pval-let PVAL_LET
  --pval-num PVAL_NUM
  --pval-spec PVAL_SPEC

The program can create A LOT OF NEW PWs based on the analyzed pw construction
patterns in the original <inpwfile>. It can 'pump up' the original <inpwfile>
by magnitudes in size, from e.g. 50k pws to e.g. 1M pws, or even more. The
data-science and word list based pw break approach is performed in three main
steps. In the first step the original pws are split into substrings of lettter
symbols, substrings of number symbols and substrings of special symbols. E.g.
the pw 'love1982!' gets splitted into 'love', '1982' and '!'. Each original pw
also gets transformed into a pw construction pattern, in the example to
'AAAA1111$', which subsequently gets further aggregated to the condensed pw
construction pattern 'A1$'. (To be interpreted as: A series of letters
followed by a series of numbers followed by a series of special characters.)
In the second step the relative cumulated frequencies of the let-, num-, and
special-substrings are computed; also the relative cumulated frequencies of
the condensed pw construction patterns. Their high-frequency outcomes below or
at the critical p-value get selected; the remaining lower frequency outcomes
are cut off (p-value=0.0 --> select 0% of the outcomes; p-value=1.0 --> select
100% of the outcomes). In the third step the selected pw element outcomes get
combined straight following the selected condensed pw construction patterns
(called 'cpattprod' inside the program). The final size of the generated
<outpwfile> is steered by the specified p-value. A high p-value creates a more
large pw list, a low p-value creates a (more) small one (pw list compression
functionality at very low p-values). The more close the p-value gets to 1.0,
the more >> over-linear will the generated <outpwfile> increase.<< As far the
<inpwfile> is not trivial simple structured, high p-vales at some point will
undenyably result in a 'never' ending pw generation job. Thus, it is highly
recommended to start with moderate p-values first (e.g. 0.5) and to switch to
interactive mode for higher p-values, in order to find individual p-values by
category, that match the trade-off between the covered pw element outcome
proportion versus the maximum acceptable generation job time/ maximum
generated file size at best. Therefor view the printouts of the relative
cumulated frequencies of every pw element category. Finally, sec officers and
sysadmins can use the generated <outpwfile> to perform a simulated pw try-out
on another unknown pw list, in order to derive a rough estimate of the own pws
at risk proportion. However, the major information gain of the tool is the
deep insight, how pws are structured and how relatively short pws can be
assembled, that none the less are pretty safe. If there is interest in this
direction, further complementary tools for -Cleansing disturbing characters
'cleansetxtfile.py' -Random sampling mega large (rockyou.txt) pw files
'samplefile.py' - Simulated pw try-out 'f1prop-in-f2.py' - can be offered.


