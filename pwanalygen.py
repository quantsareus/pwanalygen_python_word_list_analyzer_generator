#! /bin/python3


# pwanalygen
# By quantsareus.net
version="0.12.13"
# License: GPL-V3


from argparse import ArgumentParser
import numpy as np
try:
    import os
except:
    pass


################################################################################

### Argument parsing and initialization


parser= ArgumentParser()

parser.add_argument("-w", "--workdir", dest="workdir", default="/tmp")
parser.add_argument("-i", "--interactive-mode", dest="mode_interactive", default=False)
parser.add_argument("--pval", dest="pval", default="0.50")
parser.add_argument("--pval-cpatt", dest="pval_cpatt", default="")
parser.add_argument("--pval-let", dest="pval_let", default="")
parser.add_argument("--pval-num", dest="pval_num", default="")
parser.add_argument("--pval-spec", dest="pval_spec", default="")
parser.add_argument("inpwfile", type= str)
parser.add_argument("outpwfile", type= str)

parser.description="pwanalygen.py is a pw word list sec tool, that includes a sophisticated, data-science based word list analyzer and a compatible word list generator, which also builds the new frequency efficient word list following the analyzed pw construction patterns as a proof-of-concept implementation."

parser.epilog="The program can create A LOT OF NEW PWs based on the analyzed pw construction patterns in the original <inpwfile>. It can 'pump up' the original <inpwfile> by magnitudes in size, from e.g. 50k pws to e.g. 1M pws, or even more. The data-science and word list based pw break approach is performed in three main steps. In the first step the original pws are split into substrings of lettter symbols, substrings of number symbols and substrings of special symbols. E.g. the pw 'love1982!' gets splitted into 'love', '1982' and '!'. Each original pw also gets transformed into a pw construction pattern, in the example to 'AAAA1111$', which subsequently gets further aggregated to the condensed pw construction pattern 'A1$'. (To be interpreted as: A series of letters followed by a series of numbers followed by a series of special characters.) In the second step the relative cumulated frequencies of the let-, num-, and special-substrings are computed; also the relative cumulated frequencies of the condensed pw construction patterns. Their high-frequency outcomes below or at the critical p-value get selected; the remaining lower frequency outcomes are cut off (p-value=0.0 --> select 0% of the outcomes; p-value=1.0 --> select 100% of the outcomes). In the third step the selected pw element outcomes get combined straight following the selected condensed pw construction patterns (called 'cpattprod' inside the program). The final size of the generated <outpwfile> is steered by the specified p-value. A high p-value creates a more large pw list, a low p-value creates a (more) small one (pw list compression functionality at very low p-values). The more close the p-value gets to 1.0, the more >> over-linear will the generated <outpwfile> increase.<< As far the <inpwfile> is not trivial simple structured, high p-vales at some point will undenyably result in a 'never' ending pw generation job. Thus, it is highly recommended to start with moderate p-values first (e.g. 0.5) and to switch to interactive mode for higher p-values, in order to find individual p-values by category, that match the trade-off between the covered pw element outcome proportion versus the maximum acceptable generation job time/ maximum generated file size at best. Therefor view the printouts of the relative cumulated frequencies of every pw element category. Finally, sec officers and sysadmins can use the generated <outpwfile> to perform a simulated pw try-out on another unknown pw list, in order to derive a rough estimate of the own pws at risk proportion. However, the major information gain of the tool is the deep insight, how pws are structured and how relatively short pws can be assembled, that none the less are pretty safe. If there is interest in this direction, further complementary tools for  -Cleansing disturbing characters 'cleansetxtfile.py'  -Random sampling mega large (rockyou.txt) pw files 'samplefile.py'  - Simulated pw try-out 'f1prop-in-f2.py' - can be offered. Requirements: - Python3 (or higher) - numpy"         

args= parser.parse_args()

mode_interactive= bool(args.mode_interactive)
pval= args.pval
pval_cpatt= args.pval_cpatt
pval_let= args.pval_let
pval_num= args.pval_num
pval_spec= args.pval_spec

if pval_cpatt == "":
    pval_cpatt= pval 
if pval_let == "":
    pval_let= pval
if pval_num == "":
    pval_num= pval
if pval_spec == "":
    pval_spec= pval 

pval= float(pval)
pval_cpatt= float(pval_cpatt)
pval_let= float(pval_let)
pval_num= float(pval_num)
pval_spec= float(pval_spec)

workdir= args.workdir
ifl= args.inpwfile
ofl= args.outpwfile

wflpatt= workdir+ "/"+ "patt.dic"
wflcpatt= workdir +"/" +"cpatt.dic"
wfllet= workdir +"/" +"let.dic"
wflnum= workdir +"/" +"num.dic"
wflspec= workdir +"/" +"spec.dic"
wflcpattprod= workdir +"/" +"cpattprod.dic"
wflletprod= workdir +"/" +"letprod.dic"
wflnumprod= workdir +"/" +"numprod.dic"
wflspecprod= workdir +"/" +"specprod.dic"



################################################################################

### Functions


def read_inpwfile():

    """
    Reading the inpwfile separating the pwlines into the workfiles "patt.dic", "cpatt.dic", "let.dic", "num.dic" and "spec.dic"
    """   
    
    print("")
    print("")
    print("Reading the inpwfile   " +ifl +"   and pre-processing it ...")
    
    # Open for byte read
    fdrifl = open(ifl, "rb")
    # Open for utf-8 write
    fdwpatt= open(wflpatt, "w")
    fdwcpatt= open(wflcpatt, "w")
    fdwlet = open(wfllet, "w")
    fdwnum = open(wflnum, "w")
    fdwspec = open(wflspec, "w")


    byteline= b''
    lineraw= ''
    for byteline in fdrifl:

        try:
            # Drop out exotic language symbols
            lineraw= byteline.decode("ascii")
            lineraw= lineraw.strip("\n")  
        except UnicodeDecodeError:
            lineraw=""
            print("Ignored: ")
            print(byteline)
        
        
        ### Drop out special symbols first (for performance reasons): 
        
        lineletnum= lineraw
        for i in range (0, 48):
            lineletnum= lineletnum.replace(chr(i),"$")
        for i in range (58, 65):
            lineletnum= lineletnum.replace(chr(i),"$")
        for i in range (91, 97):
            lineletnum= lineletnum.replace(chr(i),"$")
        for i in range (123, 256):
            lineletnum= lineletnum.replace(chr(i),"$")
     
    
         ### Extract pw patttern (-> patt.dic):
        
        linepatt= lineletnum
        for i in range (48, 58):
            linepatt= linepatt.replace(chr(i),"1")
        for i in range (65, 91):
            linepatt= linepatt.replace(chr(i),"A")
        for i in range (97, 123):
            linepatt= linepatt.replace(chr(i),"A")
           
        if linepatt != "":
            fdwpatt.write(linepatt +"\n")
        
        
        ### Aggregate pw pattern line to condensed pw pattern line (-> cpatt.dic):
        
        linecpatt= linepatt
        
        linecpatt= linecpatt.replace("$$$$$$$$$$$$$$$$","$")
        linecpatt= linecpatt.replace("$$$$$$$$","$")
        linecpatt= linecpatt.replace("$$$$","$")
        linecpatt= linecpatt.replace("$$","$") 
        linecpatt= linecpatt.replace("$$","$")      
        
        linecpatt= linecpatt.replace("1111111111111111","1")
        linecpatt= linecpatt.replace("11111111","1")
        linecpatt= linecpatt.replace("1111","1")
        linecpatt= linecpatt.replace("11","1") 
        linecpatt= linecpatt.replace("11","1")      
        
        linecpatt= linecpatt.replace("AAAAAAAAAAAAAAAA","A")
        linecpatt= linecpatt.replace("AAAAAAAA","A")
        linecpatt= linecpatt.replace("AAAA","A")
        linecpatt= linecpatt.replace("AA","A") 
        linecpatt= linecpatt.replace("AA","A") 
        
        if linecpatt != "":
            fdwcpatt.write(linecpatt +"\n")
            
         
        ### Extract let symbols (-> let.dic):
        
        linelet= lineletnum.replace("$", chr(2) )
        
        for i in range (48, 58):
            linelet= linelet.replace(chr(i), chr(2) )

        linelet= linelet.replace(chr(2) *16, chr(2) )
        linelet= linelet.replace(chr(2) * 8, chr(2) )
        linelet= linelet.replace(chr(2) * 4, chr(2) )
        linelet= linelet.replace(chr(2) * 2, chr(2) )
        linelet= linelet.replace(chr(2) * 2, chr(2) )
        
        if linelet != "":
            lineletlist= linelet.split( chr(2) )       
            for i in range (0, len(lineletlist)):
                if lineletlist[i] != "":
                    fdwlet.write(lineletlist[i] +"\n")  
    
        
        ### Extract num symbols (-> num.dic):
        
        linenum= lineletnum.replace("$", chr(2) )
        
        for i in range (65, 91):
            linenum= linenum.replace(chr(i), chr(2) )
        for i in range (97, 123):
            linenum= linenum.replace(chr(i), chr(2) ) 
        
        linenum= linenum.replace(chr(2) *16, chr(2) )
        linenum= linenum.replace(chr(2) * 8, chr(2) )
        linenum= linenum.replace(chr(2) * 4, chr(2) )
        linenum= linenum.replace(chr(2) * 2, chr(2) )
        linenum= linenum.replace(chr(2) * 2, chr(2) )
    
        if linenum != "":
            linenumlist= linenum.split(chr(2) )
            for i in range (0, len(linenumlist)):
                if linenumlist[i] != "":
                    fdwnum.write(linenumlist[i] +"\n")

        
        ### Extract special symbols (-> spec.dic):
        
        linespec= lineraw.replace(chr(32), "")
        
        for i in range (48, 58):
            linespec= linespec.replace(chr(i), chr(2) )
        for i in range (65, 91):
            linespec= linespec.replace(chr(i), chr(2) )
        for i in range (97, 123):
            linespec= linespec.replace(chr(i), chr(2) )
    
        linespec= linespec.replace(chr(2) *16, chr(2) )
        linespec= linespec.replace(chr(2) * 8, chr(2) )
        linespec= linespec.replace(chr(2) * 4, chr(2) )
        linespec= linespec.replace(chr(2) * 2, chr(2) )
        linespec= linespec.replace(chr(2) * 2, chr(2) )
     
        if linespec != "":
            linespeclist= linespec.split(chr(2) )
            for i in range (0, len(linespeclist)):               
                if linespeclist[i] != "":
                    fdwspec.write(linespeclist[i] +"\n")
    
                  
    fdrifl.close()
    fdwlet.close()
    fdwnum.close()
    fdwspec.close()
    fdwpatt.close()
    fdwcpatt.close()
    
    print("")
    print("")
    print("Ready. Step [1] 'Reading the inpwfile and pre-processing it' has been completed.")
    
    
def select_cpatt():

    """
    Reading the workfile "cpatt.dic" and computing the relative frequencies.
    Selecting the condensed pw construction pattern cpatt up to critical p_value. 
    Writing the workfile "cpattprod.dic".
    """   
    
    global pval_cpatt
        

    if mode_interactive == True:
        
        print("")
        print("")
        print("Step [2a] Select condensed pw consruction pattern substrings by p-value")
        print("")    
        print("High p-values select many outcomes, thus create many pws in the end.")
        print("But if the p-value is too high, your machine will run into an 'endless' pw generation job.")
        print("Watch the length of the final result array, that will be printed subsequently.")    
        while True:
            print("")
            print("Set p-value input loop. Press 'Enter' to proceed without any changes.")
            print("")
            print("The current p-value for the condensed pw pattern is:   " +str(pval_cpatt))           
            inputvalue= input("New p-value (between 0.00 and 1.00)?: ")
            if inputvalue== "":
                break
            else:
                try:
                    pval_cpatt= float(inputvalue)
                except ValueError:
                    print("")
                    print("")
                    print("ERROR: wrong input.") 
        try:
            os.system("clear")
        except:
            pass            
        
    
    print("")
    print("")
    print("Processing the condensed pw construction patterns ...")


    cpatt00= np.array(np.loadtxt(wflcpatt, dtype="S32", delimiter=None ) )
    # print("")
    # print("\ncpatt00:")
    # print(cpatt00)
    # print(cpatt00.shape)


    cpatt00uniq= np.unique(cpatt00, return_counts=True)
    # print("")
    # print("\ncpatt00uniq:")
    # print(cpatt00uniq)
    # # print(cpatt00uniq.shape)


    (values, counts)= cpatt00uniq 
    cpatt01=np.array((values[0], counts[0]), dtype=[('value', 'S32'), ('count', 'int32')])
    for i in range(1, values.shape[0]):
        tupel= (values[i], counts[i])
        dtype= dtype=[('value','S32'),('count','int32')]
        cpatt01= np.append(cpatt01, np.array(tupel, dtype= dtype ) )
    # print("")
    # print("\ncpatt01:")
    # print(cpatt01)
    # print(cpatt01.shape)


    cpatt02= np.sort(cpatt01, order='count')
    cpatt02= cpatt02[range((values.shape[0] -1), -1, -1)]
    print("")
    print("\ncpatt02:")
    print(cpatt02)
    print(cpatt02.shape)


    cpatt03= np.cumsum(cpatt02['count']) /np.sum(cpatt02['count'])
    print("")
    print("\ncpatt03:")
    print(cpatt03)
    print(cpatt03.shape)

    cpatt04= cpatt03[cpatt03 < pval_cpatt] 
    # print("")
    # print("\ncpatt04:")
    # print(cpatt04)
    # print(cpatt04.shape)

    cpatt05= cpatt02[0: cpatt04.shape[0]+1]['value']
    # print("")
    # print("\ncpatt05:")
    # print(cpatt05)
    # print(cpatt05.shape)


    cpatt= np.array(cpatt05, dtype="S32")
    np.savetxt(wflcpattprod, cpatt, fmt="%s")
    print("")
    print("\ncpatt:")
    print(cpatt)
    print(cpatt.shape)
    print("")
    print("@ pval_cpatt:   " + str(pval_cpatt))
    print("")
    print("")
    print("Ready. Step [2a] 'Selecting the condensed pw construction patterns' has been completed.")    

    

def select_let():

    """
    Reading the workfile "let.dic" and computing the relative frequencies.
    Selecting the letter strings up to critical p_value. 
    Writing the workfile "letprod.dic".
    """

    global pval_let
    global count_let


    if mode_interactive == True:
        
        print("")
        print("")
        print("Step [2b] Select letter substrings by p-value")
        print("")    
        print("High p-values select many outcomes, thus create many pws in the end.")
        print("But if the p-value is too high, your machine will run into an 'endless' pw generation job.")
        print("Watch the length of the final result array, that will be printed subsequently.")    
        while True:
            print("")
            print("Set p-value input loop. Press 'Enter' to proceed without any changes.")
            print("")
            print("The current p-value for the letter substrings is:   " + str(pval_let))           
            inputvalue= input("New p-value (between 0.00 and 1.00)?: ")
            if inputvalue== "":
                break
            else:
                try:
                    pval_let= float(inputvalue)
                except ValueError:
                    print("")
                    print("")
                    print("ERROR: wrong input.") 
        try:
            os.system("clear")
        except:
            pass
                

    print("")
    print("")
    print("Processing the letter substrings ...")
    

    let00= np.array(np.loadtxt(wfllet, dtype="S32", delimiter=None ) )
    # print("")
    # print("\nlet00:")
    # print(let00)
    # print(let00.shape)
    
    
    let00uniq= np.unique(let00, return_counts=True)
    # print("")
    # print("\nlet00uniq:")
    # print(let00uniq)
    # # print(let00uniq.shape)
    
    
    (values, counts)= let00uniq 
    let01=np.array((values[0], counts[0]), dtype=[('value', 'S32'), ('count', 'int32')])
    for i in range(1, values.shape[0]):
        tupel= (values[i], counts[i])
        dtype= dtype=[('value','S32'),('count','int32')]
        let01= np.append(let01, np.array(tupel, dtype= dtype ) )
    # print("")
    # print("\nlet01:")
    # print(let01)
    # print(let01.shape)
    
    
    let02= np.sort(let01, order='count')
    let02= let02[range((values.shape[0] -1), -1, -1)]
    print("")
    print("\nlet02:")
    print(let02)
    print(let02.shape)
    
    
    let03= np.cumsum(let02['count']) /np.sum(let02['count'])
    print("")
    print("\nlet03:")
    print(let03)
    print(let03.shape)
    
    
    let04= let03[let03 < pval_let] 
    # print("")
    # print("\nlet04:")
    # print(let04)
    # print(let04.shape)
    
    
    let05= let02[0: let04.shape[0]+1]['value']
    # print("")
    # print("\nlet05:")
    # print(let05)
    # print(let05.shape)
    
    
    let= np.array(let05, dtype="S32")
    count_let= float(let.shape[0])
    np.savetxt(wflletprod, let, fmt="%s")
    print("")
    print("\nlet:")
    print(let)
    print(let.shape)
    print("")
    print("@ pval_let:   " + str(pval_let))
    print("")
    print("")
    print("Ready. Step [2b] 'Selecting the letter substrings' has been completed.")    

    
    
def select_num():

    """
    Reading the workfile "num.dic" and computing the relative frequencies.
    Selecting the num strings up to critical p_value. 
    Writing the workfile "numprod.dic".
    """

    global pval_num
    global count_num
    
    if mode_interactive == True:
        
        print("")
        print("")
        print("Step [2c] Select number substrings by p-value")
        print("")    
        print("High p-values select many outcomes, thus create many pws in the end.")
        print("But if the p-value is too high, your machine will run into an 'endless' pw generation job.")
        print("Watch the length of the final result array, that will be printed subsequently.")    
        while True:
            print("")
            print("Set p-value input loop. Press 'Enter' to proceed without any changes.")
            print("")
            print("The current p-value for the numercial substrings is:   " + str(pval_num))           
            inputvalue= input("New p-value (between 0.00 and 1.00)?: ")
            if inputvalue== "":
                break
            else:
                try:
                    pval_num= float(inputvalue)
                except ValueError:
                    print("")
                    print("")
                    print("ERROR: wrong input.") 
        try:
            os.system("clear")
        except:
            pass

    
    print("")
    print("")
    print("Processing the numerical substrings ...")
    
    
    num00= np.array(np.loadtxt(wflnum, dtype="S32", delimiter=None ) )
    # print("")
    # print("\nnum00:")
    # print(num00)
    # print(num00.shape)
    
    
    num00uniq= np.unique(num00, return_counts=True)
    # print("")
    # print("\nnum00uniq:")
    # print(num00uniq)
    # # print(num00uniq.shape)
    
    
    (values, counts)= num00uniq 
    num01=np.array((values[0], counts[0]), dtype=[('value', 'S32'), ('count', 'int32')])
    for i in range(1, values.shape[0]):
        tupel= (values[i], counts[i])
        dtype= dtype=[('value','S32'),('count','int32')]
        num01= np.append(num01, np.array(tupel, dtype= dtype ) )
    # print("")
    # print("\nnum01:")
    # print(num01)
    # print(num01.shape)
    
    
    num02= np.sort(num01, order='count')
    num02= num02[range((values.shape[0] -1), -1, -1)]
    print("")
    print("\nnum02:")
    print(num02)
    print(num02.shape)
    
    
    num03= np.cumsum(num02['count']) /np.sum(num02['count'])
    print("")
    print("\nnum03:")
    print(num03)
    print(num03.shape)
    
    
    num04= num03[num03 < pval_num] 
    # print("")
    # print("\nnum04:")
    # print(num04)
    # print(num04.shape)
    
    
    num05= num02[0: num04.shape[0]+1]['value']
    # print("")
    # print("\nnum05:")
    # print(num05)
    #print(num05.shape)
    
    
    num= np.array(num05, dtype="S32")
    count_num= float(num.shape[0])
    np.savetxt(wflnumprod, num, fmt="%s")
    print("")
    print("\nnum:")
    print(num)
    print(num.shape)
    print("")
    print("@ pval_num:   " + str(pval_num))
    print("")
    print("")
    print("Ready. Step [2c] 'Selecting number substrings' has been completed.")    
    
    

def select_spec():

    """
    Reading the workfile "spec.dic" and computing the relative frequencies.
    Selecting the spec strings up to critical p_value. 
    Writing the workfile "specprod.dic".
    """

    global pval_spec
    global count_spec

    if mode_interactive == True:
        
        print("")
        print("")
        print("Step [2d] Select special character substrings by p-value")
        print("")    
        print("High p-values select many outcomes, thus create many pws in the end.")
        print("But if the p-value is too high, your machine will run into an 'endless' pw generation job.")
        print("Watch the length of the final result array, that will be printed subsequently.")    
        while True:
            print("")
            print("Set p-value input loop. Press 'Enter' to proceed without any changes.")
            print("")
            print("The current p-value for the special character substrings is:   " + str(pval_spec))           
            inputvalue= input("New p-value (between 0.00 and 1.00)?: ")
            if inputvalue== "":
                break
            else:
                try:
                    pval_spec= float(inputvalue)
                except ValueError:
                    print("")
                    print("")
                    print("ERROR: wrong input.")
        try:
            os.system("clear")
        except:
            pass
                    

    print("")
    print("")
    print("Processing the special character substrings ...")  


    spec00= np.array(np.loadtxt(wflspec, dtype="S32", delimiter=None ) )
    # print("")
    # print("\nspec00:")
    # print(spec00)
    # print(spec00.shape)
    
    
    spec00uniq= np.unique(spec00, return_counts=True)
    # print("")
    # print("\nspec00uniq:")
    # print(spec00uniq)
    # # print(spec00uniq.shape)
    
    
    (values, counts)= spec00uniq 
    spec01=np.array((values[0], counts[0]), dtype=[('value', 'S32'), ('count', 'int32')])
    for i in range(1, values.shape[0]):
        tupel= (values[i], counts[i])
        dtype= dtype=[('value','S32'),('count','int32')]
        spec01= np.append(spec01, np.array(tupel, dtype= dtype ) )
    # print("")
    # print("\nspec01:")
    # print(spec01)
    # print(spec01.shape)
    
    
    spec02= np.sort(spec01, order='count')
    spec02= spec02[range((values.shape[0] -1), -1, -1)]
    print("")
    print("\nspec02:")
    print(spec02)
    print(spec02.shape)
    
    
    spec03= np.cumsum(spec02['count']) /np.sum(spec02['count'])
    print("")
    print("\nspec03:")
    print(spec03)
    print(spec03.shape)
    
       
    spec04= spec03[spec03 < pval_spec] 
    # print("")
    # print("\nspec04:")
    # print(spec04)
    # print(spec04.shape)
    
    
    spec05= spec02[0: spec04.shape[0]+1]['value']
    # print("")
    # print("\nspec05:")
    # print(spec05)
    # print(spec05.shape)
    
    
    spec= np.array(spec05, dtype="S32")
    count_spec= float(spec.shape[0])
    np.savetxt(wflspecprod, spec, fmt="%s")
    print("")
    print("\nspec:")
    print(spec)
    print(spec.shape)
    print("")
    print("@ pval_let:   " + str(pval_let))
    print("")
    print("")
    print("Ready. Step [2d] 'Selecting the special characters substrings' has been completed.")    


    
def gen_pws(ofl=ofl):


    """
    Initial runtime error catcher:
    Forecast of number of pws to generate.
    Wanna proceed break point to avoid endless generation job.  
    """


    global count_let
    global count_num
    global count_spec
    
        
    try:
        sum= 0.0
        fdrcpattprod= open(wflcpattprod, "r")  
        for cpattline in fdrcpattprod:
            prod= 1.0
            cpattlen= len(cpattline[2:-2])
            for i in range (0, cpattlen):
                cpattsymbol= cpattline[2:-2][i]
                # print(cpattsymbol) 
                if cpattsymbol== "A":
                    prod= prod * count_let    
                elif cpattsymbol== "1":
                    prod= prod * count_num
                elif cpattsymbol== "$":
                    fdwofl.write(specline[2:-2])
                    prod= prod * count_spec    
                sum= sum +prod
        
        sum_int= int(sum)   
        print("")
        print("")
        print("The size forecast for the number of pws to generate is:")
        print("")
        print(sum_int)
        print("")
        print("(That might be " + str(sum_int *15 /1000000) +" MB.)")
    
    except:
        print("")
        print("")
        print("The size forecast cannot be computed, because one of the previous steps has not been run fresh before.")
        print("")
     
        
    print("")
    prlinput= input("Do you want to proceed? (y/ n/ Crtl-c):")
    if prlinput != "y":
        exit()
    

    """
    Reading the workfiles "cpattprod.dic", "letprod.dic", "numprod.dic" and "specprod.dic"
    Concatenating the new pws in nested loops straight following the condensed pw construction pattern "cpattprod.dic"
    No configuration here. The selection config is performed in the previous steps. 
    Writing the final outfile "pwsgenerated.dic".
    """

    
    print("")
    print("")
    print("Generating pws ...")
  
           
    fdrchrprod= open(wflletprod, "r")
    fdrnumprod= open(wflnumprod, "r")
    fdrspecprod= open(wflspecprod, "r")
    chrlinelist= fdrchrprod.readlines()
    numlinelist= fdrnumprod.readlines()    
    speclinelist= fdrspecprod.readlines()            
    fdrchrprod.close()
    fdrnumprod.close()
    fdrspecprod.close()
    
    
    fdwofl = open(ofl, "w")
    fdrcpattprod= open(wflcpattprod, "r")         
    for cpattline in fdrcpattprod:
        # print(cpattline)
        haschr= False
        hasnum= False 
        hasspec= False
        if cpattline.find("A")!= -1:
            haschr= True
        if cpattline.find("1")!= -1:
            hasnum= True
        if cpattline.find("$")!= -1:
            hasspec= True
                
        fdrchrprod = open(wflletprod, "r")
        chrline= fdrchrprod.readline()
        while chrline in chrlinelist:
            # print(chrline)
            fdrnumprod = open(wflnumprod, "r")
            numline= fdrnumprod.readline()       
            while numline in numlinelist:
                # print(numline)
                fdrspecprod = open(wflspecprod, "r")
                specline= fdrspecprod.readline()
                while specline in speclinelist:
                    # print(specline)
                    cpattlen= len(cpattline[2:-2])
                    for i in range (0, cpattlen):
                        cpattsymbol= cpattline[2:-2][i]
                        # print(cpattsymbol) 
                        if cpattsymbol== "A":
                            fdwofl.write(chrline[2:-2])
                        elif cpattsymbol== "1":
                            fdwofl.write(numline[2:-2])
                        elif cpattsymbol== "$":
                            fdwofl.write(specline[2:-2])
                    fdwofl.write("\n")
                    
                    if hasspec== False:
                        break          
                    try:
                        specline= next(fdrspecprod)
                    except:
                        specline= chr(2)
                        fdrspecprod.close()
    
                if hasnum== False:
                    break
                try:
                    numline= next(fdrnumprod)
                except:
                    numline= chr(2)
                    fdrnumprod.close()
            
            if haschr== False:
                break
            try:
                chrline= next(fdrchrprod)
            except:
                chrline= chr(2)
                fdrchrprod.close()
    
    fdrcpattprod.close()
    fdwofl.close()           

    print("")
    print("")
    print("Ready. The final step [3] 'Generating new pws' has been completed.")
    print("")
    print("The <outpwfile> generated is   " +ofl + " .")   
    print("")

    
    
################################################################################

### Main

if mode_interactive == True:
    while True:
        print("")
        print("")
        print("")
        print("Interactive Mode")
        print("")
        print("What do you want to do?")
        print("[0]  Exit interactive mode")
        print("[1]  Reading inpwfile and pre-processing it")
        print("[2a] Selecting pw construction patterns by p-value")
        print("[2b] Selecting letter substrings by p-value")
        print("[2c] Selecting number substrings by p-value")
        print("[2d] Selecting special character substrings by p-value")
        print("[3]  Generate new pws")
        print("")
        print("Each step requires (all) the result files of the previous main number steps (in the workdir).")
        print("In step 'Generate new pws' the size forecast additionally requires the previous 4 partial steps been run fresh before.")
        print("")
        print("")
        selinput= input("Your selection: ")
        
        if selinput == "0":
            break
        elif selinput == "1":
            read_inpwfile() 
        elif selinput == "2a":
            select_cpatt()
        elif selinput == "2b":
            select_let()
        elif selinput == "2c":
            select_num()
        elif selinput == "2d":
            select_spec()
        elif selinput == "3":
            gen_pws()
        else:
            print("")
            print("")
            print("ERROR: Wrong input!")


else:
    read_inpwfile()
    select_cpatt()
    select_let()
    select_num()
    select_spec()
    gen_pws()

    
    
################################################################################

# Further development options:

# - Active support for multiple pw files



