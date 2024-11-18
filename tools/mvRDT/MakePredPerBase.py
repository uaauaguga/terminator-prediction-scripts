# -*- coding: utf-8 -*-
from __future__ import print_function
import os
import sys
import time

def deconcat(ligne, sep = '/'):
    '''
    Programmer help
         Parameters
         ----------
        a :(string)
        Returns
        -------
        output screen : string
    '''
    elems = ligne.split(sep)
    return elems
    
def Transforme_fichier_en_liste_de_ligne(mon_fichier):
    '''
    Transforms file on line list
        Parameters
        ----------
        a :root/name_file
        Returns
        -------
        a list of line read
        --------
    '''
    fichier = open(mon_fichier, "r")
    ma_liste = fichier.readlines()
    fichier.close()
    return ma_liste


def Construit_liste_de_liste(ma_liste):
    '''
    Transforms line list in list of list
        Parameters
        ----------
        a :list
        Returns
        -------
        a list of list (ldl)
        Example
        -------
        in = ['AAAA\t7175\n', 'TTTT\t2609\n']
        out = [['AAAA','7175'],['TTTT','2609']]
    '''
    ma_liste_out = []
    for line in ma_liste:
        line2 = line.rstrip()
        liste2 = line2.split()
        ma_liste_out.append(liste2)
    return ma_liste_out


def print_tronc(ldl1, int1):
    '''
    screen displays the first and last line of ldl and number of ldl
    ----------       
       Parameters
        ----------
        a :list of list (ldl)
        b : int
        Returns
        -------
        res screen:print infos
        Example
        -------
        >>> print_tronc(l1,3)
        [1, 12, 41, 'TCAACAAACGTATTACCGCCAGGTAAAGAA', 30]
        [2, 27, 56, 'CCGCCAGGTAAAGAACCCGAATCCGGTGTT', 30]
        [3, 30, 57, 'CCAGGTAAAGAACCCGAATCCGGTGTTT', 28]
        ...length ldl= 2609
        [2607, 63313, 63339, 'TCAAGACTTCTTTCTGTGCTCACTCCT', 27]
        [2608, 63321, 63350, 'TCTTTCTGTGCTCACTCCTTCTGCGCATTG', 30]
        [2609, 63325, 63354, 'TCTGTGCTCACTCCTTCTGCGCATTGTAAG', 30]
        ...
        Rq: if length ldl < int1 ERROR....Attention
    '''
    i = 0
    while i < int1:
        print(ldl1[i])
        i += 1
    print ("...length ldl=", len(ldl1))
    j = -int1
    while j < 0:
        print(ldl1[j])
        j += 1

def print_tronc2(l, int1, int2):
    '''
    screen displays the first and last line of ldl and number of ldl
    ----------       
       Parameters
        ----------
        a :list of list (ldl)
        b : int
        Returns
        -------
        res screen:print infos
        Example
        -------
        >>> print_tronc(l1,3)
        [1, 12, 41, 'TCAACAAACGTATTACCGCCAGGTAAAGAA', 30]
        [2, 27, 56, 'CCGCCAGGTAAAGAACCCGAATCCGGTGTT', 30]
        [3, 30, 57, 'CCAGGTAAAGAACCCGAATCCGGTGTTT', 28]
        ...length ldl= 2609
        [2607, 63313, 63339, 'TCAAGACTTCTTTCTGTGCTCACTCCT', 27]
        [2608, 63321, 63350, 'TCTTTCTGTGCTCACTCCTTCTGCGCATTG', 30]
        [2609, 63325, 63354, 'TCTGTGCTCACTCCTTCTGCGCATTGTAAG', 30]
        ...
        Rq: if length ldl < int1 ERROR....Attention
    '''
    i = int(int1)
    while i < int(int2):
        print(l[i])
        i += 1
    print ("...len=", len(l))
    

def print_1(l):
    
    print(l[0])
    print(l[50])
    print(l[51])
    print(l[99])
    print(l[101])
    print(l[300])
    print(l[700])
    print(l[1000])



def Cons_ldl_Data(ldl, nbr1, nbr2, nbr3, cutoff, nbr4):
    '''
    Parameters
        ----------
        1 : list of list
        2 : nbr1 = int (index column matrice positioni Pi ADN)
        3 : nbr2 = int (index column matrice positionf PF ADN)
        4 : nbr3 = int (ligne a pas prendre en compte dans la ldl)
        5 : cutoff1 = pour selectionner list ok
        6 : nbr4 = int (index column indiquand le brin 1 ou -1)
        
    Returns
        -------
        return = ldlDATAOUT
    Example
        -------

    '''
   
    ldlout =[]   
    i = nbr3
    while i < len(ldl):
            posini = int(ldl[i][nbr1])
            posfin = int(ldl[i][nbr2])
            longRSB = (posfin - posini)
            # print(longBulle)
            
            
            if longRSB >= int(cutoff):
                ldlout.append(ldl[i])
            i += 1
    return ldlout


def Cons_ldl_Data_1(ldl, nbr3):
    '''
    enleve ligne d'entete
    Parameters
        ----------
        1 : list of list
        2 : nbr3 = int (ligne a pas prendre en compte dans la ldl)
        
    Returns
        -------
        return = ldlDATAOUT
    Example
        -------

    '''
   
    ldlout =[]   
    i = nbr3
    while i < len(ldl):
        ldlout.append(ldl[i])        
        i += 1
    return ldlout


def ecrit_ldl_2(ldl, nom_fichier_out):
    '''
    writes a column LDL in a separate file with tab
        Parameters
        ----------Separation file
        a : list of list (ldl)
        b : string /absolute path/
        Returns
        -------
        writes files separate whit /tab in nom_fichier_out
       PrimaryID	IDACC	num-files	IDSeq	Genom	PiGbk	PfGbk	StrandL	StrandN	PModXPS2	tPS1	tPS2	DModXPS2NWR	CItPS1	YPredWeak	YPredStrong	YPredPS2Nothing
comptMatch [ListBST]
    '''
    print(" files write in : ", nom_fichier_out)
    fichier_out = open(nom_fichier_out, "w")
#    fichier_out.write(nom_fichier_out + "\n")
    fichier_out.write('name_Seq\tGenom\tPigbk\tPfgbk\tStrand\tIDNumSeq\tSeqNum\t')
    fichier_out.write('PModXPS\tYPredWeak\tYPredStrong\tYPredNone\tPredExclud\n')
    fichier_out.write(('\n'.join('\t'.join(str(x) for x in l) for l in ldl)))
    fichier_out.write("\n")
    fichier_out.close()


def ecrit_ldl_3(ldl, nom_fichier_out):
    '''
    writes a column LDL in a separate file with tab
        Parameters
        ----------
        a : list of list (ldl)
        b : string /absolute path/
        Returns
        -------
        writes files separate whit /tab in nom_fichier_out
    '''
    print(" files write in : ", nom_fichier_out)
    fichier_out = open(nom_fichier_out, "w")
    fichier_out.write(nom_fichier_out + "\n")
    fichier_out.write(('\n'.join('\t'.join(str(x) for x in l) for l in ldl)))
    fichier_out.write("\n")
    fichier_out.close()


def ecrit_ldl_4(ldl, nom_fichier_out):
    '''
    writes a column LDL in a separate file with tab
        Parameters
        ----------
        a : list of list (ldl)
        b : string /absolute path/
        Returns
        -------
        writes files separate whit /tab in nom_fichier_out
    '''
    print(" files write in : ", nom_fichier_out)
    fichier_out = open(nom_fichier_out, "w")
#    fichier_out.write(nom_fichier_out + "\n")
    fichier_out.write("PosGenom\tValNone\tValWeak\tValStrong\tValExcluded\tNbrSeq\t")
    fichier_out.write("VNoneNor\tVWeakNor\tVStrongNor\tVExcluNor\tValPred\n")
#    fichier_out.write("gav\tsyngav\tfeatgav\tpigav\tpfgav\tstgav\tdeltaB\tENDgav\tOrigav\n")
    fichier_out.write(('\n'.join('\t'.join(str(x) for x in l) for l in ldl)))
    fichier_out.write("\n")
    fichier_out.close()


def Zerolistmaker(n):
    listofzeros = [0] * n
    return listofzeros



"""
Created on Tue Jan 31 06:48:38 2017
Associates a prediction score at each base of a given genome
Sequences with PmodXPS+ < 0.05 (model outliers at 95%
confidence level), the YOut score was set to 1 while the YN, YW, and YS scores were
conservatively changed to zero. For other sequences (PmodXPS+ ≥ 0.05), the YOut score was set
to zero while the YN, YW, and YS scores were kept at their original YPredPS values. Then, for
each genomic strand position, mean scores YN[mean], YW[mean], YS[mean], and YOut[mean] were
calculated, respectively, as the averages of the YN, YW, YS and YOut scores obtained for all of
the 500 nt-long, same-strand sequences overlapping this position:
@author: EE

"""

debut0 = time.time()

# checks if only 2 argument was passed to call the script
if len(sys.argv) != 3:
    print("Usage: python prog-python /PATH/file1 obs")
    print("python TabPredPlusColExclud_1.py /Path/TabTexte test1")
    quit()
else:
    print("\n-- Welcome --\n")

print('----------------input_output_file-------------------')

chemin = os.getcwd()
argument = sys.argv
infichier1 = chemin + "/" + argument[1]
obs = argument[2]

# elimine le .txt
#nom_out = argument[2][0:-4]
nomsplit1 = deconcat(argument[1])
nomfile1 = nomsplit1[-1][0:-4]
print("deconcat = :", nomsplit1) 
print("nomfile1 = ", nomfile1)

nom_out_fichier1 = nomfile1 + "_" + obs + "_ColExclude.txt"
nom_out_fichier2 = nomfile1 + "_" + obs + "_GenomValTot.txt"
nom_out_fichier3 = nomfile1 + "_" + obs + "_Genom_infos.txt"
outfichier1 = chemin + "/" + nom_out_fichier1
outfichier2 = chemin + "/" + nom_out_fichier2
outfichier3 = chemin + "/" + nom_out_fichier3

print("Directory =", chemin)
print ("arguments =", argument)
print ("infichier1 =", infichier1)
print ("observation=", obs)
print ("nom_out_fichier1 =", nom_out_fichier1)
print ("outfichier1 =", outfichier1)
#print ("nom_out_fichier2 =", nom_out_fichier2)
#print ("outfichier2 =", outfichier2)
print ("nom_out_fichier3 =", nom_out_fichier3)
print ("outfichier3 =", outfichier3)

LGenom = 10000
#LGenome Ecoli , Change Lengh Genome for LT2 


ldlinfos =[]

inf01 = ["la longueur du genom = ", LGenom]
ldlinfos.append(inf01)

print("----------------load Data -------------------")

print ("_____  file1  ________")
list1 = Transforme_fichier_en_liste_de_ligne(infichier1)
# a list of lines
ldl1 = Construit_liste_de_liste(list1)
# a list of entries
#print_tronc(ldl1, 3)
inf1 = ["ldl1 = ", len(ldl1)]
ldlinfos.append(inf1)
print ("_____f1________")
ldlPred = Cons_ldl_Data_1(ldl1, 1)
# skip  first line
#print_tronc(ldlPred, 3)
inf2 = ["fichier ldlPred = ", len(ldlPred)]
ldlinfos.append(inf2)


print("------------Rajoute Col Exclud--------------")    



probl = 0
nbrOK = 0
nbrEx = 0
ldlColE =[]
ltmp = []

for l in ldlPred:
    ltmp = []

    if float(l[7]) >= 0.05:
        ltmp = [l[0], l[1], l[2], l[3], l[4], l[5], l[6], l[7], l[8], l[9], l[10], 0 ]
        ldlColE.append(ltmp)
        nbrOK += 1
    elif float(l[7]) < 0.05:
        ltmp = [l[0], l[1], l[2], l[3], l[4], l[5], l[6], l[7], 0 , 0 , 0 , 1 ]
        ldlColE.append(ltmp)
        nbrEx += 1
    else :
        probl +=1
        
inf3 = ["Probleme boucle else , nombre de ligne:", probl]
ldlinfos.append(inf3)

inf4 =["nbr>= ", nbrOK, "nbr < 0.05", nbrEx]

print_tronc(ldlColE, 10)

ldlinfos.append(inf4)
print("------------Construit Files Genoms--------------")  
inf5 = ["fichier ldlColE = ", len(ldlColE)]
ldlinfos.append(inf5)

print("------------Construit Arrays Genoms--------------")    

print(LGenom)

lPosGenom = Zerolistmaker(LGenom)
lNone = Zerolistmaker(LGenom)
lWeak = Zerolistmaker(LGenom)
lStrong = Zerolistmaker(LGenom)
lExcluded = Zerolistmaker(LGenom)
lNbrSeq = Zerolistmaker(LGenom)

print("remplie lPosGenom")

k = 0
while k < LGenom :
    lPosGenom[k] = k + 1
    k += 1

inf6 = ["list lPosGenom = ", len(lPosGenom)]
ldlinfos.append(inf6)



print("remplie lists W S N E")


print("$$$$$$$$$$$",len(ldlColE))
i = 0
while i < len(ldlColE):
#    print(i)
#    print(ldlColE[i])
    posIniSeq = int(ldlColE[i][2])
    posFinSeq = int(ldlColE[i][3])
    valW = float(ldlColE[i][8])
    valS = float(ldlColE[i][9])
    valN = float(ldlColE[i][10])
    valE = float(ldlColE[i][11])
    print(posIniSeq,posFinSeq,valW,valS,valN,valE)
    
# remplies les listes   
    for j in range((posIniSeq-1), posFinSeq):
#        print(j)
        lNone[j] = lNone[j] + valN
        lWeak[j] = lWeak[j] + valW
        lStrong[j] = lStrong[j] + valS
        lExcluded[j] = lExcluded[j] + valE
        lNbrSeq[j] = lNbrSeq[j] + 1
            
# fin boucle    
    i += 1

#print("#"*100)
#print(lNone[:10])
#print(lStrong[:10])

inf7 = ["list lNone = ", len(lNone)]
ldlinfos.append(inf7)

inf8 = ["list lWeak = ", len(lWeak)]
ldlinfos.append(inf8)

inf9 = ["list lStrong = ", len(lStrong)]
ldlinfos.append(inf9)

inf10 = ["list lExcluded = ", len(lExcluded)]
ldlinfos.append(inf10)

inf11 = ["list lNbrSeq = ", len(lNbrSeq)]
ldlinfos.append(inf11)

print("------------Construit Files Genoms--------------")  

ldlGenom = []
ldltmp = []

l = 0
while l < LGenom:
    ldltmp = [lPosGenom[l], lNone[l], lWeak[l], lStrong[l], lExcluded[l], lNbrSeq[l]]
    ldlGenom.append(ldltmp)
    ldltmp =[]
    
    
    l += 1

inf12 = ["ldlGenom = ", len(ldlGenom)]
ldlinfos.append(inf12)

print("------------Normalise_file_Genom--------------")

ldlGenomNorm = []
Prob1 = 0
comptE = 0
comptN = 0
comptW = 0
comptS = 0

for lf in ldlGenom:
    try:
        posGen = lf[0]
        VNoneNor = round(float(lf[1])/int(lf[5]), 4)
        VWeakNor = round(float(lf[2])/int(lf[5]), 4)
        VStrongNor = round(float(lf[3])/int(lf[5]), 4)
        VExcluNor = round(float(lf[4])/int(lf[5]), 4)
    except:
        pass
    
    # calcul de la valeur predictive c'est la Valeur la plus forte entre VN, VS, VW, VE.
    # ordre choisi Exclud >= None >= Weak >= Strong
    ValPred = "ND"
    
    if (VExcluNor >= VNoneNor) and (VExcluNor >= VWeakNor) and (VExcluNor>= VStrongNor) : 
        ValPred = "E"
        comptE += 1
    elif (VNoneNor > VExcluNor) and (VNoneNor >= VWeakNor) and (VNoneNor >= VStrongNor) :
        ValPred = "N"
        comptN += 1
    elif (VWeakNor > VExcluNor) and (VWeakNor > VNoneNor) and (VWeakNor >= VStrongNor):
        ValPred = "W"
        comptW += 1
    elif (VStrongNor > VExcluNor) and (VStrongNor > VNoneNor) and (VStrongNor > VWeakNor):
        ValPred = "S"
        comptS += 1
    else:
        Prob1 += 1    
    
    ltmp2 = [posGen, round(lf[1], 4), round(lf[2], 4), round(lf[3], 4), round(lf[4], 4), int(lf[5]), \
                VNoneNor, VWeakNor, VStrongNor, VExcluNor, ValPred]
    ldlGenomNorm.append(ltmp2)



print("boucle if prediction probleme:", Prob1)
print(comptE, comptN, comptW, comptS, (comptE + comptN + comptW + comptS))

infprob1 = ["InfoProb1 Normalise file genome , boucle if nombre else =", Prob1]
ldlinfos.append(infprob1)

inf13 = ["ldlGenoNorm = ", len(ldlGenomNorm)]
ldlinfos.append(inf13)
    
inf14 = ["comptE=", comptE, "comptN=", comptN, "comptW=", comptW, "comptS=", comptS, "Tot=", (comptE + comptN + comptW + comptS)]
ldlinfos.append(inf14)



print("------------ldlColeE--------------")

#for i in range(0,12):
#    print(ldlColE[i])

print("------------lPosGenom--------------")
#print_1(lPosGenom)
print("------------lNone-------------")    
#print_1(lNone)
print("------------lWeak-------------")    
#print_1(lWeak)
print("------------lStrong-------------")    
#print_1(lStrong)
print("------------lExcluded-------------")    
#print_1(lExcluded)
print("------------lNbrSeq-------------")    
#print_1(lNbrSeq)
print("------------ldl Genom-------------")    
print_tronc(ldlGenom, 10)
print("------------ldl Genom Norm-------------")    
print_tronc(ldlGenomNorm, 10)



print("---------------ldl infos-------------------")

print(ldlinfos)

print("----------------Ecriture Fichier-------------------")

#ecrit_ldl_2(ldlColE, outfichier1)
ecrit_ldl_4(ldlGenomNorm, outfichier2)
#ecrit_ldl_3(ldlinfos, outfichier3)


fin0 = time.time()

print("Total time: %.4f sec:" % (fin0 - debut0))



