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
        ----------
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
    fichier_out.write('NameReg\tNumReg\tPigbkReg\tPfgbkReg\tStrandReg\tLengReg\tPredReg\tSomPred\tDenPred\n')
#    fichier_out.write('PModXPS\tYPredWeak\tYPredStrong\tYPredNone\n')
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
    fichier_out.write(nom_fichier_out + "\n")
    fichier_out.write("Name_BP\tPosBp\tStBp\tAnType\t")
    fichier_out.write("gam\tsyngam\tfeatgam\tpigam\tpfgam\tstgam\tdeltaA\tENDgam\tOrigam\t")
    fichier_out.write("gav\tsyngav\tfeatgav\tpigav\tpfgav\tstgav\tdeltaB\tENDgav\tOrigav\n")
    fichier_out.write(('\n'.join('\t'.join(str(x) for x in l) for l in ldl)))
    fichier_out.write("\n")
    fichier_out.close()







"""
A winner-takes-it-all strategy was then used to assign each strand position to a single category
(i.e. the one with the highest mean score). For instance, strand positions for which YN[mean] >
YW[mean], YS[mean], or YOut[mean] were assigned to the ‘None’ category. The predicted None’,
‘Weak’, ‘Strong’, and ‘out-of-model’ regions were then defined as continuous successions of
strand positions falling into the corresponding category
@author: EE
"""

debut0 = time.time()

if len(sys.argv) != 4:
    print("Usage: python prog-python /PATH/file1 strand obs")
    print("python SeparefileForRev_X.py /Path/TabTexte F or R test1")
    quit()
else:
    print("\n-- Welcome --\n")

print('----------------input_output_file-------------------')

chemin = os.getcwd()
argument = sys.argv
infichier1 = chemin + "/" + argument[1]
brin = argument[2]
obs = argument[3]

# elimine le .txt
#nom_out = argument[2][0:-4]
nomsplit1 = deconcat(argument[1])
nomfile1 = nomsplit1[-1][0:-4]
print("deconcat = :", nomsplit1) 
print("nomfile1 = ", nomfile1)

nom_out_fichier1 = nomfile1 + "_" + obs + "_RegPred.txt"

nom_out_fichier3 = nomfile1 + "_" + obs + "_RegPred_infos.txt"
outfichier1 = chemin + "/" + nom_out_fichier1
outfichier3 = chemin + "/" + nom_out_fichier3

print("Directory =", chemin)
print ("arguments =", argument)
print ("infichier1 =", infichier1)
print("brin = ", brin)
print ("observation=", obs)
print ("nom_out_fichier1 =", nom_out_fichier1)
print ("outfichier1 =", outfichier1)
print ("nom_out_fichier3 =", nom_out_fichier3)
print ("outfichier3 =", outfichier3)

ldlinfos =[]

inf01 = [chemin, argument, infichier1, obs]
ldlinfos.append(inf01)

print("----------------load Data -------------------")

print ("_____  file1  ________")
list1 = Transforme_fichier_en_liste_de_ligne(infichier1)
ldl1 = Construit_liste_de_liste(list1)
print_tronc(ldl1, 3)
inf1 = ["ldl1 = ", len(ldl1)]
ldlinfos.append(inf1)
print ("_____f1________")
ldlPred = Cons_ldl_Data_1(ldl1, 1)
print_tronc(ldlPred, 3)
inf2 = ["fichier ldlPred = ", len(ldlPred)]
ldlinfos.append(inf2)

print("----------------Calcul Regions----------------")
    
ldlRegPred = []
numreg = 0
posRegI = 0
posRegF = 0
longReg = 0
predReg = 'ND'

prob1 = 0
prob2 = 0
comptPre = 0
comptDer = 0
comptMP = 0
comptDV = 0
nomregion = ""
PredValTot = 0
PredVal = 0

for l in ldlPred:
    posl = int(l[0])
    predl = str(l[10])
    
    VN = float(l[6])
    VW = float(l[7])
    VS = float(l[8])
    VE = float(l[9])
    
    if predl == "E":
        PredVal = VE
    elif predl == "S":
        PredVal = VS
    elif predl == "W":
        PredVal = VW
    elif predl == "N":
        PredVal = VN
    else:
        prob2 += 1
    
    
    
    
    
    print(predl, VN, VW, VS, VE, PredVal, PredValTot)
    if posl == 1:
        #initialisation variable
        numreg = 1
        posRegI = posl
        posRegF = posl
        longReg = 1
        predReg = predl
#        print("premiere reg" , numreg, posRegI, predReg)
        comptPre = comptPre + 1
        PredValTot = PredVal
        DenPred = 0
        
    elif posl == len(ldlPred):
        #fin de fichier
        nomregion = "Reg_" + str(numreg) + "_" + str(posRegI) + "_" + str(posRegF + 1) + "_" + brin
        DenPred = round((PredValTot / (longReg + 1)), 2)
        ldlRegPred.append([nomregion, numreg, posRegI, (posRegF + 1), brin, (longReg + 1), predReg, round(PredValTot, 2), DenPred])
#        print("derniere reg", numreg, posRegI, (posRegF + 1), (longReg + 1), predReg)
        comptDer = comptDer + 1
    
    elif predl == predReg:
        # meme variable pred
        longReg = longReg + 1
        posRegI = posRegI
        posRegF = posl
        comptMP = comptMP + 1
        predReg = predReg
        PredValTot = PredValTot + PredVal
        
    elif predl != predReg:
        # variable Pred differente
        # Calcul densité        
        DenPred = round((PredValTot / longReg), 2)
        nomregion = "Reg_" + str(numreg) + "_" + str(posRegI) + "_" + str(posRegF) + "_" + brin
        ldlRegPred.append([nomregion , numreg, posRegI, posRegF, brin, longReg, predReg, round(PredValTot, 2), DenPred])
#        print(numreg, posRegI, posRegF, longReg, predReg)
        comptDV = comptDV + 1
        numreg = numreg + 1
        posRegI = posl
        posRegF = posl
        longReg = 1
        predReg = predl
        PredValTot = PredVal
        
    else:
        prob1 += 1
        
    
#print(prob1, comptMP, comptDV)
    
print_tronc(ldlRegPred, 5)

        
inf3 = ["Probleme boucle else, nombre de ligne:", prob1]
ldlinfos.append(inf3)

inf4 =["passages differentes boucles:", "1 ligne=", comptPre, "meme pred=", comptMP, "diff pred=", comptDV, "der ligne=", comptDer ]
ldlinfos.append(inf4)

inf5 = ["longueur ldlRegPred = ", len(ldlRegPred)]
ldlinfos.append(inf5)

inf6 = ["Probleme boucle else, Cherche PredVal:", prob2]
ldlinfos.append(inf6)
#print_tronc(Fldl, 3)
print("-----------------------------------")
#print_tronc(Rldl, 3)


print("---------------ldl infos-------------------")

print(ldlinfos)

print("----------------Ecriture Fichier-------------------")

ecrit_ldl_2(ldlRegPred, outfichier1)
#ecrit_ldl_2(Rldl, outfichier2)
#ecrit_ldl_3(ldlinfos, outfichier3)


fin0 = time.time()

print("Total time: %.4f sec:" % (fin0 - debut0))



