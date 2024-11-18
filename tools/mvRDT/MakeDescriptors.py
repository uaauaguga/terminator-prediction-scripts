#-*- coding: utf-8 -*-
from __future__ import print_function
import os
import sys
import time
from Bio import SeqIO
# from Bio.Seq import Seq
from re import compile
import numpy as np


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


def Comptes_Bases_Seq(my_seq, compt1, names):
    
    maldl =[]    
    #print("my_seq =", type(my_seq))    
    seqstr = str(my_seq)    
    #print("seqstr =", type(seqstr))
    longseq = len(my_seq)
    
    nbrA = my_seq.count("A")
    nbrT = my_seq.count("T")
    nbrG = my_seq.count("G")
    nbrC = my_seq.count("C")

    PcA =round(((100 * float(nbrA))/len(my_seq)), 2)
    PcT =round(((100 * float(nbrT))/len(my_seq)), 2)
    PcG =round(((100 * float(nbrG))/len(my_seq)), 2)
    PcC =round(((100 * float(nbrC))/len(my_seq)), 2)
    
    PcApT =round((100 * float(my_seq.count("A") + my_seq.count("T")) / len(my_seq)), 2) 
    PcApG =round((100 * float(my_seq.count("A") + my_seq.count("G")) / len(my_seq)), 2)
    PcApC =round((100 * float(my_seq.count("A") + my_seq.count("C")) / len(my_seq)), 2)
    PcTpG =round((100 * float(my_seq.count("T") + my_seq.count("G")) / len(my_seq)), 2)
    PcTpC =round((100 * float(my_seq.count("T") + my_seq.count("C")) / len(my_seq)), 2)
    PcGpC =round((100 * float(my_seq.count("G") + my_seq.count("C")) / len(my_seq)), 2)
    
    nbrAA = my_seq.count("AA")
    nbrAT = my_seq.count("AT")
    nbrAG = my_seq.count("AG")
    nbrAC = my_seq.count("AC")
    
    nbrTA = my_seq.count("TA")
    nbrTT = my_seq.count("TT")
    nbrTG = my_seq.count("TG")
    nbrTC = my_seq.count("TC")
    
    nbrGA = my_seq.count("GA")
    nbrGT = my_seq.count("GT")
    nbrGG = my_seq.count("GG")
    nbrGC = my_seq.count("GC")    
    
    nbrCA = my_seq.count("CA")
    nbrCT = my_seq.count("CT")
    nbrCG = my_seq.count("CG")
    nbrCC = my_seq.count("CC")
    
#    nbrCTCT = my_seq.count("CTCT")    
#    nbrTCTT = my_seq.count("TCTT")
#    nbrCTTC = my_seq.count("CTTC")
#    nbrCTTCC = my_seq.count("CTTCC")
#    nbrCTTTC = my_seq.count("CTTTC")
    
    PcAA = round(100 * float(nbrAA)/(len(my_seq) -1), 2)
    PcAT = round(100 * float(nbrAT)/(len(my_seq) -1), 2)
    PcAG = round(100 * float(nbrAG)/(len(my_seq) -1), 2)
    PcAC = round(100 * float(nbrAC)/(len(my_seq) -1), 2)

    PcTA = round(100 * float(nbrTA)/(len(my_seq) -1), 2)
    PcTT = round(100 * float(nbrTT)/(len(my_seq) -1), 2)
    PcTG = round(100 * float(nbrTG)/(len(my_seq) -1), 2)
    PcTC = round(100 * float(nbrTC)/(len(my_seq) -1), 2)

    PcGA = round(100 * float(nbrGA)/(len(my_seq) -1), 2)
    PcGT = round(100 * float(nbrGT)/(len(my_seq) -1), 2)
    PcGG = round(100 * float(nbrGG)/(len(my_seq) -1), 2)
    PcGC = round(100 * float(nbrGC)/(len(my_seq) -1), 2)

    PcCA = round(100 * float(nbrCA)/(len(my_seq) -1), 2)
    PcCT = round(100 * float(nbrCT)/(len(my_seq) -1), 2)
    PcCG = round(100 * float(nbrCG)/(len(my_seq) -1), 2)
    PcCC = round(100 * float(nbrCC)/(len(my_seq) -1), 2)



    
    maldl.append([compt1, names, longseq, \
        nbrA, nbrT, nbrG, nbrC, \
        PcA, PcT, PcG, PcC, \
        PcApT, PcApG, PcApC, PcTpG, PcTpC, PcGpC, \
        nbrAA, nbrAT, nbrAG, nbrAC, \
        nbrTA, nbrTT, nbrTG, nbrTC, \
        nbrGA, nbrGT, nbrGG, nbrGC, \
        nbrCA, nbrCT, nbrCG, nbrCC, \
        PcAA, PcAT, PcAG, PcAC, \
        PcTA, PcTT, PcTG, PcTC, \
        PcGA, PcGT, PcGG, PcGC, \
        PcCA, PcCT, PcCG, PcCC, \
#        nbrCTCT, nbrTCTT, nbrCTTC, nbrCTTCC, nbrCTTTC, \
        seqstr])
    
    return maldl

def Recherche_Multi_RegExp(compt0, name, seqstring):
    
    TabRegExp =["(T|C)C", \
    "(T|C)C[ATGC]{9,13}(T|C)C", \
    "(T|C)C[ATGC]{9,13}(T|C)C[ATGC]{9,13}(T|C)C", \
    "(T|C)C[ATGC]{9,13}(T|C)C[ATGC]{9,13}(T|C)C[ATGC]{9,13}(T|C)C", \
    "(T|C)C[ATGC]{9,13}(T|C)C[ATGC]{9,13}(T|C)C[ATGC]{9,13}(T|C)C[ATGC]{9,13}(T|C)C", \
    "(T|C)C[ATGC]{9,13}(T|C)C[ATGC]{9,13}(T|C)C[ATGC]{9,13}(T|C)C[ATGC]{9,13}(T|C)C[ATGC]{9,13}(T|C)C"]
    
    longSeq =len(seqstring)
    comptre = 0
    myldl = [compt0, name, seqstring]
    for re in TabRegExp:
        myldlseq = []
        myseqcount= []
        mydensi = 0
        comptre += 1 
        myldlseq, myseqcount = Recherche_RegExp(re, seqstring)
        mydensi = (float(myseqcount)/float(longSeq)) * 1000
        Rmydensi = round(mydensi, 2)
        #print(compt0, name, comptre, re, myseqcount, Rmydensi)
        myldl.append([comptre, re, myseqcount, Rmydensi])
    
    return myldl

def Recherche_RegExp(pattern, sequence):
    '''
    looking for a pattern of a regular expression in a sequence;
    research with overlap.
        Parameters
        ----------
        a :Regular expression "String"
        b : String (sequence)
        Returns
        -------
        returns a: list of list :[[num1,pos_ini1,position_fin1,pattern,],[...]]
                b: number int
        Example
        -------
        >>>recherche_RegExp(((T|C)C[ATGC]{10,13}){2}, sequence)
        ldl = [[1, 12, 41, 'TCAACAAACGTATTACCGCCAGGTAAAGAA', 30],
        [2, 27, 56, 'CCGCCAGGTAAAGAACCCGAATCCGGTGTT', 30],
        [3, 30, 57, 'CCAGGTAAAGAACCCGAATCCGGTGTTT', 28]]
        count = 3
    '''
    motcompil = compile(pattern)
    count = 0
    pos = 0
    liste_out = []
    while True:
        result = motcompil.search(sequence, pos)
        if result is None:
            break
        posend = result.end()
# print posend
        pos = result.start()+1
# print pos
        count += 1
        LongPattern = len(sequence[pos-1: posend])
        liste2 = [count, pos, posend, sequence[pos-1: posend], LongPattern]
# print liste2
        liste_out.append(liste2)
    return liste_out, count

    


def Re_ecrit_ldl_1(ldl):
    '''    
    writes a LDL simple utilisable for founc ecrit_ldl()
        Parameters
        ----------
        a : list of list (ldl complexe)
        
        Returns
        -------
        a : list of list (ldl simple)
    '''
    
    outldl = []
    
    for ll in ldl:
        ldlfille = [ll[0], ll[1], \
            ll[3][2], ll[3][3], ll[4][2], ll[4][3], \
            ll[5][2], ll[5][3], ll[6][2], ll[6][3], \
            ll[7][2], ll[7][3], ll[8][2], ll[8][3], \
            ll[2]]
            
        outldl.append(ldlfille)
    
    return outldl

        


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



def Recherche_Multi_Bulle (mseq, fene1):
    '''
     "number bubble =", ctCGP1
     "number bubble == -1 =", ctCGP2
     "number bubble != -1 =", ctCGP3

        
    '''

    
    atgc_count_per_window1 = nombre_atgc_par_fenetre(mseq, fene1)
    #print("atgc_count_per_window1 " , compt0 , name)
    #print(atgc_count_per_window1[0],len(atgc_count_per_window1), atgc_count_per_window1[-1])
    atgc_pourcent_per_window1 = calcul_pourcentage(atgc_count_per_window1, fene1)
    #print(atgc_pourcent_per_window1[0],len(atgc_pourcent_per_window1), atgc_pourcent_per_window1[-1])
    CmoinsG = XmoinsY(atgc_pourcent_per_window1, 3, 2)    
    #print(CmoinsG[0], len(CmoinsG), CmoinsG[-1])
    CmoinsGTri1 = np.where(CmoinsG >= 0, CmoinsG, -1)
    #print(CmoinsGTri1[0], len(CmoinsGTri1), CmoinsGTri1[-1])
    TabCGTri1 = np.column_stack((atgc_count_per_window1, atgc_pourcent_per_window1, CmoinsGTri1))    
    #print(TabTri1[0], len(TabTri1), TabTri1[-1])
    
    
    atgc_compt_pourcent_window1 = np.column_stack((atgc_count_per_window1, atgc_pourcent_per_window1))
    #print(atgc_compt_pourcent_window1[0],len(atgc_compt_pourcent_window1), atgc_compt_pourcent_window1[-1])
    base_atgc_window1 = Concat_Base_ATGC(atgc_compt_pourcent_window1, mseq)
    #print(base_atgc_window1)

    bulleCGTri1, ctCGP1, ctCGP2, ctCGP3 = calculBulle3(TabCGTri1, fene1)   
    
    return bulleCGTri1, base_atgc_window1, ctCGP1, ctCGP2, ctCGP3



def nombre_atgc_par_fenetre(seq1, window):
    '''
    Parameters
        ----------
        a : sequence (format biopython)
        b : int, base number that you want to read in windows
    Returns
        -------
        res screen:print infos
        return array with number of A, T, G, or C in the window.
        moves is 1pb/cycle
    Example
        -------
        sreen :longueur de ma sequence 63465
        le type de sequence1 = <type 'numpy.ndarray'>
        ['G' 'T' 'G' ..., 'A' 'A' 'C']
        la taille de la fenetre= 100
        le type de atgc_count_per_window1 est: <type 'numpy.ndarray'>

        return : array (ldl)
        [31 20 23 26]
        [31 21 22 26]
        [31 20 23 26]
        [32 20 22 26]
        ...length ldl= 63366
        [29 21 25 25]
        [30 21 24 25]
        [31 21 23 25]
        [31 20 23 26]
    '''
    sequence1 = np.array(list(seq1))
    gc_count = np.zeros((len(sequence1)+1, 4), np.int)

    gc_count[1:, 0] = np.add.accumulate(sequence1 == 'A')
    gc_count[1:, 1] = np.add.accumulate(sequence1 == 'T')
    gc_count[1:, 2] = np.add.accumulate(sequence1 == 'G')
    gc_count[1:, 3] = np.add.accumulate(sequence1 == 'C')

    gc_count_per_window = gc_count[window:]-gc_count[:-window]

    return gc_count_per_window


def calcul_pourcentage(ldl, nbr):
    '''
    Parameters
        ----------
        1 : list of list <type 'numpy.ndarray'>
        2 : nbr = int, in this script
        base number that you want to read in windows (78)
    Returns
        -------
        res screen:print infos
        return b = array with % ofnumber of A, T, G, or C in the window.
        <type 'numpy.ndarray'>
    Example
        -------
            a = count
            [23 16 18 21]
            [24 16 17 21]
            [24 15 17 22]
            [24 15 17 22]
            b = pourcent
            [ 29.48717949  20.51282051  23.07692308  26.92307692]
            [ 30.76923077  20.51282051  21.79487179  26.92307692]
            [ 30.76923077  19.23076923  21.79487179  28.20512821]
            [ 30.76923077  19.23076923  21.79487179  28.20512821]
    '''
    a = ldl
    num1 = float(nbr)
    b = (a * 100) / num1
    return b


def XmoinsY(ldl, nbr1, nbr2):
    '''
    Parameters
        ----------
        1 : list of list <type 'numpy.ndarray'>
        2 : nbr1 = int (index column matrice)
        3 : nbr2 = int (index column matrice)
    Returns
        -------
        return = array <type 'numpy.ndarray'>
    Example
        -------

    '''
    ar1 = ldl
    arnbr1 = ar1[:, nbr1]
    arnbr2 = ar1[:, nbr2]
    arnbr1moinsarnbr2 = arnbr1 - arnbr2
    return arnbr1moinsarnbr2


def Concat_Base_ATGC(ldl, seq):
    '''
        
        Parameters
        ----------
        1 : list of list <type 'numpy.ndarray'>
        [[ 23.          16.          18.         ...,  20.51282051  23.07692308
        26.92307692]
        [ 24.          16.          17.         ...,  20.51282051  21.79487179
        26.92307692]
        [ 24.          15.          17.         ...,  19.23076923  21.79487179
        28.20512821]
        ..., 
        [ 22.          17.          18.         ...,  21.79487179  23.07692308
        26.92307692]
        [ 22.          17.          18.         ...,  21.79487179  23.07692308
        26.92307692]
        [ 22.          17.          18.         ...,  21.79487179  23.07692308
   26.92307692]]
        2 : seq =  ATGC "ATCCTGCA"
        Type my_seq=  <class 'Bio.Seq.Seq'>
        Type ar=  <type 'numpy.ndarray'>
        length my seq=  63442

        Returns
        -------
        return = ldl list of list
        [[1, 'G', 23.0, 16.0, 18.0, 21.0, 29.487179487179485,..., 26.923076923076923]
        [2, 'T', 24.0, 16.0, 17.0, 21.0, 30.76923076923077, ..., 26.923076923076923]
        [3, 'G', 24.0, 15.0, 17.0, 22.0, 30.76923076923077, ..., 28.205128205128204]]

        Example
        -------
         base_atgc_window1 = Concat_Base_ATGC(atgc_compt_pourcent_window1, brin_plus_seq_cir_mod)   
            '''


    ar = ldl
    my_seq = seq
    outldl = []
    
    for i in range (0, len(ar)):
        outldl.append([i+1, my_seq[i], ar[i][0],ar[i][1],ar[i][2],
                       ar[i][3], ar[i][4], ar[i][5], ar[i][6], ar[i][7]])
    return outldl


def calculBulle3(ldl, fen):
    '''
    Search CG bubbles , define initial and final point on the sequence,
    the number of C max and the number of G min in a window,
    the difference CG% max.
    ----------
    Parameters
        ----------
        1 : list of list <type 'numpy.ndarray'>
        ex :[[nbA, nbT,nbG,nbC,%A,%T,%G,%C,%C-%G][...][...]] etc...
        ex2:
        [ 24.          15.          17.          22.          30.76923077
        19.23076923  21.79487179  28.20512821   6.41025641]
        [ 24.          15.          17.          22.          30.76923077
        19.23076923  21.79487179  28.20512821   6.41025641]
        [ 24.          14.          18.          22.          30.76923077
        17.94871795  23.07692308  28.20512821   5.12820513]
        ...longueur ldl= 63366
        [ 24.          17.          17.          20.          30.76923077
        21.79487179  21.79487179  25.64102564   3.84615385]
        [ 24.          17.          18.          19.          30.76923077
        21.79487179  23.07692308  24.35897436   1.28205128]
        [ 23.          17.          19.          19.          29.48717949
        21.79487179  24.35897436  24.35897436  -1.        ]
        2 : fen, integer (lenght windows 78)
    Returns
        -------
        return a b c d
        a = array <type 'numpy.ndarray'>
        b, c ,d = int (compteur)     
        ex1 a : [[num , Pi, pf,long Bulle, somme %C-%G, diffMax%C-%G, %Cmax, %Gmin, 
        %Cmoy, %Gmoy, %Amoy, %Tmoy][...]]
        ex2 a:
        [1, 1, 16, 16, 41.025, 5.12, 0.256, 0.2051, 0.24, 0.2155, 0.3221, 0.2]
        [2, 29, 30, 2, 1794871795, etc... , 0.26923076923076922]
        [3, 32, 38, 7, 6.4102564102564124,etc....]
        
    Example
        -------
    BulleCGbrPlusTri1, ctCGP1, ctCGP2, ctCGP3 = calculBulle2(TabBrinPlusTri1)
    Rq : diffMax%C-%G = (%Cmax - %Gmin) * 100,
    '''

    ar = ldl # voir si numpy modifie ldl ? sinon faire copie
    outldl = []  # [nb,m,pi,pf,(pf-pi),s]
    pi = 0  # position initial bulle
    pf = 0  # position final bulle
    s = 0  # Somme des nombres qui se suive
    nb = 0  # represente le numero de bulle CG
    compt = 0
    compt1 = 0  # le nombre de bulle total
    compt2 = 0  # nombre de == -1
    compt3 = 0  # nombre de != -1
    maxCG = 0  # difference C-G max ex: 30,2% - 25% = 4,8 %
    maxC = 0
    minG = 0
    nbrAt = 0   # nombre C total dans la bulle
    nbrTt = 0   # nombre G total dans la bulle
    nbrGt = 0   # nombre G total dans la bulle    
    nbrCt = 0   # nombre C total dans la bulle
    
    i = 0
    while i < len(ar):
        if (ar[i][8] == -1) and (compt == 0):
            # print "rien"
            compt2 += 1
        elif (ar[i][8] != -1) and (compt > 0):
            # print "suite bulle"
            compt3 += 1
            pf += 1
            s = s + ar[i][8]
            if (ar[i][8] > maxCG):
                maxCG = ar[i][8]
            if (ar[i][3] > maxC):
                maxC = ar[i][3]
            if (ar[i][2] < minG):
                minG = ar[i][2]
            nbrCt = nbrCt + ar[i][3]
            nbrGt = nbrGt + ar[i][2]
            nbrAt = nbrAt + ar[i][0]
            nbrTt = nbrTt + ar[i][1]    
            
        elif (ar[i][8] != -1) and (compt == 0):
            # print "debut de bulle"
            nb += 1
            compt += 1
            compt1 += 1
            compt3 += 1
            pi = (i + 1)
            pf = (i + 1)
            s = ar[i][8]
            maxCG = ar[i][8]
            maxC = ar[i][3]
            minG = ar[i][2]
            nbrCt = ar[i][3]
            nbrGt = ar[i][2]
            nbrAt = ar[i][0]
            nbrTt = ar[i][1] 
            
        elif (ar[i][8] == -1) and (compt > 0):
            # print "fin de bulle"
            compt = 0
            compt2 += 1
            lg = ((pf - pi) + 1)
            outldl.append([nb, pi, pf, lg, s,
                           maxCG, maxC/fen, minG/fen, 
                           (nbrCt/lg)/fen, (nbrGt/lg)/fen,
                            (nbrAt/lg)/fen, (nbrTt/lg)/fen])
            lg = 0
            pi = 0
            pf = 0
            s = 0
            maxCG = 0
            maxC = 0
            minG = 0
            nbrCt = 0
            nbrGt = 0
            nbrAt = 0
            nbrTt = 0             
            
        else:
            print("prob fonc calculBulle3")

        i += 1
#  si pas de bulle à la fin. Rq si bulle fin et debut de la seq...meme bulle
#  genome bacteriens circulaire :-))
    if (pi != 0):
        lg = ((pf - pi) + 1)
        outldl.append([nb, pi, pf, lg, s,
                       maxCG, maxC/fen, minG/fen,
                       (nbrCt/lg)/fen, (nbrGt/lg)/fen,
                        (nbrAt/lg)/fen, (nbrTt/lg)/fen])
    return outldl, compt1, compt2, compt3


def Filtres_bulles(ldl, lgbul):
    '''
    filtre ldl sur longueur de bulle et retourne la meme ldl (moins les bulle inferieur)
    '''

    ldlout =[] 
    myldl = ldl
#    print(len(myldl))
    for ldl in myldl:
#        print(ldl[0], ldl[1], ldl[2], ldl[3], ldl[4])
        ldlbu = ldl[5]
#        print(len(ldlbu))
        ldlbuout =[]
#        print(len(ldlbuout))
        for bu in ldlbu:
            #print(bu)
            longbu = bu[3]
            #print(bu[3])
            if longbu >= lgbul:
                ldlbuout.append(bu)
#        print(len(ldlbuout))        
        ldlout.append([ldl[0], ldl[1], ldl[2], ldl[3], ldl[4], ldlbuout, ldl[6]])

    return ldlout

def Re_ecrit_ldl_Bullemere(ldl, cut):
    '''
    reecrit ldl en filtrant les bulles et en sortant une ldl 1 ligne une bulle.
    '''
    ldlout = []    
    ldltmp = Filtres_bulles(ldl, cut)
    for lf in ldltmp:
        ldlBu = lf[5]
#        seq = lf[6]
        
        #print("ldlBu = ", ldlBu)
       
        for Bu in ldlBu:
            ldlbuout =[]
#            seq1 = ""
#            if len(seq) != 0:
#                seq1 = seq[(Bu[1] - 1):Bu[2]]
            ldlbuout.append(lf[0])
            ldlbuout.append(lf[1])
            ldlbuout.append(lf[2])
            ldlbuout.append(lf[3])
            ldlbuout.append(lf[4])
            for el in Bu:
                ldlbuout.append(el)
#            print(ldlbuout)
            ldlout.append(ldlbuout)
    
    return ldlout


def CalculDataBubble(ldl, cut):
    '''
    Calcul nombre de bulle en fonction de longueur etc.... 
    '''
    
    
    ldlout = []    
    ldltmp = Filtres_bulles(ldl, cut)
    for lf in ldltmp:
        ldlBubbles = lf[5]
#       nbrBu = len(ldlBubbles)
#        print("nombre de bu = ", nbrBu)
        
        SomLongBu = 0
        SomPercGC = 0
        nbrBulong1 = 0
        SomBuLong1 = 0
        SomPercGCL1 = 0
        nbrBulong10 = 0
        SomBuLong10 = 0
        SomPercGCL10 = 0
        nbrBulong30 = 0
        SomBuLong30 = 0
        SomPercGCL30 = 0
        nbrBulong50 = 0
        SomBuLong50 = 0
        SomPercGCL50 = 0
        
        plusLongbulle = 0
        numpluslongbulle = 0
        percplusLongbulle = 0
        
        pluspercGC = 0
        numpluspercGC = 0
        longpluspercGC = 0

        for bu in ldlBubbles:
#            print(bu)
            numbu = bu[0]
            longbu = bu[3]
            PGC = bu[4]
            SomLongBu = SomLongBu + longbu
            SomPercGC = SomPercGC + PGC
            if longbu >= 1:
                nbrBulong1 += 1
                SomBuLong1 = SomBuLong1 + longbu
                SomPercGCL1 = SomPercGCL1 + PGC                    
            if longbu >= 10:
                nbrBulong10 += 1
                SomBuLong10 = SomBuLong10 + longbu
                SomPercGCL10 = SomPercGCL10 + PGC
            if longbu >= 30:
                nbrBulong30 += 1
                SomBuLong30 = SomBuLong30 + longbu
                SomPercGCL30 = SomPercGCL30 + PGC
            if longbu >= 50:
                nbrBulong50 += 1
                SomBuLong50 = SomBuLong50 + longbu
                SomPercGCL50 = SomPercGCL50 + PGC
                
            if longbu > plusLongbulle:
                numpluslongbulle = numbu
                plusLongbulle = longbu
                percplusLongbulle = PGC
            
            if PGC > pluspercGC:
                numpluspercGC = numbu                
                pluspercGC = PGC
                longpluspercGC = longbu
                
            #print(SomPercGC, SomLongBu, nbrBulong1, nbrBulong10, nbrBulong30, nbrBulong50)
            
            
        ldlout.append([lf[0], lf[1], lf[2], lf[3], lf[4], SomLongBu, SomPercGC, \
            nbrBulong1, SomBuLong1, SomPercGCL1, \
            nbrBulong10, SomBuLong10, SomPercGCL10, \
            nbrBulong30, SomBuLong30, SomPercGCL30, \
            nbrBulong50, SomBuLong50, SomPercGCL50, \
            numpluslongbulle, plusLongbulle, percplusLongbulle, \
            numpluspercGC, longpluspercGC, pluspercGC, \
            lf[6]])
       
       
    
    return ldlout
    
    
def Re_ecrit_ldl_Bullemere_2(ldl):
    '''
    reecrit ldl en sortant une ldl 1 ligne une bulle.
    '''
    ldlout = []
    comptBu = 0
    comptNoBu = 0
    comptError = 0
    
    for lf in ldl:
        ldlBu = lf[5]
        nbrBu = len(ldlBu)
        seq = lf[6]
        longSeqMat = len(seq)
  #      print("lf = ", lf)
 #       print(seq)
#        print("lenghldlBu = ", nbrBu)
        
        if nbrBu != 0:
            comptBu += 1
#            print("traite-bulle")
            numBuseq = 0
            for Bu in ldlBu:
                numBuseq += 1
                corPiBu = Bu[1]
                corPfBu = Bu[2]
                longBu = (corPfBu - corPiBu) + 1
                corPiSeqBu = (corPiBu - 1)  # dut coordonnee python
                corPfSeqBu = (corPfBu + 77) # longeur de la fenetre 78
                seqBu = seq[corPiSeqBu:corPfSeqBu]
                longSeqBu = len(seqBu)
                disPiBuFinMat = longSeqMat - corPiBu
 #               print (disPiBuFinMat)
                ldlout.append([lf[0], lf[1], lf[2], lf[3], lf[4], seq, \
                    Bu[0], Bu[1], Bu[2], Bu[3], Bu[4], Bu[5], Bu[6], Bu[7], \
                    Bu[8], Bu[9], Bu[10], Bu[11], \
                    numBuseq, corPiBu, corPfBu, longBu, longSeqBu, disPiBuFinMat, seqBu])
                    
        elif nbrBu == 0:
            comptNoBu += 1
#            print("traite-particulier")
            ldlout.append([lf[0], lf[1], lf[2], lf[3], lf[4], seq, \
                0, 0, 0, 0, 0, 0, 0, 0, \
                0, 0, 0, 0, \
                0, 0, 0, 0, 0, 0, "-"])
        else:
            print("Error")
            comptError += 1
    
    
#    print("comptBu = ", comptBu)
#    print("comptNoBu = ", comptNoBu)
#    print("comptError = ", comptError)
        
    return ldlout    


def Transforme_col_ldl_en_dic_de_ldl (ma_liste_mere,num_lf_col):
    '''fonction
    entree une liste de liste (2D),un numero index liste fille qui donnerons la cle
    la valeur serra la liste de liste fille ayant la mem cle.
    sortie un dictionnaire 
    '''
    mon_dico = {}
    
    for i, liste_fille in enumerate(ma_liste_mere):
        ma_liste_val= []
        
        for j, elm_liste_fille in enumerate(liste_fille):
            if j == num_lf_col:
                ma_cle = elm_liste_fille
              
        mon_dico[ma_cle] = ma_liste_val

    for cle, valeur in mon_dico.items():
        ma_cle2 = cle
        #print ma_cle2
        ma_valeur2 = valeur
        #print ma_valeur2
        for k, liste_fille2 in enumerate(ma_liste_mere):
            #print liste_fille2
            if liste_fille2[num_lf_col] == ma_cle2:
                #print "oui"
                ma_valeur2.append(liste_fille2)    
              
        mon_dico[ma_cle2] = ma_valeur2
      
    return mon_dico


def Affiche_Dictionnaire (mon_dico):
    '''fonction
    entree un dictionnaire
    sortie ecran du dictionnaire 
    '''
    for cle, valeur in mon_dico.items():
        print ("La cle {} contient la valeur {}.".format(cle,valeur))


def Paramettre_Dictionnaire (mon_dico):

    long_dico = len(mon_dico)
    print (" la longueur du dictionnaire= " , long_dico)
    for cle, valeur in mon_dico.items():
        print ("pour la cle" , cle)
        print ("la longeur de la valeur =", len(valeur))

    
def Select_Mat_Bulle_Long(dico):
    '''fonction
    entree un dictionnaire
    sortie ecran du dictionnaire 
    '''
    ldlout =[]
    
    for cle, valeur in dico.items():
#       numMat = cle
        lfselect = []
        longselect = 0
        for lf in valeur:
            longtmp = lf[9]
            if longtmp >= longselect:
                lfselect = lf
                longselect = longtmp
#        print(numMat, longselect, lfselect)
        ldlout.append(lfselect)
    
    return ldlout
            

def Recherche_RegExp_Seq_SelectBulle(ldl):
    
    ldlout = []
    
    for lf in ldl:
        numseq = lf[0]
        nameseq = lf[1]
        numBuseq = lf[6]
        nameBuseq = nameseq + "_" + str(numBuseq)
        seqBu = lf[24]
#        print(numseq, nameBuseq, seqBu)
        ldlseqregexp = Recherche_Multi_RegExp_Bu_1(numseq, nameBuseq,seqBu)
        ldlout.append(ldlseqregexp)
    
    return ldlout



def Recherche_Multi_RegExp_Bu_1(compt0, name, seqstring):
    
    TabRegExp =["(T|C)C", \
    "(T|C)C[ATGC]{9,13}(T|C)C", \
    "(T|C)C[ATGC]{9,13}(T|C)C[ATGC]{9,13}(T|C)C", \
    "(T|C)C[ATGC]{9,13}(T|C)C[ATGC]{9,13}(T|C)C[ATGC]{9,13}(T|C)C", \
    "(T|C)C[ATGC]{9,13}(T|C)C[ATGC]{9,13}(T|C)C[ATGC]{9,13}(T|C)C[ATGC]{9,13}(T|C)C", \
    "(T|C)C[ATGC]{9,13}(T|C)C[ATGC]{9,13}(T|C)C[ATGC]{9,13}(T|C)C[ATGC]{9,13}(T|C)C[ATGC]{9,13}(T|C)C"]
    
    longSeq =len(seqstring)
    comptre = 0
    myldl = [compt0, name, seqstring]
    for re in TabRegExp:
        myldlseq = []
        myseqcount= []
        mydensi = 0
        comptre += 1 
        myldlseq, myseqcount = Recherche_RegExp(re, seqstring)
        mydensi = (float(myseqcount)/float(longSeq)) * 1000
        Rmydensi = round(mydensi, 2)
        #print(compt0, name, comptre, re, myseqcount, Rmydensi)
        myldl.append([comptre, re, myseqcount, Rmydensi])
    
    return myldl    
    

def Re_ecrit_ldl_2(ldl):
    '''    
    writes a LDL simple utilisable for founc ecrit_ldl()
        Parameters
        ----------
        a : list of list (ldl complexe)
        
        Returns
        -------
        a : list of list (ldl simple)
    '''
    
    outldl = []
    
    for ll in ldl:
        ldlfille = [ll[0], ll[1], \
            ll[3][2], ll[3][3], ll[4][2], ll[4][3], \
            ll[5][2], ll[5][3], ll[6][2], ll[6][3], \
            ll[7][2], ll[7][3], ll[8][2], ll[8][3], \
            ll[2]]
            
        outldl.append(ldlfille)
    
    return outldl


def Fait_ldl_Fasta_1(ldl):
    
    outldl =[]
    
    for ll in ldl:
        namefasta = ">numseq_" + str(ll[0]) + "_" + ll[1]
        seqBu = ll[14]
        outldl.append([namefasta])
        outldl.append([seqBu])

    return outldl


def Fait_ldl_Fasta_2(ldl):
    
    outldl =[]
    
    for ll in ldl:    
        namefasta = ">numseq_" + str(ll[0]) + "_" + ll[1]
        seqBu = ll[14]
        if seqBu != "-":
            outldl.append([namefasta])
            outldl.append([seqBu])

    return outldl



def Fait_Tab_Var_1(DB, DR, DC, DBML, DRML, DCT):
    '''
    Make LDL whit var for tab Pls-Da
    DB = Data Bubbles 
    DR = Datas RegExp
    DCT = Datas Compts seq Matrices
    DBML = Datas Bubbles Max Lengh
    DRML = Datas Reg Bubbles Max Lengh
    DCTT = Datas Compt Triplet Seq    
    '''
    ldlout=[]
    
    i = 0
    while i < len(DB):
        ldlout.append([DB[i][0], DB[i][1], DB[i][2], DB[i][3], DB[i][4],\
            DB[i][5], DB[i][6], DB[i][7], DB[i][8], DB[i][9],\
            ((float(DB[i][8])/DB[i][2]) * 1000), ((DB[i][9]/DB[i][2]) * 1000),\
            DB[i][10], DB[i][11], DB[i][12], DB[i][13], DB[i][14],\
            DB[i][15], DB[i][16], DB[i][17], DB[i][18], DB[i][19],\
            DB[i][20], DB[i][21], DB[i][22], DB[i][23], DB[i][24], DB[i][25],\
#
            DR[i][0], DR[i][1], DR[i][2], DR[i][3],DR[i][4],\
            DR[i][5], DR[i][6], DR[i][7], DR[i][8],DR[i][9],\
            DR[i][10], DR[i][11], DR[i][12], DR[i][13],DR[i][14],\
#           DR[i][0], DR[i][1], DR[i][2], DR[i][3],DR[i][4],\
#
            DC[i][0], DC[i][1], DC[i][2], DC[i][3], DC[i][4],\
            DC[i][5], DC[i][6], DC[i][7], DC[i][8], DC[i][9],\
            DC[i][10], DC[i][11], DC[i][12], DC[i][13], DC[i][14],\
            DC[i][15], DC[i][16], DC[i][17], DC[i][18], DC[i][19],\
            DC[i][20], DC[i][21], DC[i][22], DC[i][23], DC[i][24],\
            DC[i][25], DC[i][26], DC[i][27], DC[i][28], DC[i][29],\
            DC[i][30], DC[i][31], DC[i][32], DC[i][33], DC[i][34],\
            DC[i][35], DC[i][36], DC[i][37], DC[i][38], DC[i][39],\
            DC[i][40], DC[i][41], DC[i][42], DC[i][43], DC[i][44],\
            DC[i][45], DC[i][46], DC[i][47], DC[i][48], DC[i][49],\
#
            DBML[i][0], DBML[i][1], DBML[i][2], DBML[i][3], DBML[i][4],\
            DBML[i][5], DBML[i][6], DBML[i][7], DBML[i][8], DBML[i][9],\
            DBML[i][10], DBML[i][11], DBML[i][12], DBML[i][13],((DBML[i][14] - DBML[i][15]) * 100), DBML[i][14],\
            DBML[i][15], DBML[i][16], DBML[i][17], DBML[i][18], DBML[i][19],\
            DBML[i][20], DBML[i][21], DBML[i][22], DBML[i][23], DBML[i][24],\
            
            DRML[i][0], DRML[i][1], DRML[i][2], DRML[i][3], DRML[i][4],\
            DRML[i][5], DRML[i][6], DRML[i][7], DRML[i][8], DRML[i][9],\
            DRML[i][10], DRML[i][11], DRML[i][12], DRML[i][13], DRML[i][14],\
            
#
            DCT[i][0], DCT[i][1], DCT[i][2], DCT[i][3], DCT[i][4],\
            DCT[i][5], DCT[i][6], DCT[i][7], DCT[i][8], DCT[i][9],\
            DCT[i][10], DCT[i][11], DCT[i][12], DCT[i][13], DCT[i][14],\
            DCT[i][15], DCT[i][16], DCT[i][17], DCT[i][18], DCT[i][19],\
            DCT[i][20], DCT[i][21], DCT[i][22], DCT[i][23], DCT[i][24],\
            DCT[i][25], DCT[i][26], DCT[i][27], DCT[i][28], DCT[i][29],\
            DCT[i][30], DCT[i][31], DCT[i][32], DCT[i][33], DCT[i][34],\
            DCT[i][35], DCT[i][36], DCT[i][37], DCT[i][38], DCT[i][39],\
            DCT[i][40], DCT[i][41], DCT[i][42], DCT[i][43], DCT[i][44],\
            DCT[i][45], DCT[i][46], DCT[i][47], DCT[i][48], DCT[i][49],\
            DCT[i][50], DCT[i][51], DCT[i][52], DCT[i][53], DCT[i][54],\
            DCT[i][55], DCT[i][56], DCT[i][57], DCT[i][58], DCT[i][59],\
            DCT[i][60], DCT[i][61], DCT[i][62], DCT[i][63], DCT[i][64],\
            DCT[i][65], DCT[i][66], DCT[i][67], DCT[i][68], DCT[i][69],\
            DCT[i][70], DCT[i][71], DCT[i][72], DCT[i][73], DCT[i][74],\
            DCT[i][75], DCT[i][76], DCT[i][77], DCT[i][78], DCT[i][79],\
            DCT[i][80], DCT[i][81], DCT[i][82], DCT[i][83], DCT[i][84],\
            DCT[i][85], DCT[i][86], DCT[i][87], DCT[i][88], DCT[i][89],\
            DCT[i][90], DCT[i][91], DCT[i][92], DCT[i][93], DCT[i][94],\
            DCT[i][95], DCT[i][96], DCT[i][97], DCT[i][98], DCT[i][99],\
            DCT[i][100], DCT[i][101], DCT[i][102], DCT[i][103], DCT[i][104],\
            DCT[i][105], DCT[i][106], DCT[i][107], DCT[i][108], DCT[i][109],\
            DCT[i][110], DCT[i][111], DCT[i][112], DCT[i][113], DCT[i][114],\
            DCT[i][115], DCT[i][116], DCT[i][117], DCT[i][118], DCT[i][119],\
            DCT[i][120], DCT[i][121], DCT[i][122], DCT[i][123], DCT[i][124],\
            DCT[i][125], DCT[i][126], DCT[i][127], DCT[i][128], DCT[i][129],\
            DCT[i][130], DCT[i][131],
            
            ])
            
        i += 1
    
    return ldlout
    
    
def ecrit_ldl_tabtot(ldl, nom_fichier_out):
    '''
    writes a column LDL in a separate fGCile with tab
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
    fichier_out.write("num-seq\tnames-seq\tlongseq\tlongwin\tnbrBubble\tSomLongBu\tSomSurfCG\t")
    fichier_out.write("nbrBulong1\tSomBuLong1\tSomSurfCGL1\tDenLen1Tem\tDenSurf1Tem\t")
    fichier_out.write("nbrBulong10\tSomBuLong10\tSomSurfCGL10\t")
    fichier_out.write("nbrBulong30\tSomBuLong30\tSomSurfCGL30\t")
    fichier_out.write("nbrBulong50\tSomBuLong50\tSomSurfCGL50\t")
    fichier_out.write("NumBuPlusLongBulle\tLongPlusLongBulle\tSurfPlusLongBulle\t")
    fichier_out.write("NumBuPlusSurfCG\tLongPlusSurfCG\tSurfPlusSurfCG\t")
    fichier_out.write("Seq\t")
    
    fichier_out.write("num-seq\tnames-seq\t")
    fichier_out.write("Nb_Reg1\tDen_Reg1\tNb_Reg2\tDen_Reg2\t")
    fichier_out.write("Nb_Reg3\tDen_Reg3\tNb_Reg4\tDen_Reg4\t")
    fichier_out.write("Nb_Reg5\tDen_Reg5\tNb_Reg6\tDen_Reg6\t")
    fichier_out.write("Seq\t")

    fichier_out.write("num-seq\tnames-seq\tlongseq\tnbrA\tnbrT\tnbrG\tnbrC\t")
    fichier_out.write("PcA\tPcT\tPcG\tPcC\t")
    fichier_out.write("PcApT\tPcApG\tPcApC\tPcTpG\tPcTpC\tPcGpC\t")
    fichier_out.write("nbrAA\tnbrAT\tnbrAG\tnbrAC\t")
    fichier_out.write("nbrTA\tnbrTT\tnbrTG\tnbrTC\t")
    fichier_out.write("nbrGA\tnbrGT\tnbrGG\tnbrGC\t")
    fichier_out.write("nbrCA\tnbrCT\tnbrCG\tnbrCC\t")
    fichier_out.write("PcAA\tPcAT\tPcAG\tPcAC\t")
    fichier_out.write("PcTA\tPcTT\tPcTG\tPcTC\t")
    fichier_out.write("PcGA\tPcGT\tPcGG\tPcGC\t")
    fichier_out.write("PcCA\tPcCT\tPcCG\tPcCC\t")
    fichier_out.write("Seq\t")   
    
    fichier_out.write("num-seq\tnames-seq\tlongseq\tlongwin\tnbrBubble\tseqMat\tnumBuseq\t")
    fichier_out.write("PiBu\tPfBu\tLenBub\tSurfCGBu\tdiffCGMax%Bu\t%CMaxBu\t%GMinBu\tdiffCGmoy%\t%CMoyBu\t%GMoyBu\t%AMoyBu\t%TMoyBu\t")
    fichier_out.write("numBuseq\tPiBu\tPfBu\tLenghBu\tLenghSeqBu\tDisPiBuFinMat\tSeqBu\t")    
    
    fichier_out.write("num-seq\tnames-seq\t")
    fichier_out.write("Nb_Reg1Bu\tDen_Reg1Bu\tNb_Reg2Bu\tDen_Reg2Bu\t")
    fichier_out.write("Nb_Reg3Bu\tDen_Reg3Bu\tNb_Reg4Bu\tDen_Reg4Bu\t")
    fichier_out.write("Nb_Reg5Bu\tDen_Reg5Bu\tNb_Reg6Bu\tDen_Reg6Bu\t")
    fichier_out.write("Seq\t")
   
    fichier_out.write("num-seq\tnames-seq\tlongseq\t")
    for cod1 in ListCodons:
        nameCodon = ("nbr" + cod1)
        fichier_out.write((nameCodon + "\t"))
    for cod2 in ListCodons:
        PercCodon = ("Pc" + cod2)
        fichier_out.write((PercCodon + "\t"))
    fichier_out.write("Seq\n")   
   
   
    fichier_out.write(('\n'.join('\t'.join(str(x) for x in l) for l in ldl)))
    fichier_out.write("\n")
    fichier_out.close()
    print("end write")        


def Fait_Tab_Var_Mod1(DB, DR, DCT, DBML, DRML):
    '''
    Make LDL whit var for tab Pls-Da
    DB = Data Bubbles 
    DR = Datas RegExp
    DCT = Datas Compts seq Matrices
    DBML = Datas Bubbles Max Lengh
    DRML = Datas Bubbles Max Lengh    
    '''
    ldlout=[]
    
    i = 0
    while i < len(DB):
        ldlout.append([DB[i][0], DB[i][1], DB[i][8], DB[i][9],\
            DB[i][20], DB[i][21],\
#
            DR[i][0], DR[i][1], DR[i][3], DR[i][7],DR[i][9],\
#
            DCT[i][0], DCT[i][1], DCT[i][7], DCT[i][8], DCT[i][9],\
            DCT[i][10], DCT[i][12], DCT[i][15], DCT[i][43], DCT[i][48],\
#
            DBML[i][0], DBML[i][1], DBML[i][10], DBML[i][11],\
            
            DRML[i][0], DRML[i][1], DRML[i][2], DRML[i][4]])
            
        i += 1
    
    return ldlout


def Recherche_Codons_1(seq, numseq, nameseq, List_C):
    
    '''
   researche triplets in sequence c est un count que l on fait.
        Parameters
        ----------
        a : sequence (type biopython)
        b : num seq
        c : name seq
        -------
       returne l avec numseq, names, nombre par codons, seq
        
    '''
    
#    myseq = str(seq) seq type bio python
    
    ldlout=[]
    ldl_C_nbrC=[]
    
    namecodon = "XXX"
    nbrcodon = 0
    
    for Cod in List_C:
        namecodon = "nbr" + Cod
        nbrcodon = seq.count(Cod)
        ldl_C_nbrC.append([namecodon, nbrcodon])
    
    ldlout= [numseq, nameseq, len(seq), ldl_C_nbrC, str(seq)]
    
    return ldlout


def Reecrit_et_calcul__ldl_codonsmere(ldl):
    '''    
    writes a LDL simple utilisable for founc ecrit_ldl()
        Parameters
        ----------
        a : list of list (ldl complexe)
        
        Returns
        -------
        a : list of list (ldl simple)
    '''
#PcA =round(((100 * float(nbrA))/len(my_seq)), 2)    
    outldl = []
    
    for ll in ldl:
        
        lout =[]        
        numseq = ll[0]
        nameseq = ll[1]
        lenseq=int(ll[2])
        ldlCodNbr=ll[3]
        seq = ll[4]
        
        lout =[numseq, nameseq, lenseq]
        
        i = 0
# reecrit le nombre (dans l ordre de la list des codons initial)        
        while i < len(ldlCodNbr):
#            print(i, ldlCodNbr[i][0], ldlCodNbr[i][1])
            nbrCodi = int(ldlCodNbr[i][1])
            lout.append(nbrCodi)
            i += 1

# calcul et ecrit la densite pour 1000 bases)

        k = 0
# reecrit le Pourcentage en fonction de la longueur
#        (dans l ordre de la list des codons initial)        
        while k < len(ldlCodNbr):
#            print(i, ldlCodNbr[i][0], ldlCodNbr[i][1])
            nbrCodk = float(ldlCodNbr[k][1])
            PercentNbrCod = round(((nbrCodk / lenseq) * 100), 2)
            lout.append(PercentNbrCod)
            k += 1
        
        lout.append(seq)
        
        outldl.append(lout)
    
    return outldl




def Re_ecrit_ldl_Bullemere_3(ldl):
    '''
    reecrit ldl  bulles meres (eleimine bulle de C=G)
    '''
    ldlout = []    
    for lf in ldl:
        numseq = lf[0]
        nameseq = lf[1]
        lenSeq = lf[2]
        lenWin = lf[3]
        nbBu = 0
        ldlBu = lf[5]
#        print(ldlBu)
        ADN = lf[6]
        ldlBuout = []
        if len(ldlBu) == 0:
            ldlout.append([numseq, nameseq,lenSeq,lenWin,nbBu,ldlBuout,ADN])
        elif len(ldlBu) != 0:
            numBu = 0
            for Bu in ldlBu:
               
                piBu = Bu[1]
                pfBu = Bu[2]
                lenBu = Bu[3]
                SomBu = Bu[4]
                diffMax = Bu[5]
                Cmax = Bu[6]
                Gmin = Bu[7]
                Cmoy = Bu[8]
                Gmoy = Bu[9]
                Amoy = Bu[10]
                Tmoy = Bu[11]
                
                if diffMax != 0 :
                    numBu += 1
                    nbBu += 1
                    ldlBuout.append([numBu, piBu, pfBu, lenBu, SomBu, diffMax, Cmax, Gmin,Cmoy,Gmoy, Amoy, Tmoy])
                    
            ldlout.append([numseq, nameseq,lenSeq,lenWin,nbBu,ldlBuout,ADN])
        else:
            print("Erreur dans Re_ecrit_ldl_Bullemere_3")
            
    return ldlout


def Fait_Tab_Var_OPLSDA(DB, DR, DC, DBML, DRML, DCT):
    '''
    Make LDL whit var for tab Pls-Da
    DB = Data Bubbles 
    DR = Datas RegExp
    DCT = Datas Compts seq Matrices
    DBML = Datas Bubbles Max Lengh
    DRML = Datas Reg Bubbles Max Lengh
    DCTT = Datas Compt Triplet Seq    
    '''
    ldlout=[]
    
    i = 0
    while i < len(DB):
        ldlout.append([DB[i][0], DB[i][1], DB[i][2], DB[i][3], DB[i][4], 'Empty',\
            DB[i][8], DB[i][9],\
            ((float(DB[i][8])/DB[i][2]) * 1000), ((DB[i][9]/DB[i][2]) * 1000),\
            DB[i][20], DB[i][21],\
#
            DR[i][3],\
            DR[i][5], DR[i][7],DR[i][9],\
            DR[i][11], DR[i][13],\
#           DR[i][0], DR[i][1], DR[i][2], DR[i][3],DR[i][4],\
#
#            DC[i][0], DC[i][1], DC[i][2], DC[i][3], DC[i][4],\
            DC[i][7], DC[i][8], DC[i][9],\
            DC[i][10], DC[i][33], DC[i][34],\
            DC[i][35], DC[i][36], DC[i][37], DC[i][38], DC[i][39],\
            DC[i][40], DC[i][41], DC[i][42], DC[i][43], DC[i][44],\
            DC[i][45], DC[i][46], DC[i][47], DC[i][48],\
#
#            DBML[i][0], DBML[i][1], DBML[i][2], DBML[i][3], DBML[i][4],\
#            DBML[i][5], DBML[i][6], DBML[i][7], DBML[i][8], DBML[i][9],\
            DBML[i][11], DBML[i][12], DBML[i][13],((DBML[i][14] - DBML[i][15]) * 100), DBML[i][14],\
            DBML[i][15], DBML[i][16], DBML[i][17],\
#            DBML[i][20], DBML[i][21], DBML[i][22], DBML[i][23], DBML[i][24],\
            
            DRML[i][2],  DRML[i][4],\
            DRML[i][6],  DRML[i][8],\
            DRML[i][10], DRML[i][12],\
            
#
#            DCT[i][0], DCT[i][1], DCT[i][2], DCT[i][3], DCT[i][4],\
#            DCT[i][5], DCT[i][6], DCT[i][7], DCT[i][8], DCT[i][9],\
#            DCT[i][10], DCT[i][11], DCT[i][12], DCT[i][13], DCT[i][14],\
#           DCT[i][15], DCT[i][16], DCT[i][17], DCT[i][18], DCT[i][19],\
#            DCT[i][20], DCT[i][21], DCT[i][22], DCT[i][23], DCT[i][24],\
#            DCT[i][25], DCT[i][26], DCT[i][27], DCT[i][28], DCT[i][29],\
#            DCT[i][30], DCT[i][31], DCT[i][32], DCT[i][33], DCT[i][34],\
#            DCT[i][35], DCT[i][36], DCT[i][37], DCT[i][38], DCT[i][39],\
#            DCT[i][40], DCT[i][41], DCT[i][42], DCT[i][43], DCT[i][44],\
#            DCT[i][45], DCT[i][46], DCT[i][47], DCT[i][48], DCT[i][49],\
#            DCT[i][50], DCT[i][51], DCT[i][52], DCT[i][53], DCT[i][54],\
#            DCT[i][55], DCT[i][56], DCT[i][57], DCT[i][58], DCT[i][59],\
#            DCT[i][60], DCT[i][61], DCT[i][62], DCT[i][63], DCT[i][64],\
#            DCT[i][65], DCT[i][66],
            DCT[i][67], DCT[i][68], DCT[i][69],\
            DCT[i][70], DCT[i][71], DCT[i][72], DCT[i][73], DCT[i][74],\
            DCT[i][75], DCT[i][76], DCT[i][77], DCT[i][78], DCT[i][79],\
            DCT[i][80], DCT[i][81], DCT[i][82], DCT[i][83], DCT[i][84],\
            DCT[i][85], DCT[i][86], DCT[i][87], DCT[i][88], DCT[i][89],\
            DCT[i][90], DCT[i][91], DCT[i][92], DCT[i][93], DCT[i][94],\
            DCT[i][95], DCT[i][96], DCT[i][97], DCT[i][98], DCT[i][99],\
            DCT[i][100], DCT[i][101], DCT[i][102], DCT[i][103], DCT[i][104],\
            DCT[i][105], DCT[i][106], DCT[i][107], DCT[i][108], DCT[i][109],\
            DCT[i][110], DCT[i][111], DCT[i][112], DCT[i][113], DCT[i][114],\
            DCT[i][115], DCT[i][116], DCT[i][117], DCT[i][118], DCT[i][119],\
            DCT[i][120], DCT[i][121], DCT[i][122], DCT[i][123], DCT[i][124],\
            DCT[i][125], DCT[i][126], DCT[i][127], DCT[i][128], DCT[i][129],\
            DCT[i][130],
            
            ])
            
        i += 1
    
    return ldlout
    
    
def ecrit_ldl_tabOPLSDA(ldl, nom_fichier_out):
    '''
    writes a column LDL in a separate fGCile with tab
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
    fichier_out.write("num-seq\tnames-seq\tlongseq\tlongwin\tnbrBubble\tTempDenDG\t")
    fichier_out.write("SomBuLong1\tSomSurfCGL1\tDenLen1Tem\tDenSurf1Tem\t")
#    fichier_out.write("nbrBulong10\tSomBuLong10\tSomSurfCGL10\t")
#    fichier_out.write("nbrBulong30\tSomBuLong30\tSomSurfCGL30\t")
#    fichier_out.write("nbrBulong50\tSomBuLong50\tSomSurfCGL50\t")
    fichier_out.write("LongPlusLongBulle\tSurfPlusLongBulle\t")
#    fichier_out.write("NumBuPlusSurfCG\tLongPlusSurfCG\tSurfPlusSurfCG\t")
#    fichier_out.write("Seq\t")
    
#    fichier_out.write("num-seq\tnames-seq\t")
    fichier_out.write("Den_Reg1\tDen_Reg2\t")
    fichier_out.write("Den_Reg3\tDen_Reg4\t")
    fichier_out.write("Den_Reg5\tDen_Reg6\t")
#    fichier_out.write("Seq\t")

#    fichier_out.write("num-seq\tnames-seq\tlongseq\tnbrA\tnbrT\tnbrG\tnbrC\t")
    fichier_out.write("PcA\tPcT\tPcG\tPcC\t")
#    fichier_out.write("PcApT\tPcApG\tPcApC\tPcTpG\tPcTpC\tPcGpC\t")
#    fichier_out.write("nbrAA\tnbrAT\tnbrAG\tnbrAC\t")
#    fichier_out.write("nbrTA\tnbrTT\tnbrTG\tnbrTC\t")
#    fichier_out.write("nbrGA\tnbrGT\tnbrGG\tnbrGC\t")
#    fichier_out.write("nbrCA\tnbrCT\tnbrCG\tnbrCC\t")
    fichier_out.write("PcAA\tPcAT\tPcAG\tPcAC\t")
    fichier_out.write("PcTA\tPcTT\tPcTG\tPcTC\t")
    fichier_out.write("PcGA\tPcGT\tPcGG\tPcGC\t")
    fichier_out.write("PcCA\tPcCT\tPcCG\tPcCC\t")
#    fichier_out.write("Seq\t")   
    
#    fichier_out.write("num-seq\tnames-seq\tlongseq\tlongwin\tnbrBubble\tseqMat\tnumBuseq\t")
    fichier_out.write("diffCGMax%Bu\t%CMaxBu\t%GMinBu\tdiffCGmoy%\t%CMoyBu\t%GMoyBu\t%AMoyBu\t%TMoyBu\t")
#    fichier_out.write("numBuseq\tPiBu\tPfBu\tLenghBu\tLenghSeqBu\tDisPiBuFinMat\tSeqBu\t")    
    
#    fichier_out.write("num-seq\tnames-seq\t")
    fichier_out.write("Nb_Reg1Bu\tNb_Reg2Bu\t")
    fichier_out.write("Nb_Reg3Bu\tNb_Reg4Bu\t")
    fichier_out.write("Nb_Reg5Bu\tNb_Reg6Bu\t")
#    fichier_out.write("Seq\t")
   
#    fichier_out.write("num-seq\tnames-seq\tlongseq\t")
#    for cod1 in ListCodons:
#        nameCodon = ("nbr" + cod1)
#        fichier_out.write((nameCodon + "\t"))
    for cod2 in ListCodons:
        PercCodon = ("Pc" + cod2)
        fichier_out.write((PercCodon + "\t"))
    fichier_out.write("\n")   
   
   
    fichier_out.write(('\n'.join('\t'.join(str(x) for x in l) for l in ldl)))
    fichier_out.write("\n")
    fichier_out.close()
    print("end write")        


        
"""
Scripts for detection of C>G bubbles and calculation of sequence descriptors
@author: EE
"""
debut1 = time.time()

# checks if only 2 argument was passed to call the script
if len(sys.argv) !=3 :
    print("Usage: python prog-python /PATH/file obs")
    quit()
else:
    print("\n-- Welcome --\n")

print('----------------input_output_file-------------------')

chemin = os.getcwd()
argument = sys.argv
infichier1 = chemin + "/" + argument[1]
obs = argument[2]
# elimine le .gbk
#nom_out = argument[2][0:-4]
nomsplit1 = deconcat(argument[1])
nomfasta = nomsplit1[-1][0:-6]
print("deconcat = :", nomsplit1) 
print("nomfasta = ", nomfasta)

nom_out_fichier15 = nomfasta + "_" + obs + "_VarOPLSDA.txt"
outfichier15 = chemin + "/" + nom_out_fichier15


print("Directory =", chemin)
print ("arguments =", argument)
print ("infichier1 =", infichier1)
print ("observation=", obs)

print ("nom_out_fichier =", nom_out_fichier15)
print ("outfichier =", outfichier15)

ldlinfos =[]
inf01 = [chemin, argument, infichier1, obs]
ldlinfos.append(inf01)

compt0 = 0
ldlcomptemere = []
ldlregexpmere = []
ldlbasemere = []
ldlbullemere =[]
ldlcodonsmere =[]

ListCodons = ['TTT', 'TTC', 'TTA', 'TTG', 'CTT', 
  'CTC', 'CTA', 'CTG', 'ATT', 'ATC', 
  'ATA', 'ATG', 'GTT', 'GTC', 'GTA', 
  'GTG', 'TAT', 'TAC', 'TAA', 'TAG', 
  'CAT', 'CAC', 'CAA', 'CAG', 'AAT', 
  'AAC', 'AAA', 'AAG', 'GAT', 'GAC', 
  'GAA', 'GAG', 'TCT', 'TCC', 'TCA', 
  'TCG', 'CCT', 'CCC', 'CCA', 'CCG', 
  'ACT', 'ACC', 'ACA', 'ACG', 'GCT', 
  'GCC', 'GCA', 'GCG', 'TGT', 'TGC', 
  'TGA', 'TGG', 'CGT', 'CGC', 'CGA', 
  'CGG', 'AGT', 'AGC', 'AGA', 'AGG', 
  'GGT', 'GGC', 'GGA', 'GGG'] 


print("----------------Seq Fasta-------------------")
# recupere sequence ADN fichier Genbank

for seq_record in SeqIO.parse(infichier1, "fasta"):
    
#    print("name =", seq_record.id)
#    print("seq_record.seq = ", seq_record.seq)
#    print("seq =",(repr(seq_record.seq)))
    myseq = seq_record.seq
    compt0 = compt0 + 1    
    name = seq_record.id
    seqstring = str(myseq)
    #print("seq num %s = seq num %i" %(name,compt0)) 
    #print(type(myseq))
############ Compt ###############

    ldlcomptefille = Comptes_Bases_Seq(myseq, compt0, name) 
    ldlcomptemere.append(ldlcomptefille[0])

############ RegEXp #############
    
    ldlregexpfille = Recherche_Multi_RegExp(compt0, name, seqstring)
    #print("ldlregexpfille =", ldlregexpfille)
    ldlregexpmere.append(ldlregexpfille)


############## Bubble #################

    window1 = 78
    
    ldlbullefille , ldlbasefille, NbBulle, Nbneg, NbDifneg = Recherche_Multi_Bulle (myseq, window1)     
    # ldl bullefille = num-Bubble\tpi\tpf\tLenBub\tSumCG%\tdiffCGMax%\t%CMax\t%GMin\t%CMoy\t%GMoy\t%AMoy\t%TMoy\n
    ldlbasemere.append([compt0, name, len(myseq), ldlbasefille])
    #print_tronc(ldlbasefille, 5)
    #print(ldlbullefille[0], NbBulle, Nbneg, NbDifneg)
    ldlbullemere.append([compt0, name, len(myseq), (Nbneg + NbDifneg), NbBulle, ldlbullefille, seqstring])


############ codons ################"

    ldlcodonsfille = Recherche_Codons_1(myseq, compt0, name, ListCodons)
#    print(ldlcodonsfille)
    ldlcodonsmere.append(ldlcodonsfille)

###############
#   fin de boucle for #
############### 


ldlCodons = Reecrit_et_calcul__ldl_codonsmere(ldlcodonsmere)

#print_tronc(ldlCodons,4)


############prepare pour ecrire et traite ##############

ldlRegEXP = Re_ecrit_ldl_1(ldlregexpmere)

################### Bulle mere ########################

#print_tronc(ldlbullemere, 5)
print("----------------------------------------------------------------")
ldlbullemere2 = Re_ecrit_ldl_Bullemere_3(ldlbullemere)
#print_tronc(ldlbullemere2, 5)


############ filtre Bubble ##########

longbul = 1 # taille de la bulle en fenetre (prob 0.056)

ldlbullemerefiltre1 = Filtres_bulles(ldlbullemere2, longbul)

########## Fait tableau Paramettre 1 ligne une Bulle ##################

cutoff = 1

ldlParametreBulle = Re_ecrit_ldl_Bullemere(ldlbullemere2, cutoff)
#print_tronc(ldlParametreBulle, 3)

################# Tableau somme Data bulle par sequence #############"

ldlDataBubble = CalculDataBubble(ldlbullemere2, cutoff)

################# Tableau cherche RegExp pour la plus grande Bulle #############"
print("----------16---------")
#print_tronc(ldlbullemere2,3)
ldlParametreBulle2 = Re_ecrit_ldl_Bullemere_2(ldlbullemere2)

########### filtre une Bulle par Matrice ################
print("----------17--------")
DicParametreBulle2 = Transforme_col_ldl_en_dic_de_ldl(ldlParametreBulle2, 0)
#Paramettre_Dictionnaire(DicParametreBulle2)
# Affiche_Dictionnaire(DicParametreBulle2)

ldlMatBulleLaPlusLongue = Select_Mat_Bulle_Long(DicParametreBulle2)
#print_tronc(ldlMatBulleLaPlusLongue, 20)

########### recherche regexp dans bulle la pluslongue ################
print("---------18----------")
ldlMatRegExpBulleLaPlusLongue = Recherche_RegExp_Seq_SelectBulle(ldlMatBulleLaPlusLongue)
#print_tronc(ldlMatRegExpBulleLaPlusLongue, 5)
ldlRegExpBuLaplusLongue = Re_ecrit_ldl_2(ldlMatRegExpBulleLaPlusLongue)
#print_tronc(ldlRegExpBuLaplusLongue, 5)
print("---------19-----------")
########### Fasta seq bulle la plus longue ################

ldlFastaMatBuLaplusLongue = Fait_ldl_Fasta_1(ldlRegExpBuLaplusLongue)
#print_tronc(ldlFastaMatBuLaplusLongue, 6)


ldlFastaBuLaplusLongue = Fait_ldl_Fasta_2(ldlRegExpBuLaplusLongue)
#print_tronc(ldlFastaBuLaplusLongue, 6)

print("-----------20----------")
########### Fait Tableau Variable PlsDa ################

ldlTabVarPlsda1 = Fait_Tab_Var_1(ldlDataBubble, ldlRegEXP, ldlcomptemere, \
                        ldlMatBulleLaPlusLongue, ldlRegExpBuLaplusLongue, ldlCodons)
#print_tronc(ldlTabVarPlsda1, 5)

ldlTabVarPlsMod1 = Fait_Tab_Var_Mod1(ldlDataBubble, ldlRegEXP, ldlcomptemere, \
                        ldlMatBulleLaPlusLongue, ldlRegExpBuLaplusLongue)

ldlTabVarOPLSDA = Fait_Tab_Var_OPLSDA(ldlDataBubble, ldlRegEXP, ldlcomptemere, \
                        ldlMatBulleLaPlusLongue, ldlRegExpBuLaplusLongue, ldlCodons)



#print_tronc(ldlTabVarPlsMod1, 5)


########### Infos ################


inf02 = ["nombre de seq analyse pour comptage =  ", compt0]
ldlinfos.append(inf02)

inf03 = ["longueur ldlcomptemere = ", len(ldlcomptemere)]
ldlinfos.append(inf03)

inf04 = ["longueur ldlregexmere = ", len(ldlregexpmere)]
ldlinfos.append(inf04)

inf05 = ["longueur ldlRegEXP = ", len(ldlRegEXP)]
ldlinfos.append(inf05)

inf06 = ["longueur ldlbasemere = ", len(ldlbasemere)]
ldlinfos.append(inf06)

inf07 = ["longueur ldlbullemere = ", len(ldlbullemere)]
ldlinfos.append(inf07)

inf07bis = ["longueur ldlbullemere2 = ", len(ldlbullemere2)]
ldlinfos.append(inf07bis)

inf08 = ["longueur ldlbullemerefiltre = ", len(ldlbullemerefiltre1)]
ldlinfos.append(inf08)

inf09 = ["longueur ldlParametreBulle = ", len(ldlParametreBulle)]
ldlinfos.append(inf09)

inf10 = ["longueur ldlDataBubble = ", len(ldlDataBubble)]
ldlinfos.append(inf10)

inf11 = ["longueur ldlParametreBulle_2 = ", len(ldlParametreBulle2)]
ldlinfos.append(inf11)

inf12 = ["longueur DicoParametreBulle_2 = ", len(DicParametreBulle2)]
ldlinfos.append(inf12)

inf13 = ["longueur ldlMatBulleLaPlusLongue = ", len(ldlMatBulleLaPlusLongue)]
ldlinfos.append(inf13)

inf14 = ["longueur ldlMatRegExpBulleLaPlusLongue = ", len(ldlMatRegExpBulleLaPlusLongue)]
ldlinfos.append(inf14)

inf15 = ["longueur ldlRegExpBuLaplusLongue = ", len(ldlRegExpBuLaplusLongue)]
ldlinfos.append(inf15)

inf16 = ["longueur ldlFastaMapBuLaplusLongue = ", len(ldlFastaMatBuLaplusLongue)]
ldlinfos.append(inf16)

inf17 = ["longueur ldlFastaBuLaplusLongue = ", len(ldlFastaBuLaplusLongue)]
ldlinfos.append(inf17)

inf18 = ["longueur ldlTabVarPlsda1 = ", len(ldlTabVarPlsda1)]
ldlinfos.append(inf18)

inf19 = ["longueur ldlTabVarPlsMod1 = ", len(ldlTabVarPlsMod1)]
ldlinfos.append(inf19)

inf20 = ["longueur ldlComptecodons = ", len(ldlcodonsmere)]
ldlinfos.append(inf20)

inf21 = ["longueur ldlCodons = ", len(ldlCodons)]
ldlinfos.append(inf21)



print("----------------------")

fin1 = time.time()
print ("Time Prog before write %.4f secondes:" % (fin1 - debut1))

inf01 = ["Time prog before write in seconds", (fin1 - debut1)]
ldlinfos.append(inf01)


print("----------------Ecriture Fichier-------------------")

ecrit_ldl_tabOPLSDA(ldlTabVarOPLSDA, outfichier15)

fin1bis = time.time()
print("Time Prog after write %.4f secondes:" % (fin1bis - debut1))

inf01bis = ["Time after write in seconds", (fin1bis - debut1)]
ldlinfos.append(inf01bis)

#ecrit_ldl_3(ldlinfos, outfichier3)

print("----------------Fin -------------------")
print(ldlinfos)

