import sys
import re
import os
import mosaic
from Bio import AlignIO
import pandas as pd


def specfunc(name):
    specieslist_ca = [
        ['(CCDS|hg19|Homo_sapiens)', 'Hom'],
        ['(ENSPTR|Pan|Pan_troglodytes)', 'Pan'],
        ['(ENSGGO|Gor|Gorilla_gorilla)', 'Gor'],
        ['(ENSPPY|Pon|Pongo_abelii)', 'Pon'],
        ['(ENSMMU|Mac|rheMac|Macaca_mulatta)', 'Mac'],
        ['(ENSCJA|Cal|Callithrix_jacchus)', 'Cal'],
        ['(ENSOGA|Oto|Otolemur_garnettii)', 'Oto'],
        ['(ENSFCA|Fel|Felis_catus)', 'Fel'],
        ['(ENSBTA|Bos|BOVIN|Bos_taurus)', 'Bos'],
        ['(ENSECA|Equ|HORSE|Equus_caballus)', 'Equ']
    ]
    for xx in specieslist_ca:
        if re.search(xx[0], name, re.IGNORECASE):
            return xx[1]
    return 'Unknown'


def specfunc2(name):
    # Same as above, but it could be different!
    specieslist_ca = [
        ['(CCDS|hg19|Homo_sapiens)', 'Hom'],
        ['(ENSPTR|Pan|Pan_troglodytes)', 'Pan'],
        ['(ENSGGO|Gor|Gorilla_gorilla)', 'Gor'],
        ['(ENSPPY|Pon|Pongo_abelii)', 'Pon'],
        ['(ENSMMU|Mac|rheMac|Macaca_mulatta)', 'Mac'],
        ['(ENSCJA|Cal|Callithrix_jacchus)', 'Cal'],
        ['(ENSOGA|Oto|Otolemur_garnettii)', 'Oto'],
        ['(ENSFCA|Fel|Felis_catus)', 'Fel'],
        ['(ENSBTA|Bos|BOVIN|Bos_taurus)', 'Bos'],
        ['(ENSECA|Equ|HORSE|Equus_caballus)', 'Equ']
    ]
    for xx in specieslist_ca:
        if re.search(xx[0], name, re.IGNORECASE):
            return xx[1]
    return 'Unknown'


def getMSAlist_for_me(testfile, methodnames):
    MSAset1 = []
    MSAset2 = []
    methods = []
    for tt in ['AA', 'DNA']:
        for xx in methodnames:
            fname1 = 'aligned_%s_%s/%s' % (tt, xx, testfile)
            print "Looking for", fname1
            if os.path.exists(fname1):
                print "\tFound it!"
                with open(fname1) as fopen:
                    MSA = AlignIO.read(fopen, 'fasta')
                    if tt == 'AA':
                        MSAset1.append(MSA)
                        methods.append(xx)
                    else:
                        MSAset2.append(MSA)

    return MSAset1, MSAset2, methods


cutoffs_perID = {'Pan': .8247692861,
                 'Gor': .7730707713,
                 'Pon': .7627883286,
                 'Mac': .7527148696,
                 'Cal': .7396314922,
                 'Oto': .7039398156,
                 'Fel': .70,
                 'Bos': .70,
                 'Equ': .70
}

cutoffs_bit = {'Pan': 500,
               'Gor': 500,
               'Pon': 500,
               'Mac': 500,
               'Cal': 500,
               'Oto': 500,
               'Fel': 500,
               'Bos': 500,
               'Equ': 500
}

specorder = ['Hom', 'Pan', 'Gor',
             'Pon', 'Mac', 'Cal',
             'Oto', 'Fel', 'Bos',
             'Equ']

specnames = ['Hom', 'Pan', 'Gor',
             'Pon', 'Mac', 'Cal',
             'Oto', 'Fel', 'Bos',
             'Equ']

methodnames = ['inpara', 'blat', 'multiz', 'OMA']

if __name__ == '__main__':
    with open(sys.argv[1]) as f1:
        testfile_list = [xx.strip() for xx in f1.readlines() if xx.strip()]

    for optrule in ['pairwise', 'toref']:
        for edgefunc in ['perID', 'bitscore']:

            cutoffs = cutoffs_perID if edgefunc == 'perID' else cutoffs_bit

            for testfile in testfile_list:
                print '*********** Starting %s *******************' % testfile
                print "reading in first MSA set"
                MSAsetAA1, MSAsetDNA1, methods1 = \
                    getMSAlist_for_me(testfile, methodnames[:2])

                print "reading in second MSA set"
                MSAsetAA2, MSAsetDNA2, methods2 = \
                    getMSAlist_for_me(testfile, methodnames[2:])

                print "making first multiMSA"
                multiMSA1 = mosaic.MultiMSA(MSAsetAA1, MSAlist2=MSAsetDNA1,
                                            specfunc1=specfunc,
                                            methodnames1=methods1,
                                            methodnames2=methods1)

                print "making second multiMSA"
                multiMSA2 = mosaic.MultiMSA(MSAsetAA2, MSAlist2=MSAsetDNA2,
                                            specfunc1=specfunc2,
                                            methodnames1=methods2,
                                            methodnames2=methods2)

                print "Catenating multiMSAs"
                multiMSA = multiMSA1 + multiMSA2
                print multiMSA

                print "Running MOSAIC"
                optimizedcluster = mosaic.Mosaic(multiMSA, optrule=optrule,
                                                 ref='Hom',
                                                 ignoregaps=True,
                                                 # ignored if not perID&toref
                                                 edgefunc=edgefunc,
                                                 speccutoffs=cutoffs)

                dirname1 = '../unalignedAA_MOSAIC_%s_%s' % (optrule, edgefunc)
                dirname2 = dirname1.replace('AA', 'DNA')

                if not os.path.exists(dirname1):
                    os.mkdir(dirname1)
                if not os.path.exists(dirname2):
                    os.mkdir(dirname2)

                filename1 = '%s/%s' % (dirname1, testfile)
                filename2 = '%s/%s' % (dirname2, testfile)

                print "Writing unaligned files"
                optimizedcluster.write_unaligned(filename1=filename1,
                                                 filename2=filename2,
                                                 inclspec=True, inclmet=True,
                                                 specorder=['Hom', 'Pan',
                                                            'Gor',
                                                            'Pon', 'Mac',
                                                            'Cal',
                                                            'Oto', 'Fel',
                                                            'Bos',
                                                            'Equ']
                )

                dirname1_aligned = dirname1.replace('unaligned', 'aligned')
                dirname2_aligned = dirname2.replace('unaligned', 'aligned')

                if not os.path.exists(dirname1_aligned):
                    os.mkdir(dirname1_aligned)
                if not os.path.exists(dirname2_aligned):
                    os.mkdir(dirname2_aligned)

                filename1_aligned = filename1.replace('unaligned', 'aligned')
                filename2_aligned = filename2.replace('unaligned', 'aligned')

                print "\nAligning sequences"
                print "Output file 1 =", filename1_aligned
                print "Output file 2 =", filename2_aligned
                optimizedcluster.align(filename1_aligned,
                                       filename2=filename2_aligned)
