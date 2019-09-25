import sys, os, re, numpy as np, subprocess
from collections import defaultdict
from scipy.optimize import curve_fit
try:
    from .modules.configure import xrange, uopen, externals, logger
except :
    from modules.configure import xrange, uopen, externals, logger

def func_powerlaw(x, m, c ):
    return x**m * c

def splitGFF(gff, folder, prefix) :
    if not os.path.isdir(folder) :
        try :
            os.makedirs(folder)
        except FileExistsError :
            raise FileExistsError('Fail to create the output folder for GFF files')
    prev, fout = None, None
    with uopen(gff) as fin :
        for line in fin :
            header, line = line.split(':', 1)
            if header != prev :
                prev = header
                if fout :
                    fout.close()
                fout = uopen(os.path.join(folder, '{0}.{1}.gff.gz'.format(prefix.rsplit('/', 1)[-1], header)), 'w')
                fout.write('#!gff-version 3\n#!annotation-source PEPPA from enterobase.warwick.ac.uk\n')
            fout.write(line)
    fout.close()
    logger('GFF files are saved under folder {0}'.format(folder))

def getOrtho(fname) :
    ids = set([])
    groups = defaultdict(dict)
    with uopen(fname) as fin :
        for line in fin :
            part = line.strip().split('\t')
            tag = {'CDS':2, 'pseudogene':1}.get(part[1], 0)
            if tag > 0 :
                ID = re.findall(r'ID=([^;]+);', part[8])[0]
                if ID in ids : continue
                
                genome = part[0].split(':', 1)[0].replace('_genomic', '')

                ortho = re.findall(r'inference=ortholog_group:([^;]+)', part[8])[0]
                orthos = [ o.split(':')[1] for o in ortho.split(',') ]

                groups[genome].update({ortho : max(tag, groups[genome].get(ortho, 0)) for ortho in orthos})
                ids.add(ID)
    return groups

def writeCurve(prefix, groups, pseudogene=True, n_iter=300) :
    prefix2 = '{0}.{1}'.format(prefix, ['CDS_content', 'gene_content'][int(pseudogene)])
    gtype = ['CDSs', 'genes'][int(pseudogene)]
    encode = {}
    count = {}
    for grp in groups.values() :
        for g in grp.keys() :
            if g not in encode :
                encode[g] = len(encode)
    mat = [ set([encode[g] for g, c in grp.items() if pseudogene or c > 1]) for grp in groups.values() ]
    mat = [np.array(list(m)) for m in mat]
    ids, cnts = np.unique(np.concatenate(mat), return_counts=True)
    x = np.arange(len(mat))+1

    curves = np.zeros([len(mat), n_iter, 2], dtype=int)
    popts = np.zeros([n_iter, 2, 3])
    for ite in np.arange(n_iter) :
        for n in np.arange(5) :
            try :
                np.random.shuffle(mat)
                genes = np.zeros(len(encode), dtype=int)
                for id, m in enumerate(mat) :
                    genes[m] += 1
                    curves[id, ite, :] = (np.sum(genes>=1), np.sum(genes>=(id+1)) )
                popts[ite,   0,:2]=curve_fit(func_powerlaw, x, curves[:, ite, 0], maxfev=3000)[0]
                popts[ite,   0, 2]=(x[-1]+1)**popts[ite,  0,  0]*popts[ite,  0,  1] - (x[-1])**popts[ite,  0,  0]*popts[ite,  0,  1]
                popts[ite,   1,:2]=curve_fit(func_powerlaw, x[1:], curves[1:, ite, 0] - curves[:-1, ite, 0], maxfev=3000)[0]
                popts[ite,   1, 2]=(x[-1]+1)**popts[ite,  1,  0]*popts[ite,  1,  1]
                break
            except:
                continue
    popts[:, 1, 0] *= -1
    popt_sum = np.zeros([2, 3, 3])
    for i in np.arange(2) :
        for j in np.arange(3) :
            popt = np.sort(popts[:, i, j])
            popt_sum[i, j, :] = popt[int(popt.size*0.025)], popt[int(popt.size*0.5)], popt[int(popt.size*0.975)]
    with open('{0}.curve'.format(prefix2), 'w') as fout :
        fout.write('#! No. genomes: {0}\n'.format(len(groups)))
        fout.write('#! Ave. {1} per genome: {0:.03f}\n'.format(np.mean([m.size for m in mat]), gtype))
        fout.write('#! No. pan {1}: {0}\n'.format(len(encode), gtype))
        fout.write('#! No. core {1}: {0}\n'.format(curves[-1, 0, 1], gtype))
        fout.write('#! gamma (Heaps\' law model in DOI: 10.1016/j.mib.2008.09.006): {2:.03f}  CI95%: ({1:.03f} - {3:.03f}) {0}\n'.format(['Closed pan-genome! (gamma < 0)', 'Open pan-genome (gamma >= 0)'][popt_sum[0, 0, 1]>=0], *popt_sum[0, 0]))
        fout.write('#! alpha (Power\' law model in DOI: 10.1016/j.mib.2008.09.006): {2:.03f}  CI95%: ({1:.03f} - {3:.03f}) {0}\n'.format(['Closed pan-genome! (alpha > 1)', 'Open pan-genome (alpha <= 1)'][popt_sum[1, 0, 1]<=1], *popt_sum[1, 0]))
        fout.write('#No. genome\t(Pan-genome) Median\t2.5%\t97.5%\t|\t(Core-genome) Median\t2.5%\t97.5%\n')
        summary = np.zeros([len(mat), 2, 5], dtype=int)
        for id, curve in enumerate(curves) :
            pan  = np.sort(curve.T[0])
            core = np.sort(curve.T[1])
            pan_s =  [  pan[int(pan.size*0.025)],   pan[int(pan.size*0.5)],   pan[int(pan.size*0.975)]  ]
            core_s = [ core[int(core.size*0.025)], core[int(core.size*0.5)], core[int(core.size*0.975)] ]
            fout.write('{0}\t{1[1]}\t{1[0]}\t{1[2]}\t|\t{2[1]}\t{2[0]}\t{2[2]}\n'.format(id+1, pan_s, core_s))
    logger('Curves for {1} are saved in {0}.curve'.format(prefix2, ['CDS', 'all genes'][int(pseudogene)]))
    return popt_sum, summary

def writeMatrix(prefix, ortho) :
    genomes = sorted(ortho.keys())
    genes = defaultdict(int)
    for gs in ortho.values() :
        for g, n in gs.items() :
            genes[g] += int(n>0)
    genes = [g[0] for g in sorted(genes.items(), key=lambda x:(-x[1], x[0]))]
    with open('{0}.gene_content.matrix'.format(prefix), 'w') as fout :
        fout.write('\t'.join(['#']+genomes)+'\n')
        for g in genes :
            mat = [ ['-', 'p', 'c'][ortho[genome].get(g, 0)] for genome in genomes ]
            fout.write('\t'.join([g] + mat)+'\n')
    logger('Gene content matrix is saved in {0}.gene_content.matrix'.format(prefix))


def writeTree(prefix, ortho) :
    genes = sorted({g for gs in ortho.values() for g in gs})
    with open('{0}.gene_content.fas'.format(prefix), 'w') as fout :
        for genome, content in ortho.items() :
            fout.write('>{0}\n{1}\n'.format(genome, ''.join([['A','T'][int(g in content)] for g in genes ])))
    subprocess.Popen('{fasttree} -quiet -nt {0}.gene_content.fas > {0}.gene_content.nwk'.format(prefix, **externals), shell=True).wait()
    logger('Gene content tree is saved in {0}.gene_content.nwk'.format(prefix))

def PEPPA_parser(args) :
    param = arg_parser(args)
    if param.split :
        splitGFF(param.gff, param.split, param.prefix)
    if param.matrix or param.tree or param.curve :
        ortho = getOrtho(param.gff)
        if param.matrix :
            writeMatrix(param.prefix, ortho)
        if param.tree :
            writeTree(param.prefix, ortho)
        if param.curve :
            if param.curve % 2 > 0 :
                writeCurve(param.prefix, ortho, True)
            if param.curve >= 2 :
                writeCurve(param.prefix, ortho, False)
    
    
def arg_parser(a) :
    import argparse
    parser = argparse.ArgumentParser(description='''
PEPPA_parser.py 
(1) reads xxx.PEPPA.gff file
(2) split it into individual GFF files
(3) draw a present/absent matrix
(4) create a tree based on gene presence
(5) draw rarefraction curves of all genes and only intact CDSs
''', formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-g', '--gff', help='[REQUIRED] generated PEPPA.gff file from PEPPA.py.', required=True)
    parser.add_argument('-p', '--prefix', help='[Default: Same prefix as GFF input] Prefix for all outputs.', default=None)
    parser.add_argument('-s', '--split', help='[optional] A folder for splitted GFF files. ', default=None)
    parser.add_argument('-m', '--matrix', help='[Default: False] Flag to generate the gene present/absent matrix', default=False, action='store_true')
    parser.add_argument('-t', '--tree', help='[Default: False] Flag to generate the gene present/absent tree', default=False, action='store_true')
    parser.add_argument('-c', '--curve', help='[Default: 0] Choose to generate a rarefraction curve. \n0: No rarefraction curve. \n1: use all genes/pseudogenes. \n2: use only intact CDS. \n3: two curves for intact CDS only and for all genes/pseudogenes.', default=0, type=int)
    params = parser.parse_args(a)
    
    assert params.split or params.matrix or params.tree or params.curve, 'At least one type of output needs to be specified. '
    if not params.prefix : 
        params.prefix = params.gff.rsplit('.gff', 1)[0]
    
    return params
if __name__ == '__main__' :
    groups = PEPPA_parser(sys.argv[1:])
    
