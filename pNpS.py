from __future__ import print_function
import egglib
import sys
import gffutils
import os
import cyvcf
import numpy as np
import argparse
import gzip
# import random

STANDARD_CODON_TABLE = {
    "TTT": "F", "CTT": "L", "ATT": "I", "GTT": "V",
    "TTC": "F", "CTC": "L", "ATC": "I", "GTC": "V",
    "TTA": "L", "CTA": "L", "ATA": "I", "GTA": "V",
    "TTG": "L", "CTG": "L", "ATG": "M", "GTG": "V",

    "TCT": "S", "CCT": "P", "ACT": "T", "GCT": "A",
    "TCC": "S", "CCC": "P", "ACC": "T", "GCC": "A",
    "TCA": "S", "CCA": "P", "ACA": "T", "GCA": "A",
    "TCG": "S", "CCG": "P", "ACG": "T", "GCG": "A",

    "TAT": "Y", "CAT": "H", "AAT": "N", "GAT": "D",
    "TAC": "Y", "CAC": "H", "AAC": "N", "GAC": "D",
    "TAA": "*", "CAA": "Q", "AAA": "K", "GAA": "E",
    "TAG": "*", "CAG": "Q", "AAG": "K", "GAG": "E",

    "TGT": "C", "CGT": "R", "AGT": "S", "GGT": "G",
    "TGC": "C", "CGC": "R", "AGC": "S", "GGC": "G",
    "TGA": "*", "CGA": "R", "AGA": "R", "GGA": "G",
    "TGG": "W", "CGG": "R", "AGG": "R", "GGG": "G",
}


def get_degeneracy(codon_table):

    codon_degen = {}
    nucs = ['A', 'T', 'G', 'C']

    for codon in codon_table:
        degen_list = []
        for i in range(len(codon)):
            degen = 0
            for j in nucs:
                if i == 0:
                    mod_codon = j + codon[1:]
                elif i == 1:
                    mod_codon = codon[0] + j + codon[2]
                elif i == 2:
                    mod_codon = codon[0:2] + j

                if codon_table[codon] == codon_table[mod_codon]:
                    degen += 1

            degen_list.append(degen)
            for site in range(len(degen_list)):
                if degen_list[site] == 1:
                    degen_list[site] = 0

        codon_degen[codon] = degen_list

    return codon_degen


def get_seqid_str(ref_genome, chrom):

    """ Create the str to give to the limit parameter in the
    gffutils """

    seqid_str = None

    with egglib.io.fasta_iter(ref_genome) as f:
        for item in f:
            if item.name == chrom:
                length = item.ls
                seqid_str = chrom + ':1-%s' % length
                break

    return seqid_str


def find_overlapping_genes():
    """
    Use pybedtools intersect to get list of overlapping genes.

    """
    pass


def extract_cds_coords(db, seqid_str):
    """
    Extract the CDS coordinates for each protein in the GFF file

    :param db: gff database instance
    :param seqid_str: chr:start-end

    :return: dictionary with chromosome_genename as keys and each gene has a dictionary as a value
    containing the transcripts for that gene with the transcript id a key and the CDS
    coordinated as a list of tuple pairs containing the start and end coordinates for
    the CDS block (start, end)
    """

    cds_coords = {}
    for gene in db.features_of_type(featuretype='gene', limit=seqid_str):
        if gene['gene_biotype'][0] == 'protein_coding':
            gene_name = gene['ID'][0]
            gene_feat = db[gene_name]
            cds_coords[gene_name] = {}

            for i in db.children(gene_feat, featuretype='mRNA', order_by='start'):
                transcript_name = i['ID'][0]
                transcript_feat = db[transcript_name]
                cds_coords[gene_name][transcript_name] = []

                for j in db.children(transcript_feat, featuretype='CDS', order_by='start'):
                    cds_coords[gene_name][transcript_name].append((j.start, j.end))

                cds_coords[gene_name][transcript_name].append(gene.seqid)  # add chr info
                cds_coords[gene_name][transcript_name].append(gene.strand)  # add strand info
        else:
            continue

    #print(gene_count)
    return cds_coords


def extract_gene(vcf, gff, gene_name):

    """ Extract a gene from a VCF and return a cyvcf reader instance"""

    gene_feat = gff[gene_name]

    start = gene_feat.start
    end = gene_feat.end
    chrom = gene_feat.seqid

    try:
        return vcf.fetch(chrom, start, end)

    except ValueError:
        return None


def call_to_bases(vcf_site, site_pos, align, sample_list):

    """Set the ALign instance sites acording to the genotype information in
    the VCF file"""

    for sample in sample_list:
        index_1 = align.find(sample + '_1', index=True)
        index_2 = align.find(sample + '_2', index=True)
        call = vcf_site.genotype(sample)

        if 'PGT' in vcf_site.FORMAT:
            phased_gt = vcf_site.genotype(sample)['PGT']

        else:
            phased_gt = None

        if phased_gt is not None and vcf_site.ALT[0] is not None:

            allele_1 = phased_gt.split('|')[0]
            allele_2 = phased_gt.split('|')[1]

            if allele_1 == '0' and allele_2 == '1':
                base_1 = vcf_site.REF
                base_2 = vcf_site.ALT[0]
            elif allele_1 == '1' and allele_2 == '0':
                base_1 = vcf_site.ALT[0]
                base_2 = vcf_site.REF
            elif allele_1 == '1' and allele_2 == '1':
                base_1 = vcf_site.ALT[0]
                base_2 = vcf_site.ALT[0]

            # print(vcf_site)
            # print(vcf_site.ALT)
            # print(base_1, base_2)

            align.set(index_1, site_pos, base_1)
            align.set(index_2, site_pos, base_2)

        else:
            if len(vcf_site.REF) > 1:  # some monomorphic site with indel like ref allele in VCF

                align.set(index_1, site_pos, call.gt_bases.split('/')[0][0])
                align.set(index_2, site_pos, call.gt_bases.split('/')[1][0])

            else:
                align.set(index_1, site_pos, call.gt_bases[0])
                align.set(index_2, site_pos, call.gt_bases[2])


def rc_align(align):

    """ Reverse complement an egglib align"""

    for index in range(0, align.ns):

        seq = align.get_sequence(index).str()
        rc_seq = egglib.tools.rc(seq)

        align.set_sequence(index, rc_seq)

    return align


def find_prem_stop(cds_align):
    """Detect if there is a premature stop codon in a cds alignment"""

    for codon in egglib.tools.iter_stops(cds_align):  # check for premature stop codons
        if codon[1] < cds_align.ls - 3:

            return True
        else:
            return False


def extract_cds_align(vcf, min_dp, max_dp, sample_list, gff, gene_id,  cds_dict, filtered=True):

    """
    Iterator that goes site by site through vcf gene region and returns
    a CDS align and positions of sites in that align
    Returns:
        object: list
    """

    cds_dict_unique = {}  # Keep only only one of mulitple transcripts when they have exactly the same CDS coordinates
    transcript_coords = []
    for tran_coord in cds_dict[gene_id]:
        if cds_dict[gene_id][tran_coord] not in transcript_coords:
            transcript_coords.append(cds_dict[gene_id][tran_coord])
            cds_dict_unique[tran_coord] = cds_dict[gene_id][tran_coord]

    gene_cds_aligns = []

    for transcript in cds_dict_unique:
        cds_pos_list = []

        for coords in cds_dict[gene_id][transcript]:
            if type(coords) is tuple:
                cds_pos_list += range(coords[0], coords[1] + 1)

        len_bp = len(cds_pos_list)
        cds_align = egglib.Align(nsit=len_bp)

        for sample in sample_list:

            cds_align.add_sample(sample + '_1', data='N'*len_bp)
            cds_align.add_sample(sample + '_2', data='N'*len_bp)

        align_pos = -1
        indel_end = 0

        gene_region = extract_gene(vcf, gff, gene_id)

        if gene_region is None:  # deal with CDS on scaffolds which are returned as None
            continue

        for site in gene_region:

            if site.POS not in cds_pos_list:  # only extract genic sites in CDS blocks
                continue

            align_pos += 1

            if site.REF == 'N':
                # ref_N += 1
                continue

            if site.is_indel and site.aaf != 0.0:
                if len(site.REF) > len(site.ALT):
                    indel_end = site.POS + len(site.REF) - 1
                    continue
                # indels += 1
                else:
                    continue

            if len(site.ALT) >= 1 and site.ALT[-1] == '*':  # SNPs at spanning deletion
                # spanning_deletion += 1
                continue

            all_dp = [x['DP'] for x in site.samples]  # get the genotype depths

            if None in all_dp:  # Exclude sites where a genotype DP field is set to '.'
                # low_call_rate += 1
                continue

            mean_dp = np.mean(all_dp)

            if mean_dp < min_dp or mean_dp > max_dp:  # depth filter
                # extreme_depth += 1
                continue

            if 0 in all_dp or site.call_rate < 1.0:  # only consider sites where all samples have coverage
                # low_call_rate += 1
                continue

            if site.FILTER is not None:
                if site.FILTER == "REPEAT" or "REPEAT" in site.FILTER:  # exclude sites in repeat regions
                    # repeat_sites += 1
                    continue

            if site.is_monomorphic:
                if site.POS > indel_end:
                    call_to_bases(site, align_pos, cds_align, sample_list)
                    indel_end = 0
                    continue
                else:
                    continue

            if site.is_indel and site.aaf == 0.0:
                if site.POS > indel_end:
                    call_to_bases(site, align_pos, cds_align, sample_list)
                    indel_end = 0
                    # valid_sites += 1
                    continue
                else:
                    continue

            if site.is_snp:

                if len(site.ALT) > 1:
                    # multiallelic_snp += 1
                    continue

                if site.aaf == 1.0 or site.aaf == 0.0:  # only want SNPs polymorphic in our sample
                    if site.POS > indel_end:
                        call_to_bases(site, align_pos, cds_align, sample_list)
                        continue
                    else:
                        continue

                if filtered:
                    if site.FILTER == 'PASS' and site.POS > indel_end:
                        call_to_bases(site, align_pos, cds_align, sample_list)
                        continue
                    else:
                        # failed_snp += 1
                        continue
                else:
                    if site.POS > indel_end:
                        call_to_bases(site, align_pos, cds_align, sample_list)
                    else:
                        continue
            else:
                problem_site = site.CHROM + '\t' + str(site.POS)
                error_message = 'Could not assign ' + problem_site + ' to site type'
                sys.exit(error_message)

        if cds_dict[gene_id][transcript][-1] == '-':  # reverse-complement when on '-' strand
            cds_align = rc_align(cds_align)
            cds_pos_list.reverse()

        if find_prem_stop(cds_align) or cds_align.ls % 3 != 0:
            print('premature stop detected in or CDS not a multiple of 3 ', 'gene:', gene_id, 'transcript:', transcript)
            continue
        else:
            cds_align_nostop = cds_align.extract(0, cds_align.ls - 3)  # drop the stop codon at the end
            gene_cds_aligns.append((cds_align_nostop, cds_pos_list[0:-3]))
        # cds_align.to_fasta(gene_id + '_' + transcript + '.fas')

    if len(gene_cds_aligns) == 0:
        return None
    else:
        return gene_cds_aligns


def max_dict_val(d):
    """ return the key that has the highest value"""
    v = list(d.values())
    k = list(d.keys())
    return k[v.index(max(v))]


def degeneracy_over_all_transcripts(cds, codon_degen_dict):

    """ Identifies the 0-fold and 4-fold sites in a codon alignment
    and returns the list of site positions in the alignment when considering
    all transcripts for a gene (works also for genes with a single transcript).

    The first parameter aln_list is a list of Align instances, each item being the
    CDS region for each transript at a gene that has more.

     :return an Align instance containing 4-fold and an instance containing 0-fold site positions at a gene"""

    fourfold_pos_lists = [[] for tran in cds]  # list for each transcript
    zerofold_pos_lists = [[] for tran in cds]
    twofold_pos_lists = [[] for tran in cds]
    conflict_pos_lists = [[] for tran in cds]
    missing_pos_lists = [[] for tran in cds]

    tran_index = -1
    # tran_num = 0
    for transcript_cds_aln, genome_pos in cds:
        # print genome_pos
        # tran_num += 1
        # print transcript_cds_aln, genome_pos[0], genome_pos[-1]
        # transcript_cds_aln.to_fasta(fname=gene_id + '_' + str(tran_num) + '.fas')
        tran_index += 1
        start = 0
        for codon_aln in transcript_cds_aln.slider(3, 3):  # take codon by codon
            pos = genome_pos[start:start + 3]
            sites = {0: [], 1: [], 2: []}
            for sample in codon_aln:  # loop over samples
                codon_seq = sample.sequence.str()
                for i in range(3):
                    try:
                        sites[i].append(codon_degen_dict[codon_seq][i])
                    except KeyError:
                        sites[i].append('N')

            # print(pos, sites)

            for j in range(3):
                if len(set(sites[j])) != 1:  # site is not completely 4-fold, 2-fold or 0-fold degenerate across samples
                    if 'N' not in sites[j]:
                        conflict_pos_lists[tran_index].append(pos[j])
                    else:
                        missing_pos_lists[tran_index].append(pos[j])
                elif len(set(sites[j])) == 1:
                    if sites[j][0] == 0:
                        zerofold_pos_lists[tran_index].append(pos[j])
                    elif sites[j][0] == 2:
                        twofold_pos_lists[tran_index].append(pos[j])
                    elif sites[j][0] == 4:
                        fourfold_pos_lists[tran_index].append(pos[j])
                    elif sites[j][0] == 'N':
                        missing_pos_lists[tran_index].append(pos[j])  # missing data or conflicting degeneracy

            start += 3

        # print '0fd:', zerofold_pos_lists
        # print '2fd:', twofold_pos_lists
        # print '4fd:', fourfold_pos_lists
        # print 'Xfd:', conflict_pos_lists
        # print 'Nfd:', missing_pos_lists

    if len(cds) > 1:
        zerofold_sites = []
        fourfold_sites = []
        other_sites = []

        for zero_site in zerofold_pos_lists:
            zerofold_sites += zero_site

        for four_site in fourfold_pos_lists:
            fourfold_sites += four_site

        for two_site in twofold_pos_lists:
            other_sites += two_site

        for conflict_site in conflict_pos_lists:
            other_sites += conflict_site

        for missing_site in missing_pos_lists:
            other_sites += missing_site

        fourfold_pos_filtered = set(fourfold_sites).difference(set(zerofold_sites + other_sites))
        zerofold_pos_filtered = set(zerofold_sites).difference(set(fourfold_sites + other_sites))

        # print len(fourfold_pos_filtered), len(zerofold_pos_filtered)

        fourfold_align_list = []
        zerofold_align_list = []

        for cds_align, genome_pos in cds:
                fourfold_sites_to_extract = fourfold_pos_filtered.intersection(set(genome_pos))
                fourfold_sites_to_extract_index = [genome_pos.index(i) for i in list(fourfold_sites_to_extract)]
                fourfold_aln = cds_align.extract(fourfold_sites_to_extract_index)
                fourfold_align_list.append(fourfold_aln)
                fourfold_pos_filtered = fourfold_pos_filtered - fourfold_sites_to_extract

                zerofold_sites_to_extract = zerofold_pos_filtered.intersection(set(genome_pos))
                zerofold_sites_to_extract_index = [genome_pos.index(i) for i in list(zerofold_sites_to_extract)]
                zerofold_aln = cds_align.extract(zerofold_sites_to_extract_index)
                zerofold_align_list.append(zerofold_aln)
                zerofold_pos_filtered = zerofold_pos_filtered - zerofold_sites_to_extract

        four_fold_align = egglib.tools.concat(*fourfold_align_list)
        zero_fold_align = egglib.tools.concat(*zerofold_align_list)

        # print four_fold_align.ls, zero_fold_align.ls

    else:

        fourfold_pos_filtered = fourfold_pos_lists[0]

        fourfold_sites_to_extract_index = [cds[0][1].index(i) for i in fourfold_pos_filtered]
        four_fold_align = cds[0][0].extract(fourfold_sites_to_extract_index)

        zerofold_pos_filtered = zerofold_pos_lists[0]
        zerofold_sites_to_extract_index = [cds[0][1].index(i) for i in zerofold_pos_filtered]
        zero_fold_align = cds[0][0].extract(zerofold_sites_to_extract_index)

        # print len(fourfold_pos_filtered), len(zerofold_pos_filtered)
        # print four_fold_align.ls, zero_fold_align.ls

    return four_fold_align, zero_fold_align


def find_longest_transcript(cds_coord, db):

    """ Retrieve the longest transcript for each gene and return the CDS coordinates
    """

    cds_coord_longest_tran = {}

    for gene in cds_coord.keys():
        if len(cds_coord[gene].keys()) == 1:
            tran_id = cds_coord[gene].keys()[0]
            tran_dict = dict()
            tran_dict[tran_id] = cds_coord[gene][tran_id]
            cds_coord_longest_tran[gene] = tran_dict
            # print cds_coord_longest_tran
        else:
            tran_len_dict = {}
            for tran_id in cds_coord[gene].keys():
                transcript_feat = db[tran_id]
                transcript_len = 0
                for exon in db.children(transcript_feat, featuretype='exon', order_by='start'):
                    transcript_len += (exon.end - exon.start + 1)

                tran_len_dict[tran_id] = transcript_len
                longest_tran = max_dict_val(tran_len_dict)
                tran_dict = dict()
                tran_dict[longest_tran] = cds_coord[gene][longest_tran]
                cds_coord_longest_tran[gene] = tran_dict

    return cds_coord_longest_tran


# def choose_random_transcript(cds_coord, db):
#
#     cds_coord_random_tran = {}
#
#     for gene in cds_coord.keys():
#         if len(cds_coord[gene].keys()) == 1:
#             tran_id = cds_coord[gene].keys()[0]
#             cds_coord_random_tran[gene + '_' + tran_id] = cds_coord[gene][tran_id]
#             # print cds_coord_longest_tran
#         else:
#             random_tran = random.sample(cds_coord_random_tran[gene].keys(), 1)
#             print cds_coord_random_tran[gene].keys(), random_tran


def calc_polystats(align):

    cs = egglib.stats.ComputeStats()
    cs.add_stat('S')
    cs.add_stat('thetaW')
    cs.add_stat('Pi')
    cs.add_stat('D')
    cs.add_stat('ls')

    pol = cs.process_align(align)

    sample_size = align.ns

    delta_pi = calc_delta_pi(pol, sample_size)

    pol['delta_pi'] = delta_pi

    return pol


def calc_delta_pi(pol_stats, n):
    """
    Calculation of delta_pi from Langley et al. (2014) PLoS Genet 10(7): e1004457.

    :param pol_stats: Dictionary containing polymorphism stats needed to calculate delta_pi
    :param n: Sample size
    :return: float

    """
    S = pol_stats['S']
    k = pol_stats['Pi']

    if S > 0:
        assert isinstance(S, int)
        delta_pi = (k / float(S)) - (1 / sum(1.0 / i for i in range(1, n)))
    else:
        delta_pi = 'NA'

    return delta_pi


def main():

    parser = argparse.ArgumentParser(description="Program to calculate population genetic statistics from coding "
                                                 "regions in a VCF file")
    parser.add_argument('-i', '--in', dest='vcf', required=True, help="All-sites VCF file produced by GATK "
                                                                      "GenotypeGVCF")
    parser.add_argument('-g', '--gff', dest='gff', help="Reference genome annotation file in GFF3 format")
    parser.add_argument('-r', '--ref_genome', dest='ref_genome', help="Reference Genome in fasta format")
    parser.add_argument('-o', '--out', required=True, dest='outfile', help="Outfile to write results")
    parser.add_argument('-d', '--min_mean_dp', required=True, dest='min_dp', type=int,
                        help="Minimum mean depth across genotypes")
    parser.add_argument('-D', '--max_mean_dp', required=True, dest='max_dp', type=int,
                        help="Maximum mean depth across genotypes")
    parser.add_argument('-c', '--chromosome', required=True, dest='chrom', type=str,
                        help="Chromosome to analyse.")
    parser.add_argument('-m', '--method', required=True, dest='method', type=int,
                        help="Choice of transcript for degnerate sites:"
                             "( 1.) Longest transcript (-m 1)."
                             "( 2.) Choose method of that checks that sites are 0-fold or 4-fold degenerate across all"
                             "transcripts (-m 2).")

    args = parser.parse_args()

    vcf_file = cyvcf.Reader(filename=args.vcf)
    
    chrom_region = get_seqid_str(args.ref_genome, args.chrom)
    codon_degen_dict = get_degeneracy(STANDARD_CODON_TABLE)

    samples = vcf_file.samples

    min_depth = args.min_dp
    max_depth = args.max_dp

    gff_file = args.gff
    gff_str = ''
    with gzip.open(gff_file, 'r') as gff_infile:
        for line in gff_infile:
            if line.startswith(args.chrom):
                gff_str += line

    if os.path.isfile(gff_file[0:-6] + args.chrom + '.' + 'db'):
        gff_db = gffutils.FeatureDB(gff_file[0:-6] + args.chrom + '.' + 'db', keep_order=True)
    else:
        gffutils.create_db(gff_str, gff_file[0:-6] + args.chrom + '.' + 'db', merge_strategy="create_unique",
                           from_string=True)
        gff_db = gffutils.FeatureDB(gff_file[0:-6] + args.chrom + '.' + 'db', keep_order=True)

    cds_coords_dict = extract_cds_coords(gff_db, chrom_region)

    if args.method == 2:
        longest_trans_dict = find_longest_transcript(cds_coords_dict, gff_db)


    with open(args.outfile, 'w', 0) as outfile:
        print('gene', 'transcript', 'chr', 'fourfold_sites', 'fourfold_S', 'theta4', 'pi4', 'TajD4', 'delta_pi4',
              'zerofold_sites', 'zerofold_S', 'theta0', 'pi0', 'TajD0', 'delta_pi0', 'theta0_theta4', 'pi0_pi4',
              sep='\t', file=outfile)

        genes_processed = 0
        for gene_cds in cds_coords_dict:
            genes_processed += 1
            if args.method == 1:
                cds_alns = extract_cds_align(vcf_file, min_depth, max_depth, samples, gff_db, gene_cds, cds_coords_dict)
                transcript = 'NA'

            elif args.method == 2:
                cds_alns = extract_cds_align(vcf_file, min_depth, max_depth, samples, gff_db, gene_cds,
                                             longest_trans_dict)
                transcript = longest_trans_dict[gene_cds].keys()[0]

            else:
                sys.exit("Need to specify -m 1 or -m 2")

            if cds_alns is None:
                continue

            alns = degeneracy_over_all_transcripts(cds_alns, codon_degen_dict)

            fourfold_polystats = calc_polystats(alns[0])
            zerofold_polystats = calc_polystats(alns[1])

            if fourfold_polystats['ls'] is None or zerofold_polystats['ls'] is None:
                print('No called sites in', gene_cds)
                continue

            if fourfold_polystats['D'] is None:
                fourfold_polystats['D'] = 'NA'

            if zerofold_polystats['D'] is None:
                zerofold_polystats['D'] = 'NA'

            pi4 = fourfold_polystats['Pi'] / float(fourfold_polystats['ls'])
            pi0 = zerofold_polystats['Pi'] / float(zerofold_polystats['ls'])

            if pi4 == 0.0:
                pi0_pi4 = 'NA'
            else:
                pi0_pi4 = pi0 / pi4

            theta4 = fourfold_polystats['thetaW'] / float(fourfold_polystats['ls'])
            theta0 = zerofold_polystats['thetaW'] / float(zerofold_polystats['ls'])

            if theta4 == 0.0:
                theta0_theta4 = 'NA'
            else:
                theta0_theta4 = theta0 / theta4

            print(gene_cds, transcript, args.chrom, fourfold_polystats['ls'], fourfold_polystats['S'],
                  theta4, pi4, fourfold_polystats['D'], fourfold_polystats['delta_pi'],
                  zerofold_polystats['ls'], zerofold_polystats['S'],
                  theta0, pi0, zerofold_polystats['D'],  zerofold_polystats['delta_pi'],
                  theta0_theta4, pi0_pi4, sep='\t', file=outfile)

    print("Total Protein coding Genes processed", genes_processed)

if __name__ == '__main__':
    main()
