from __future__ import print_function
import egglib
import sys
import gffutils
import os
import cyvcf
import numpy as np
import argparse
import gzip


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

    seq = egglib.Align()

    with egglib.io.fasta_iter(ref_genome) as f:
        for item in f:
            if item.name == chrom:
                seq.add_sample(item.name, item.sequence.str())
                length = item.ls
                seqid_str = chrom + ':1-%s' % length
                break

    return seqid_str, seq


def extract_cds_coords(db, seqid_str, transcript_id):
    """
    Extract the CDS coordinates for each protein in the GFF file

    :param db: gff database instance
    :param seqid_str: chr:start-end

    :return: dictionary with gene name as keys and each gene has a dictionary as a value
    containing the transcripts for that gene with the transcript id a key and the CDS
    coordinates as a list of tuple pairs containing the start and end coordinates for
    the CDS block (start, end)
    """

    cds_coords = {}

    for gene in db.features_of_type(featuretype='gene', limit=seqid_str):
        try:
            gene_type = gene['gene_biotype'][0]
        except KeyError:
            gene_type = gene['biotype'][0] # ensembl files use biotype

        if gene_type == 'protein_coding':
            # try:
            name = gene['ID'][0]
            # except KeyError:
            #     name = gene['gene_id'][0]

            gene_feat = db[name]

            cds_coords[name] = {}

            for i in db.children(gene_feat, featuretype=transcript_id, order_by='start'):
                transcript_name = i['ID'][0]
                transcript_feat = db[transcript_name]
                cds_coords[name][transcript_name] = []

                for j in db.children(transcript_feat, featuretype='CDS', order_by='start'):
                    cds_coords[name][transcript_name].append((j.start, j.end))

                cds_coords[name][transcript_name].append(gene.seqid)  # add chr info
                cds_coords[name][transcript_name].append(gene.strand)  # add strand info
                try:
                    cds_coords[name][transcript_name].append(gene['Name'][0])  # append id for gene
                except KeyError:
                    cds_coords[name][transcript_name].append(gene['gene_id'][0])

        else:
            continue

    return cds_coords


def extract_ref_cds(chrom_seq, cds_coords):
    """
    Extract the CDS sequence for a gene transcript from the reference sequence

    :param chrom_seq:
    :param cds_coords:
    :return: CDS Align instance
    """

    cds_pos_list = []
    for coords in cds_coords:
        if type(coords) is tuple:
            cds_pos_list += range(coords[0], coords[1] + 1)

    cds_pos_list_zero = [site - 1 for site in cds_pos_list]

    cds_ref_seq = chrom_seq.extract(cds_pos_list_zero)

    if cds_coords[-2] == '-':
        cds_ref_seq = rc_align(cds_ref_seq)

    return cds_ref_seq


def extract_gene(vcf, gff, gene_name):

    """ Extract a gene from a VCF and return a cyvcf reader instance"""

    gene_feat = gff[gene_name]

    start = gene_feat.start
    end = gene_feat.end
    chrom = gene_feat.seqid

    # try:
    return vcf.fetch(chrom, start, end)

    # except ValueError:
    #     return None


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


def find_prem_stop(cds_ref_seq, gene_id, tran_id):

    """Detect if there is a premature stop codon in a cds alignment"""

    cds_ref_seq_nostop = cds_ref_seq.extract(0, cds_ref_seq.ls - 3)

    if egglib.tools.has_stop(cds_ref_seq_nostop):
        print('Premature stop codon detected  gene: %s transcript: %s' % (gene_id, tran_id))
        return True
    else:
        return False


def check_cds_incomplete(cds_ref_seq, gene_id, tran_id):

    stop_list = ['TGA', 'TAA', 'TAG']

    cds_seq = egglib.SequenceView(cds_ref_seq, 0, outgroup=False).str()

    if cds_seq[0:3] != 'ATG' or cds_seq[-3:] not in stop_list:
        print('No start, or no stop, codon detected gene: %s transcript: %s' % (gene_id, tran_id))
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

    gene_cds_aligns = []

    #for transcript in cds_dict[gene_id]:
    transcript = cds_dict[gene_id].keys()[0]
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

    if cds_dict[gene_id][transcript][-2] == '-':  # reverse-complement when on '-' strand
        cds_align = rc_align(cds_align)
        cds_pos_list.reverse()

    cds_align_nostop = cds_align.extract(0, cds_align.ls - 3)  # drop the stop codon at the end

    gene_cds_aligns.append((cds_align_nostop, cds_pos_list[0:-3]))

    return cds_align_nostop, cds_pos_list[0:-3]


def max_dict_val(d):
    """ return the key that has the highest value"""
    v = list(d.values())
    k = list(d.keys())
    return k[v.index(max(v))]


def extract_degenerate_sites(cds, codon_degen_dict):

    """ Identifies the 0-fold and 4-fold sites in a codon alignment
    and returns the list of site positions in the alignment when considering
    all transcripts for a gene (works also for genes with a single transcript).

    The first parameter aln_list is a list of Align instances, each item being the
    CDS region for each transript at a gene that has more.

     :return an Align instance containing 4-fold and an instance containing 0-fold site positions at a gene"""

    fourfold_pos_list = []
    zerofold_pos_list = []

    transcript_cds_aln = cds[0]
    genome_pos = cds[1]
    start = 0
    for codon_aln in transcript_cds_aln.slider(3, 3):  # take codon by codon
        codon_pol = calc_polystats(codon_aln)
        if codon_pol['S'] > 1:  # exclude codons with more than one snp
            start += 3
            continue

        if codon_pol['ls'] < 3:  # exclude any codon with missing data
            start += 3
            continue

        pos = genome_pos[start:start + 3]
        sites = {0: [], 1: [], 2: []}
        for sample in codon_aln:  # loop over samples
            codon_seq = sample.sequence.str()
            for i in range(3):
                # try:
                sites[i].append(codon_degen_dict[codon_seq][i])
                # except KeyError:
                #     sites[i].append('N')

        #print(pos, sites)

        for j in range(3):
            if len(set(sites[j])) == 1:  # Check that site is completely 4-fold, 2-fold or 0-fold degenerate across samples
                if sites[j][0] == 0:
                    zerofold_pos_list.append(pos[j])
                elif sites[j][0] == 4:
                    fourfold_pos_list.append(pos[j])
                else:
                    continue
            else:
                continue

        start += 3

    fourfold_sites_to_extract_index = [cds[1].index(i) for i in fourfold_pos_list]
    four_fold_align = cds[0].extract(fourfold_sites_to_extract_index)

    zerofold_sites_to_extract_index = [cds[1].index(i) for i in zerofold_pos_list]
    zero_fold_align = cds[0].extract(zerofold_sites_to_extract_index)

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


def choose_transcripts_from_list(cds_coord, ortholog_dict):
    cds_coord_from_list = {}

    for gene in cds_coord.keys():
        tran_id = cds_coord[gene].keys()[0]
        gene_id = cds_coord[gene][tran_id][-1]
        if gene_id in ortholog_dict.keys():
            ortho_tran = ortholog_dict[gene_id]
            tran_dict = dict()
            tran_dict[ortho_tran] = cds_coord[gene][tran_id]
            cds_coord_from_list[gene] = tran_dict

    return cds_coord_from_list


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


def calc_sfs(aln):

    """ Calculates the folded sfs for an alignment

    The sfs is returned as a list and also contains the counts of monomorphic
    sites as the first element"""

    sfs_folded = [0 for i in range(int(aln.ns / 2) + 1)]
    for site_index in range(aln.ls):
        site_freq = egglib.stats.SiteFrequency.make_from_dataset(aln, site_index)
        allele_freq = site_freq.freqs(no_pop=True)
        if len(allele_freq) == 1:
            sfs_folded[0] += 1
        else:
            minor_allele_freq = min(allele_freq)
            sfs_folded[minor_allele_freq] += 1

    return sfs_folded


def extract_ww_ss(aln):
    site_indices = []
    for site_index in range(aln.ls):
        col = aln.column(site_index,ingroup=True, outgroup=False)
        if len(set(map(chr, col))) == 1:
            site_indices.append(site_index)
        else:
            alleles = list(set(map(chr, col)))
            if alleles[0] == 'A' or alleles[0] == 'T':
                if alleles[1] == 'A' or alleles[1] == 'T':
                    site_indices.append(site_index)
            elif alleles[0] == 'G' or alleles[0] == 'C':
                if alleles[1] == 'G' or alleles[1] == 'C':
                    site_indices.append(site_index)

    ww_ss_aln = aln.extract(site_indices)

    return ww_ss_aln



def main():

    parser = argparse.ArgumentParser(description="Program to calculate population genetic statistics from coding "
                                                 "regions in a VCF file")
    parser.add_argument('-i', '--in', dest='vcf', required=True, help="All-sites VCF file produced by GATK "
                                                                      "GenotypeGVCF")
    parser.add_argument('-g', '--gff', dest='gff', help="Reference genome annotation file in GFF3 format")
    parser.add_argument('-t', '--type_gff', dest='type', help="Specify whether GFF is Ensembl or NCBI")
    parser.add_argument('-r', '--ref_genome', dest='ref_genome', help="Reference Genome in fasta format")
    parser.add_argument('-o', '--out', required=True, dest='outfile', help="Outfile to write results")
    parser.add_argument('-d', '--min_mean_dp', required=True, dest='min_dp', type=int,
                        help="Minimum mean depth across genotypes")
    parser.add_argument('-D', '--max_mean_dp', required=True, dest='max_dp', type=int,
                        help="Maximum mean depth across genotypes")
    parser.add_argument('-c', '--chromosome', required=True, dest='chrom', type=str,
                        help="Chromosome to analyse.")
    parser.add_argument('-m', '--method', required=True, dest='method', type=int,
                        help="Choice of transcript for degenerate sites:"
                             "( 1.) Longest transcript (-m 1)."
                             "( 2.) Choose transcript from the an ortholog list file specified")
    parser.add_argument('-f', '--ortholog_file', dest='orthologs', type=str, help="Tabe delimited list of orthologs "
                                                                                  "between species")
    parser.add_argument('-s', '--species', dest='species', type=int, help="Specify the column number for the"
                                                                          "species identifier in the ortholog list file"
                                                                          "Required when -m 2 option is used")

    args = parser.parse_args()

    if args.method == 3 and args.orthologs is None:
        sys.exit("\nNeed to specify an ortholog list file with -f option\n")

    if args.method == 3 and args.species is None:
        sys.exit("\nSpecify the column number for the species identifier in the ortholog list file\n")

    vcf_file = cyvcf.Reader(filename=args.vcf)

    chrom = get_seqid_str(args.ref_genome, args.chrom)

    chrom_region = chrom[0]
    chrom_seq = chrom[1]

    codon_degen_dict = get_degeneracy(STANDARD_CODON_TABLE)

    samples = vcf_file.samples

    min_depth = args.min_dp
    max_depth = args.max_dp

    gff_file = args.gff
    gff_str = ''
    with gzip.open(gff_file, 'r') as gff_infile:
        for line in gff_infile:
            gff_col = line.split('\t')
            if gff_col[0] == args.chrom:
                gff_str += line

    if os.path.isfile(gff_file[0:-6] + args.chrom + '.' + 'db'):
        gff_db = gffutils.FeatureDB(gff_file[0:-6] + args.chrom + '.' + 'db', keep_order=True)
    else:
        gffutils.create_db(gff_str, gff_file[0:-6] + args.chrom + '.' + 'db', merge_strategy="create_unique",
                           from_string=True)
        gff_db = gffutils.FeatureDB(gff_file[0:-6] + args.chrom + '.' + 'db', keep_order=True)

    if args.type == 'Ensembl':
        transcript_str = 'transcript'
    elif args.type == 'NCBI':
        transcript_str = 'mRNA'
    else:
        sys.exit('Need to specify -t Ensembl or -t NCBI')

    cds_coords_dict = extract_cds_coords(gff_db, chrom_region, transcript_str)

    if args.method == 1:
        longest_trans_dict = find_longest_transcript(cds_coords_dict, gff_db)

    elif args.method == 2 and args.orthologs:
        species_col_num = args.species - 1
        with open(args.orthologs) as infile:
            ortholog_list_dict = {}
            ortholog_num_dict = {}
            for line in infile:
                col = line.rstrip().split(' ')
                ortholog = col[species_col_num:species_col_num + 3]
                ortholog_list_dict[ortholog[2].split('=')[1]] = ortholog[1]

                ortholog_num_dict[ortholog[2].split('=')[1]] = col[-1]

        ortho_trans_dict = choose_transcripts_from_list(cds_coords_dict, ortholog_list_dict)

    with open(args.outfile, 'w', 0) as outfile:
        if args.method == 1:
            print('gene', 'gene_id', 'transcript', 'chr', 'start', 'end', 'fourfold_sites', 'fourfold_S', 'theta4',
                  'pi4', 'TajD4', 'delta_pi4', 'zerofold_sites', 'zerofold_S', 'theta0', 'pi0', 'TajD0', 'delta_pi0',
                  'theta0_theta4', 'pi0_pi4', 'fourfold_sfs', 'zerofold_sfs', 'fourfold_sites_ww_ss', 'fourfold_S_ww_ss'
                  , 'theta4_ww_ss', 'pi4_ww_ss', 'TajD4_ww_ss', 'delta_pi4_ww_ss', 'zerofold_sites_ww_ss',
                  'zerofold_S_ww_ss', 'theta0_ww_ss', 'pi0_ww_ss', 'TajD0_ww_ss', 'delta_pi0_ww_ss',
                  'theta0_theta4_ww_ss', 'pi0_pi4_ww_ss', 'fourfold_sfs_ww_ss', 'zerofold_sfs_ww_ss', sep='\t',
                  file=outfile)
        elif args.method == 2:
            print('gene', 'gene_id', 'transcript', 'chr', 'start', 'end', 'fourfold_sites', 'fourfold_S', 'theta4',
                  'pi4', 'TajD4', 'delta_pi4', 'zerofold_sites', 'zerofold_S', 'theta0', 'pi0', 'TajD0', 'delta_pi0',
                  'theta0_theta4', 'pi0_pi4', 'fourfold_sfs', 'zerofold_sfs', 'fourfold_sites_ww_ss', 'fourfold_S_ww_ss'
                  , 'theta4_ww_ss', 'pi4_ww_ss', 'TajD4_ww_ss', 'delta_pi4_ww_ss', 'zerofold_sites_ww_ss',
                  'zerofold_S_ww_ss', 'theta0_ww_ss', 'pi0_ww_ss', 'TajD0_ww_ss', 'delta_pi0_ww_ss',
                  'theta0_theta4_ww_ss', 'pi0_pi4_ww_ss', 'fourfold_sfs_ww_ss', 'zerofold_sfs_ww_ss', 'ortholog_num',
                  sep='\t', file=outfile)

        genes_processed = 0
        for gene_cds in cds_coords_dict:
            genes_processed += 1

            if args.method == 1:

                transcript = longest_trans_dict[gene_cds].keys()[0]

                cds_ref_seq = extract_ref_cds(chrom_seq, longest_trans_dict[
                    gene_cds][transcript])

                if check_cds_incomplete(cds_ref_seq, gene_cds,
                                        transcript):  # check CDS has valid start codon and stop codon
                    # cds_ref_seq.to_fasta('%s_%s.fas' % (gene_cds, transcript))
                    continue

                if find_prem_stop(cds_ref_seq, gene_cds, transcript):  # check for premature stop
                    # cds_ref_seq.to_fasta('%s_%s.fas' % (gene_cds, transcript))
                    continue

                cds_alns = extract_cds_align(vcf_file, min_depth, max_depth, samples, gff_db, gene_cds,
                                             longest_trans_dict)

                # cds_alns[0].to_fasta('%s_%s.fas' % (gene_cds, transcript))

            elif args.method == 2:

                if gene_cds not in ortho_trans_dict.keys():
                    # print("%s not in otholog list" % gene_cds)
                    continue

                transcript = ortho_trans_dict[gene_cds].keys()[0]
                cds_ref_seq = extract_ref_cds(chrom_seq, transcript)

                if check_cds_incomplete(cds_ref_seq, gene_cds, transcript):
                    # cds_ref_seq.to_fasta('%s_%s.fas' % (gene_cds, transcript))
                    continue

                if find_prem_stop(cds_ref_seq, gene_cds, transcript):
                    # cds_ref_seq.to_fasta('%s_%s.fas' % (gene_cds, transcript))
                    continue

                cds_alns = extract_cds_align(vcf_file, min_depth, max_depth, samples, gff_db, gene_cds,
                                             ortho_trans_dict)

                gene_id = ortho_trans_dict[gene_cds][transcript][-1]
                ortholog_num = ortholog_num_dict[gene_id]

            else:
                sys.exit("Need to specify -m 1 or -m 2")

            alns = extract_degenerate_sites(cds_alns, codon_degen_dict)

            fourfold_polystats = calc_polystats(alns[0])
            zerofold_polystats = calc_polystats(alns[1])

            if fourfold_polystats['ls'] is None or zerofold_polystats['ls'] is None:
                print('No called sites in', gene_cds)
                continue

            fourfold_sfs_str = ','.join(str(s) for s in calc_sfs(alns[0]))
            zerofold_sfs_str = ','.join(str(s) for s in calc_sfs(alns[1]))

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

            # WW and SS stats

            fourfold_ww_ss_aln = extract_ww_ss(alns[0])
            zerofold_ww_ss_aln = extract_ww_ss(alns[1])

            fourfold_polystats_ww_ss = calc_polystats(fourfold_ww_ss_aln)
            zerofold_polystats_ww_ss = calc_polystats(zerofold_ww_ss_aln)

            fourfold_sfs_ww_ss_str = ','.join(str(s) for s in calc_sfs(fourfold_ww_ss_aln))
            zerofold_sfs_ww_ss_str = ','.join(str(s) for s in calc_sfs(zerofold_ww_ss_aln))

            if fourfold_polystats_ww_ss['D'] is None:
                fourfold_polystats_ww_ss['D'] = 'NA'

            if zerofold_polystats_ww_ss['D'] is None:
                zerofold_polystats_ww_ss['D'] = 'NA'

            pi4_ww_ss = fourfold_polystats_ww_ss['Pi'] / float(fourfold_polystats_ww_ss['ls'])
            pi0_ww_ss = zerofold_polystats_ww_ss['Pi'] / float(zerofold_polystats_ww_ss['ls'])

            if pi4_ww_ss == 0.0:
                pi0_pi4_ww_ss = 'NA'
            else:
                pi0_pi4_ww_ss = pi0_ww_ss / pi4_ww_ss

            theta4_ww_ss = fourfold_polystats_ww_ss['thetaW'] / float(fourfold_polystats_ww_ss['ls'])
            theta0_ww_ss = zerofold_polystats_ww_ss['thetaW'] / float(zerofold_polystats_ww_ss['ls'])

            if theta4_ww_ss == 0.0:
                theta0_theta4_ww_ss = 'NA'
            else:
                theta0_theta4_ww_ss = theta0_ww_ss / theta4_ww_ss

            gene_feat = gff_db[gene_cds]
            chrom = gene_feat.seqid
            start = gene_feat.start
            end = gene_feat.end

            if args.method == 1:
                print(gene_cds, cds_coords_dict[gene_cds][transcript][-1], transcript, chrom, start, end,
                      fourfold_polystats['ls'], fourfold_polystats['S'], theta4, pi4, fourfold_polystats['D'],
                      fourfold_polystats['delta_pi'], zerofold_polystats['ls'], zerofold_polystats['S'], theta0,
                      pi0, zerofold_polystats['D'],  zerofold_polystats['delta_pi'], theta0_theta4, pi0_pi4,
                      fourfold_sfs_str, zerofold_sfs_str, fourfold_polystats_ww_ss['ls'], fourfold_polystats_ww_ss['S'],
                      theta4_ww_ss, pi4_ww_ss, fourfold_polystats_ww_ss['D'], fourfold_polystats_ww_ss['delta_pi'],
                      zerofold_polystats_ww_ss['ls'], zerofold_polystats_ww_ss['S'], theta0_ww_ss, pi0_ww_ss,
                      zerofold_polystats_ww_ss['D'], zerofold_polystats_ww_ss['delta_pi'], theta0_theta4_ww_ss,
                      pi0_pi4_ww_ss, fourfold_sfs_ww_ss_str, zerofold_sfs_ww_ss_str, sep='\t', file=outfile)
            elif args.method == 2:
                print(gene_cds, cds_coords_dict[gene_cds][transcript][-1], transcript, chrom, start, end,
                      fourfold_polystats['ls'], fourfold_polystats['S'], theta4, pi4, fourfold_polystats['D'],
                      fourfold_polystats['delta_pi'], zerofold_polystats['ls'], zerofold_polystats['S'], theta0,
                      pi0, zerofold_polystats['D'],  zerofold_polystats['delta_pi'], theta0_theta4, pi0_pi4,
                      fourfold_sfs_str, zerofold_sfs_str, fourfold_polystats_ww_ss['ls'], fourfold_polystats_ww_ss['S'],
                      theta4_ww_ss, pi4_ww_ss, fourfold_polystats_ww_ss['D'], fourfold_polystats_ww_ss['delta_pi'],
                      zerofold_polystats_ww_ss['ls'], zerofold_polystats_ww_ss['S'], theta0_ww_ss, pi0_ww_ss,
                      zerofold_polystats_ww_ss['D'], zerofold_polystats_ww_ss['delta_pi'], theta0_theta4_ww_ss,
                      pi0_pi4_ww_ss, fourfold_sfs_ww_ss_str, zerofold_sfs_ww_ss_str, ortholog_num, sep='\t',
                      file=outfile)

    print("Total Protein coding Genes processed", genes_processed)

if __name__ == '__main__':
    main()
