import os
import subprocess
import pandas as pd
import numpy as np


def getOverlap(a,b):
    return max(0,min(a[1],b[1]) - max(a[0], b[0])+1)



def map_function_over_all_files(data_path,
                                file_end,
                                func_to_loop, **kwargs):
    '''
    Loop a function over all files with a certain file-end in the data_path.
    For example:

        for fqgz_file in os.listdir(fastq_path):
            if fqgz_file.endswith(".fq.gz"):
                unpack_fq_file(fastq_path, fqgz_file)

    This can be generalized by map():

        files = [fastq_path + fn for fn in os.listdir(fastq_path) if fn.endswith('fq.gz)]
        list(map(unpack_fq_file, files))
    '''

    files = [data_path + fn for fn in os.listdir(data_path) if fn.endswith(file_end)]

    if kwargs:
        [func_to_loop(a, **kwargs) for a in files]
    else:
        list(map(func_to_loop, files))

    return


def prefilter_vcf(vcf_file):
    """
    Only mutations that pass the call-specific filter are considered further on.
    :param vcf_file: Output VCF file of Manta
    """

    try:
        subprocess.run("grep 'PASS\t' "+ vcf_file +" > " + vcf_file + '_PASS', shell=True, check=True)
    except:
        print(vcf_file)

    return

def add_vcf_info(vcf_pass_file, sample_dict):
    '''
    :param vcf_pass_file: One VCF file after filtering with prefilter_vcf()
    :param sample_dict: Dictionary containing meta information for the samples
    :return: 
    '''
    # get the meta info into a proper format
    vcf = pd.read_csv(vcf_pass_file, sep='\t', header=None)
    vcf.columns = ['chr', 'pos', 'id', 'ref', 'alt', 'quality', 'filter', 'info', 'format', 'sample1']
    vcf['Type'] = ''
    vcf.loc[vcf.id.str.contains('DEL'), 'Type'] = 'Deletion'
    vcf.loc[vcf.id.str.contains('INS'), 'Type'] = 'Insertion'
    vcf.loc[vcf.id.str.contains('DUP'), 'Type'] = 'Tandem Duplication'
    vcf.loc[vcf.id.str.contains('BND'), 'Type'] = 'Translocation'

    vcf['Zygosity'] = ''
    vcf.loc[(vcf.sample1.str.contains('1/1')), 'Zygosity'] = 'Homozygous'
    vcf.loc[~(vcf.sample1.str.contains('1/1')), 'Zygosity'] = 'Heterozygous'
    vcf['File'] = vcf_pass_file.split('/')[-1].split('_')[0]
    vcf['Treatment'] = sample_dict[int(vcf_pass_file.split('/')[-1].split('.')[0])]

    return vcf


def start_analysis(vcf_path):
    '''
    Preprocess VCF files, join to one DataFrame, add meta information, and save as CSV
    :param vcf_path: Path with all VCF files
    '''

    # For each VCF file in the path get only the PASS-SV
    map_function_over_all_files(data_path=vcf_path, file_end='vcf', func_to_loop=prefilter_vcf)

    # get the sample info and save into the sample_dict
    met = pd.read_csv(f'{vcf_path}/sample_info.csv', sep=',')
    met['Name'] = met.Generation + '_' + met.Strain + '_' + met.Treatment + '_' + met.Sex
    sample_dict = {}
    for i in range(len(met)):
        sample_dict[met.iloc[i].CCG] = met.iloc[i].Name

    # load the files,  append to one df, and add meta information
    all_vcf = pd.DataFrame()
    for v in os.listdir(vcf_path):
        if v.endswith('PASS'):
            try:
                tmp = add_vcf_info(f'{vcf_path}{v}', sample_dict)
                all_vcf = all_vcf.append(tmp)
            except:
                print(v)

    # get the end position and chromosome
    all_vcf['end_pos'] = np.nan
    all_vcf['end_chr'] = np.nan
    all_vcf.loc[all_vcf.Type == 'Translocation', 'end_chr'] = all_vcf.loc[all_vcf.Type == 'Translocation'].alt.str.split(r'\[|]|:', expand=True)[1]
    all_vcf.loc[all_vcf.Type == 'Translocation', 'end_pos'] = all_vcf.loc[all_vcf.Type == 'Translocation'].alt.str.split(r'\[|]|:', expand=True)[2]

    all_vcf.loc[all_vcf.Type == 'Deletion', 'end_pos'] = all_vcf.loc[all_vcf.Type == 'Deletion']['info'].str.split(r'\=|;', expand=True)[1]
    all_vcf.loc[all_vcf.Type == 'Deletion', 'end_chr'] = all_vcf.loc[all_vcf.Type == 'Deletion']['chr']
    all_vcf.loc[all_vcf.Type == 'Insertion', 'end_pos'] = all_vcf.loc[all_vcf.Type == 'Insertion']['info'].str.split(r'\=|;', expand=True)[1]
    all_vcf.loc[all_vcf.Type == 'Insertion', 'end_chr'] = all_vcf.loc[all_vcf.Type == 'Insertion']['chr']
    all_vcf.loc[all_vcf.Type == 'Tandem Duplication', 'end_pos'] = all_vcf.loc[all_vcf.Type == 'Tandem Duplication']['info'].str.split(r'\=|;', expand=True)[1]
    all_vcf.loc[all_vcf.Type == 'Tandem Duplication', 'end_chr'] = all_vcf.loc[all_vcf.Type == 'Tandem Duplication']['chr']
    all_vcf.end_pos = all_vcf.end_pos.astype(int)
    all_vcf.to_csv(f'{vcf_path}/all_vcf.csv', index=False)
    return all_vcf




def make_int(all_vcf, col):
	all_vcf.loc[all_vcf[col].isna(), col] = 0
	all_vcf[col] = all_vcf[col].astype(int)
	return all_vcf

def get_translocation_type(tmp):
    '''

    :param tmp: one ALT field (tmp = filtered_vcf.alt.iloc[i])
    :return: type (1,2,3,4 or -1)
    '''
    tmptype = -1 # no translocation
    insert = np.nan # some translocations have an insertion between their break points
    if '[' in tmp:
        tmpa = tmp.split('[')
        if tmpa[0] == '':
            tmptype = 4  # [p[t
            if len(tmpa[2]) > 1:  # insert
                insert = tmpa[2][:-1] # last is reference
        else:
            tmptype = 1  # t[p[
            if len(tmpa[0]) > 1:  # insert
                insert = tmpa[0][1:] # first is reference
    elif ']' in tmp:
        tmpa = tmp.split(']')
        if tmpa[0] == '':
            tmptype = 3  # ]p]t
            if len(tmpa[2]) > 1:  # insert
                insert = tmpa[2][:-1] # last is reference
        else:
            tmptype = 2  # t]p]
            if len(tmpa[0]) > 1:  # insert
                insert = tmpa[0][1:] # first is reference
    return tmptype, insert



def mark_indel_background(P0, F1, all_vcf):
    '''
    Loop through all F1 INDELs and mark those that overlap with INDELs in any P0 sample
    :param P0: Subset of all_vcf with only P0 INDELs
    :param F1: Subset of all_vcf with only F1 INDELs
    :param all_vcf: The whole dataset
    :return: 
    '''
    for _, F1row in F1.iterrows():
        # Max_CI contains the maximum of the confidence intervals around the break points
        # The position of some structural variants are hard to pinpoint to a position due to sequence homology
        CI_area_F1 = F1row.Max_CI 
        current_pos = [F1row.pos - CI_area_F1, F1row.end_pos + CI_area_F1]
        current_chr = F1row.chr
        current_ind = F1row.name

        P0_tmp = P0[(P0.end_chr == current_chr)]
        for _, P0row in P0_tmp.iterrows():
            CI_area_P0 = P0row.Max_CI
            P0_pos = [P0row.pos - CI_area_P0, P0row.end_pos + CI_area_P0]
            # mark the INDEL if it overlaps with any other INDEL
            if ((getOverlap(current_pos, P0_pos) > 0)):
                all_vcf.loc[current_ind, 'Marked'] = 'x'
                break
            if P0_pos[0] > current_pos[1]:  # if the P0_pos is higher than F1_endpos they can not match anymore
                break

    return all_vcf


def mark_translocation_background(P0, F1, all_vcf):
    '''
    Loop through all F1 Translocations are mark those that overlap with Translocations in any P0 sample
    :param P0: Subset of all_vcf with only P0 Translocations
    :param F1: Subset of all_vcf with only F1 Translocations
    :param all_vcf: The whole dataset
    :return: 
    '''
    for _, F1row in F1.iterrows():

        # get the highest absolute number of the confidence intervals and use it as the uncertainty region for both break points
        transloc_area_F1 = F1row.Max_CI
        current_start_pos = [F1row.pos -transloc_area_F1, F1row.pos +transloc_area_F1]
        current_end_pos = [F1row.end_pos -transloc_area_F1, F1row.end_pos +transloc_area_F1]
        current_start_chr = F1row.chr
        current_end_chr = F1row.end_chr
        current_ind = F1row.name

        P0_tmp = P0[(P0.chr == current_start_chr)&(P0.end_chr == current_end_chr)]
        for _,P0row in P0_tmp.iterrows():  
            transloc_area_P0 = P0row.Max_CI 
            P0_start_pos = [P0row.pos -transloc_area_P0, P0row.pos +transloc_area_P0]
            P0_end_pos = [P0row.end_pos -transloc_area_P0, P0row.end_pos +transloc_area_P0]
            # only mark a translocation if both the start and end break point overlap
            if ((getOverlap(current_start_pos ,P0_start_pos ) >0) & (getOverlap(current_end_pos ,P0_end_pos) > 0)):
                all_vcf.loc[current_ind, 'Marked'] = 'x'
                break
            if P0_start_pos[0] > current_start_pos[1]: # if the P0_pos is higher than F1_endpos they can not match anymore
                break

    return all_vcf




def loop_repeat_overlap(all_vcf, reps_tmp_filtered, current_pos, current_ind, repeat_col):
    '''
    Loop over all potential repeats to annotate whether and if so what repeat is overlapping with 
    the current structural variant
    
    :param all_vcf: The whole dataset
    :param reps_tmp_filtered: Pandas dataframe with the potentially relevant repeat regions
    :param current_pos: List with 2 elements [start_pos, end_pos], either the start and end of an INDEL
                        or one flank with a CI of a translocation
    :param current_ind: The index in all_vcf of the current variant
    :param repeat_col: Start_Repeat for INDELs or the left flank of translocations
    :return: all_vcf with repeat annotation for the current variant
    '''
    for _, reps_row in reps_tmp_filtered.iterrows():
        if (getOverlap(current_pos,
                       [reps_row.begin, reps_row.end]) > 0):
            all_vcf.loc[current_ind, repeat_col] = reps_row['class']
            break
    return all_vcf



def overlap_indels_with_repeats(all_vcf, reps_dict):
    '''
    Loop through all structural variants to annotate whether the INDEL lies within a repeat region.
    :param all_vcf: The whole mutation dataset with structural variants for all samples
    :return: 
    '''

    tmp = all_vcf[(all_vcf.Type != 'Translocation')]
   

    # loop over all INDELs
    for _, tmpr in tmp.iterrows():
        current_pos = [tmpr.pos, tmpr.end_pos]
        current_chr = tmpr.chr
        current_ind = tmpr.name

        reps_tmp = reps_dict[current_chr]

        reps_tmp_filtered = reps_tmp[(reps_tmp.begin <= current_pos[1]) & (reps_tmp.end >= current_pos[0])]
        all_vcf = loop_repeat_overlap(all_vcf=all_vcf,
                                      reps_tmp_filtered=reps_tmp_filtered,
                                      current_pos=current_pos,
                                      current_ind=current_ind,
                                      repeat_col='Start_Repeat')

    return all_vcf



def overlap_translocation_with_repeats(all_vcf, reps):
    '''
    Loop through all structural variants to annotate whether a break point from a translocation lies within a repeat region.
    :param all_vcf: The whole mutation dataset with structural variants for all samples
    :return: 
    '''
    all_vcf['Start_Repeat'] = np.nan
    all_vcf['End_Repeat'] = np.nan
    tmp = all_vcf[(all_vcf.Type == 'Translocation')]
    
    # loop over all translocations
    for _,tmprow in tmp.iterrows():
        # structural variants sometimes have a uncertain break point region and get called with a confidence interval
        # add the confidence interval to each side of the current position to get a more stringent filter setting
        region_around_translocation = tmprow.Max_CI 
        current_start_pos = [tmprow.pos - region_around_translocation, tmprow.pos + region_around_translocation]
        current_end_pos = [tmprow.end_pos - region_around_translocation, tmprow.end_pos + region_around_translocation]
        current_start_chr = tmprow.chr
        current_end_chr = tmprow.end_chr
        current_ind = tmprow.name

        # loop over all potentially overlapping repeat regions for the left flank of the translocation
        reps_tmp = reps_dict[current_start_chr]
        reps_tmp_filtered = reps_tmp[(reps_tmp.begin <= current_start_pos[1]) & (
                reps_tmp.end >= current_start_pos[0])]
        all_vcf = loop_repeat_overlap(all_vcf=all_vcf, 
                              reps_tmp_filtered=reps_tmp_filtered, 
                              current_pos=current_start_pos, 
                              current_ind=current_ind, 
                              repeat_col='Start_Repeat')

        # loop over all potentially overlapping repeat regions for the right flank of the translocation
        reps_tmp = reps_dict[current_end_chr]
        reps_tmp_filtered = reps_tmp[(reps_tmp.begin <= current_end_pos[1]) & (reps_tmp.end >= current_end_pos[0])]
        all_vcf = loop_repeat_overlap(all_vcf=all_vcf, 
                              reps_tmp_filtered=reps_tmp_filtered, 
                              current_pos=current_end_pos, 
                              current_ind=current_ind, 
                              repeat_col='End_Repeat')

    return all_vcf



def process_all_vcfs(all_vcf, spath, reps_dict):
    '''
    Add homology length information and confidence intervals around break points. Mark background mutations
    and mutations within repeat elements.
    :param all_vcf: The whole mutation dataset with structural variants for all samples. DataFrame.
    :param spath: 
    :param reps_dict: Dictionary with all ce11 repeat regions
    :return: The updated all_vcf DataFrame
    '''
    all_vcf.File = all_vcf.File.str.split('.', expand=True)[0]

    # get the homology length information
    all_vcf['HOMLEN'] = all_vcf['info'].str.split(';HOMLEN=', expand=True)[1].str.split(';', expand=True)[0].values
    all_vcf.loc[all_vcf.HOMLEN.isna(), 'HOMLEN'] = 0
    all_vcf.HOMLEN = all_vcf.HOMLEN.astype(int)
    
    # get the confidence intervals around the break points
    all_vcf['Start_Pos_CI_Start'] = \
    all_vcf['info'].str.split(';CIPOS=', expand=True)[1].str.split(';', expand=True)[0].str.split(',', expand=True)[
        0].values
    all_vcf['Start_Pos_CI_End'] = \
    all_vcf['info'].str.split(';CIPOS=', expand=True)[1].str.split(';', expand=True)[0].str.split(',', expand=True)[
        1].values
    all_vcf['End_Pos_CI_Start'] = \
    all_vcf['info'].str.split(';CIEND=', expand=True)[1].str.split(';', expand=True)[0].str.split(',', expand=True)[
        0].values
    all_vcf['End_Pos_CI_End'] = \
    all_vcf['info'].str.split(';CIEND=', expand=True)[1].str.split(';', expand=True)[0].str.split(',', expand=True)[
        1].values
    all_vcf = make_int(all_vcf, 'Start_Pos_CI_Start')
    all_vcf = make_int(all_vcf, 'Start_Pos_CI_End')
    all_vcf = make_int(all_vcf, 'End_Pos_CI_Start')
    all_vcf = make_int(all_vcf, 'End_Pos_CI_End')

    # add replication data
    sinfo = pd.read_csv(f'{spath}sample_info_withMapping.csv', index_col=0)
    sinfo = sinfo[['Rep', 'Mapped_Reads']]
    sinfo.index = sinfo.index.astype(str)
    all_vcf = all_vcf.set_index('File')
    all_vcf = all_vcf.join(sinfo)

    all_vcf = all_vcf.reset_index()
    # get all Translocation types
    all_vcf['Translocation_Type'] = np.nan
    all_vcf['Insert'] = np.nan
    for i in range(len(all_vcf)):
        translocation_type, insert = get_translocation_type(all_vcf.alt.iloc[i])
        if translocation_type == -1:  # no translocation
            continue
        all_vcf.loc[i, 'Translocation_Type'] = translocation_type
        all_vcf.loc[i, 'Insert'] = insert

    all_vcf['Max_CI'] = abs \
        (all_vcf[['Start_Pos_CI_Start', 'Start_Pos_CI_End', 'End_Pos_CI_Start', 'End_Pos_CI_End']]).max(axis=1)

    P0 = all_vcf.loc[~all_vcf.Treatment.str.contains('F1')]
    P0 = P0.sort_values('pos')
    P0 = P0.drop_duplicates(subset=['chr', 'pos', 'end_pos'], keep='last')
    F1 = all_vcf.loc[all_vcf.Treatment.str.contains('F1')]
    P0_indels = P0[P0.Type != 'Translocation']
    F1_indels = F1[F1.Type != 'Translocation']
    P0_trans = P0[P0.Type == 'Translocation']
    F1_trans = F1[F1.Type == 'Translocation']

    # mark an INDEL if the F1 indel is overlapping with any of the P0 indels, regardless of strain
    all_vcf = mark_indel_background(P0_indels, F1_indels, all_vcf)
    # mark a Translocation if both breakpoints of the F1 Translocation are overlapping with any of the P0 Translocations, 
    # regardless of strain
    all_vcf = mark_translocation_background(P0_trans, F1_trans, all_vcf)

    # overlap with known repeats
    all_vcf = overlap_translocation_with_repeats(all_vcf, reps_dict)
    all_vcf = overlap_indels_with_repeats(all_vcf, reps_dict)
    # remove Translocations that have the same repeat class at the start and the end break point
    all_vcf = all_vcf[all_vcf.Start_Repeat != all_vcf.End_Repeat]

    # mark the background translocations in the treated P0 as well
    P0_treat = all_vcf[(all_vcf.Type == 'Translocation') & (all_vcf.Treatment.str.contains('treat')) &
                       (all_vcf.Treatment.str.contains('P0'))]
    P0_ctrl = all_vcf[(all_vcf.Type == 'Translocation') & (~all_vcf.Treatment.str.contains('treat')) &
                      (all_vcf.Treatment.str.contains('P0'))]
    all_vcf = mark_translocation_background(P0_ctrl, P0_treat, all_vcf)

    # filter INDELs that lie within a repeat --> Start_Repeat in this case
    all_vcf = all_vcf[
        ((all_vcf.Type != 'Translocation') & (all_vcf.Start_Repeat.isna())) | (all_vcf.Type == 'Translocation')]

    all_vcf.to_csv(f'{spath}all_vcf_with_Repeat_and_Background.csv')
    return all_vcf


def get_repeat_data(repeat_data_file):
    '''
    
    :param repeat_data_file: File containing all repeat elements of the ce11 genome. 
    Download repeats from https://hgdownload.soe.ucsc.edu/goldenPath/ce11/bigZips/chromOut.tar.gz
    The file contains 4 columns: chromosome, begin, end, class
    With class indicating the repeat class
    :return: A dictionary containing the repeat information split into the 6 chromosomes.
    '''
    # 
    reps = pd.read_csv(repeat_data_file)
    # save the repeat data in a dictionary for speed-up
    reps_dict = {}
    reps_dict['I'] = reps[reps.chromosome == 'chrI']
    reps_dict['II'] = reps[reps.chromosome == 'chrII']
    reps_dict['III'] = reps[reps.chromosome == 'chrIII']
    reps_dict['IV'] = reps[reps.chromosome == 'chrIV']
    reps_dict['V'] = reps[reps.chromosome == 'chrV']
    reps_dict['X'] = reps[reps.chromosome == 'chrX']
    
    return reps_dict



if __name__ == "__main__":

	vcf_path = ''
	all_vcf = start_analysis(vcf_path)

	reps_dict = get_repeat_data(f'{vcf_path}all_repeats_ce11.csv')
	all_vcf = process_all_vcfs(all_vcf, vcf_path, reps_dict)




