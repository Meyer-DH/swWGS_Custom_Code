import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np

def getOverlap(a,b):
    return max(0,min(a[1],b[1]) - max(a[0], b[0])+1)

def mark_translocation_dups(filtered_vcf):
    '''
    filtered_vcf contains breakpoints, so not directly translocations.
    Breakpoints are often duplicated (called in both directions).
    If a breakpoint has 
        - the same start and end as the reverse of an other break point 
        - or the direction is the same and the breakpoints lie within the confidence intervals.
        - or one breakpoint overlaps, while the other is within the same repeat class. 
    it is the same translocation.
    :param filtered_vcf: 
    :return: dup_index, a list with the indeces of all duplicated translocation break points.
    '''
    filtered_vcf = filtered_vcf[filtered_vcf.Type == 'Translocation']
    filtered_vcf = filtered_vcf.sort_values(by='Translocation_Type')
    dup_index = []
    # loop over all translocations
    for _, tmp in filtered_vcf.iterrows():
        if tmp.name in dup_index:
            continue
        region_around_translocation = tmp.Max_CI
        current_start_pos = [tmp.pos - region_around_translocation, 
                             tmp.pos + region_around_translocation]
        current_end_pos = [tmp.end_pos - region_around_translocation, 
                           tmp.end_pos + region_around_translocation]
        current_start_chr = tmp.chr
        current_end_chr = tmp.end_chr
        current_ind = tmp.name
        current_startrepeat = tmp.Start_Repeat
        current_endrepeat = tmp.End_Repeat

        # only filter translocations from the same File, i.e. 'index'
        # check whether it is the same translocation just annotated the other way around, i.e.
        # start position == end position
        rev_tmp = filtered_vcf[(filtered_vcf.chr == current_end_chr) &
                               (filtered_vcf.end_chr == current_start_chr) & 
                               (filtered_vcf['index'] == tmp['index'])]
        dups = []
        for _, rev_tmpr in rev_tmp.iterrows():
            transloc_area_rev = rev_tmpr.Max_CI
            rev_tmpr_start_pos = [rev_tmpr.pos - transloc_area_rev, rev_tmpr.pos + transloc_area_rev]
            rev_tmpr_end_pos = [rev_tmpr.end_pos - transloc_area_rev, rev_tmpr.end_pos + transloc_area_rev]
            
            # if the left break point == right break point (within the region given by Max_CI) filter
            if ((getOverlap(current_start_pos, rev_tmpr_end_pos) > 0) &
                    (getOverlap(current_end_pos, rev_tmpr_start_pos) > 0) &
                    (current_start_chr == rev_tmpr.end_chr) &
                    (current_end_chr == rev_tmpr.chr) & 
                    (tmp['index'] == rev_tmpr['index'])):
                dups.append(rev_tmpr.name)

        # it can also happen that the same translocation is called twice in the same orientation
        # this is mostly the case for translocation with a high area of uncertainty, i.e. Max_CI
        forw_tmp = filtered_vcf[(filtered_vcf.chr == current_start_chr) &
                                (filtered_vcf.end_chr == current_end_chr) & 
                                (filtered_vcf['index'] == tmp['index'])]
        for _, forw_tmpr in forw_tmp.iterrows():
            transloc_area_rev = forw_tmpr.Max_CI
            forw_tmpr_start_pos = [forw_tmpr.pos - transloc_area_rev, forw_tmpr.pos + transloc_area_rev]
            forw_tmpr_end_pos = [forw_tmpr.end_pos - transloc_area_rev, forw_tmpr.end_pos + transloc_area_rev]
            forw_tmpr_name = forw_tmpr.name
            if ((getOverlap(current_start_pos, forw_tmpr_start_pos) > 0) &
                    (getOverlap(current_end_pos, forw_tmpr_end_pos) > 0) &
                    (current_start_chr == forw_tmpr.chr) &
                    (current_end_chr == forw_tmpr.end_chr) & 
                    (tmp['index'] == forw_tmpr['index']) & 
                    (current_ind != forw_tmpr_name)):
                dups.append(forw_tmpr_name)

        # The above still does not filter if the start (or end) position is the same, while the other position is different,
        # but lies within in the same repeat class, this is probably just a mapping problem and artifact.
        # Filter out those translocations, that have the same start pos and end repeat, respective end pos and start repeat
        forw_tmp = filtered_vcf[(filtered_vcf['index'] == tmp['index'])]
        for _, forw_tmpr in forw_tmp.iterrows():
            transloc_area_rev = forw_tmpr.Max_CI

            forw_tmpr_start_pos = [forw_tmpr.pos - transloc_area_rev, forw_tmpr.pos + transloc_area_rev]
            forw_tmpr_end_pos = [forw_tmpr.end_pos - transloc_area_rev, forw_tmpr.end_pos + transloc_area_rev]
            forw_tmpr_End_Repeat = forw_tmpr.End_Repeat
            forw_tmpr_Start_Repeat = forw_tmpr.Start_Repeat
            forw_tmpr_chr = forw_tmpr.chr
            forw_tmpr_end_chr = forw_tmpr.end_chr
            forw_tmpr_name = forw_tmpr.name
            # Start positions(left flanks) the same, end positions different, but within the same repeat class
            if ((getOverlap(current_start_pos, forw_tmpr_start_pos) > 0)&
                    (current_start_chr == forw_tmpr_chr)&
                    (current_endrepeat == forw_tmpr_End_Repeat)&
                    (tmp['index'] == forw_tmpr['index'])&
                    (current_ind != forw_tmpr_name)&
                    (forw_tmpr_name not in dups)):
                dups.append(forw_tmpr_name)
            # End positions(right flanks) the same, start positions different, but within the same repeat class
            elif ((getOverlap(current_end_pos, forw_tmpr_end_pos) > 0)&
                  (current_end_chr == forw_tmpr_end_chr)&
                  (current_startrepeat == forw_tmpr_Start_Repeat)&
                  (tmp['index'] == forw_tmpr['index'])&
                  (current_ind != forw_tmpr_name) & 
                  (forw_tmpr_name not in dups)):
                dups.append(forw_tmpr_name)
            # 2 break points overlap, while the other 2 are within the same repeat class
            elif ((getOverlap(current_start_pos, forw_tmpr_end_pos) > 0) &
                  (current_start_chr == forw_tmpr_end_chr)&
                  (current_endrepeat == forw_tmpr_Start_Repeat)&
                  (tmp['index'] == forw_tmpr['index'])&
                  (current_ind != forw_tmpr_name)&
                  (forw_tmpr_name not in dups)):
                dups.append(forw_tmpr_name)
            elif ((getOverlap(current_end_pos, forw_tmpr_start_pos) > 0)&
                  (current_end_chr == forw_tmpr_chr)&
                  (current_startrepeat == forw_tmpr_End_Repeat)&
                  (tmp['index'] == forw_tmpr['index'])&
                  (current_ind != forw_tmpr_name)&
                  (forw_tmpr_name not in dups)):
                dups.append(forw_tmpr_name)

        if len(dups) > 0:
            dup_index += dups
    return dup_index


def filter_all_vcf(all_vcf, outdir):
    # only get novel F1 mutations that are not appearing in P0
    filtered_vcf = all_vcf[(all_vcf.Treatment.str.contains('F1')) &
                           (all_vcf.Marked.isna())]

    # focus on heterozygous mutations
    filtered_vcf = filtered_vcf[(filtered_vcf.Zygosity == 'Heterozygous')]

    meta = filtered_vcf.Treatment.str.split('_', expand=True)
    meta.columns = ['Generation', 'Strain', 'Treated', 'Sex']
    filtered_vcf = filtered_vcf.join(meta)
    filtered_vcf.to_csv(f'{outdir}all_vcf_F1_Heterozygous_filtered.csv')

    # filter out duplicate breakpoints
    dup_index = mark_translocation_dups(filtered_vcf)
    filtered_vcf = filtered_vcf.loc[~filtered_vcf.index.isin(dup_index)]
    # add Parent and Replicate information
    filtered_vcf['Parent'] = filtered_vcf.Rep.str.split('-', expand=True)[0]
    filtered_vcf['Replicate'] = filtered_vcf.Rep.str.split('-', expand=True)[1]


    filtered_vcf=filtered_vcf[filtered_vcf.Type=='Translocation']

    filtered_vcf.to_csv(f'{outdir}Translocations_F1_Heterozygous_filtered_withoutBreakPointDups.csv')
    return filtered_vcf



def count_mutations(tmp_df, count_name = 'Number of new F1 mutations'):
    '''
    Count the amount of translocations per file.
    :param tmp_df: a filtered_vcf DataFrame
    :param count_name: How to call the new column with the count information
    :return: The DataFrame mcounts with the structural variant counts 
    (split by Type, i.e. Insertion, Deletion, Tandem Duplication, Translocation)
    '''
    # count mutations
    mcounts = \
        tmp_df.groupby(['Strain', 'Treated', 'Sex', 'Type', 'Rep', 'index', 'Mapped_Reads', 'Zygosity']).count()['chr']
    mcounts = mcounts.reset_index()
    mcounts['Parent'] = mcounts.Rep.str.split('-', expand=True)[0]
    mcounts['Name'] = mcounts.Strain + '_' + mcounts.Sex
    mcounts.columns = ['Strain', 'Treated', 'Sex', 'Type', 'Rep', 'File', 'Mapped_Reads', 'Zygosity',
                       count_name, 'Parent', 'Name']

    mcounts['Log_Mapped_Reads'] = np.log(mcounts.Mapped_Reads)

    return mcounts


def plot_mcounts(tmp_df, outdir,order,new_labels, count_name = 'Number of new F1 mutations', fileend='pdf', w=0.5):
    ax = sns.boxplot(x='Name', 
                     y=count_name, 
                     data=tmp_df,
                     color='white', 
                     width=w)
    ax = sns.swarmplot(x='Name', 
                       y=count_name, 
                       data=tmp_df, 
                       order=order, 
                       dodge=True)
    ax.set_xticklabels(new_labels)
    plt.xlabel('Strain')
    plt.tight_layout()
    plt.savefig(f'{outdir}Translocations_{count_name}.{fileend}')
    plt.close()


# count mutations, correct for mapping rate, and plot
def count_plot(tmp_df, outdir, order, new_labels):
    count_name = 'Number of new F1 translocations'


    mcounts = count_mutations(tmp_df, count_name=count_name)
    mcounts_translocations = mcounts[mcounts.Type == 'Translocation']
    
    w = 0.25 * len(tmp_df.Strain.unique())  # adjust the size/width of the boxes according to the number of boxes
    plot_mcounts(mcounts_translocations, outdir=outdir, order=order, new_labels=new_labels,
                 count_name=count_name, w=w)

    mcounts.to_csv(f'{outdir}Translocations_{count_name}_mcounts.csv')

    return 


# load ce11.fa
def load_fasta_into_dict(fasta_file = 'ce11_genome.fa'):
    fasta = {}
    with open(fasta_file, 'r') as file_one:
        for line in file_one:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                chr_name = line[1:]
                if chr_name not in fasta:
                    fasta[chr_name] = []
                continue
            sequence = line
            fasta[chr_name].append(sequence)

    for ch in fasta:
        fasta[ch] = ''.join(fasta[ch])
    return fasta


def complement_seq(seq):
    reverse_dict = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
    return ''.join([reverse_dict[s] for s in seq])




def calc_microhomology(vcf_df, fasta, region_size=8):
    '''
    Get region around start and end of the break points and build a 2D array with zeros and add +1 every time the base is the same
    :return: A region_size*region_size matrix around the break points. Each bin containing the percentage base similarity (microhomology) pooled over all translocations in vcf_df.
    '''
    # will contain the raw data for the heatmap
    array_microhomology = np.zeros((region_size, region_size))

    # loop over all translocations
    for _,vcf_df_row in vcf_df.iterrows():
        start_chr = vcf_df_row.chr
        start_pos = vcf_df_row.pos
        end_chr = vcf_df_row.end_chr
        end_pos = vcf_df_row.end_pos
        translocation_type = vcf_df_row.Translocation_Type

        # compare the left and right flank dependent on the translocation type
        if translocation_type == 1:
            left_flank = ''.join(reversed(fasta[start_chr][start_pos - (int(region_size / 2) - 1):start_pos + int(
                region_size / 2) + 1]))  # 3 left from pos, 4 right
            right_flank = fasta[end_chr][end_pos - int(region_size / 2):end_pos + int(
                region_size / 2)]  # 4 left from the break point and 3 right
        elif translocation_type == 2:
            left_flank = ''.join(reversed(fasta[start_chr][start_pos - (int(region_size / 2) - 1):start_pos + int(
                region_size / 2) + 1]))  # 3 left from pos, 4 right
            right_flank = complement_seq(''.join(reversed(fasta[end_chr][end_pos - (int(region_size / 2) - 1):end_pos + int(
                region_size / 2) + 1])))  # 3 left from pos, 4 right and REVERSED
        elif translocation_type == 3:
            left_flank = ''.join(reversed(fasta[start_chr][start_pos - (int(region_size / 2) - 1):start_pos + int(
                region_size / 2) + 1]))  # 3 left from pos, 4 right
            right_flank = fasta[end_chr][end_pos - int(region_size / 2):end_pos + int(
                region_size / 2)]  # 4 left from the break point and 3 right
        elif translocation_type == 4:
            left_flank =  complement_seq(fasta[start_chr][start_pos - int(region_size / 2):start_pos + int(
                region_size / 2)])  # left is reversed and complemented, and then reversed again for the plot
            right_flank =fasta[end_chr][end_pos - int(region_size / 2):end_pos + int(
                region_size / 2)] # right is normal

        for i, bl in enumerate(left_flank): 
            for j, br in enumerate(right_flank):
                if bl == br:
                    array_microhomology[i][j] += 1
    # divide the count array by the total amount of translocations to get a value between 0 and 1
    return array_microhomology/len(vcf_df)



def calc_random_microhomology(vcf_df, fasta, region_size=8):
    '''
    Same as calc_microhomology() but for a random permutation.
    Get region around start and end of the random break points and build a 2D array with zeros and add +1 every time the base is the same
    :return: A region_size*region_size matrix around the random break points. Each bin containing the percentage base similarity (microhomology) pooled over all random regions
    '''
    array_microhomology = np.zeros((region_size, region_size))
    chroms = ['I', 'II', 'III', 'IV', 'V', 'X']
    chrom_len_dict_ce11 = {'I': 15072434,
                           'II': 15279421,
                           'III': 13783801,
                           'IV': 17493829,
                           'V': 20924180,
                           'X': 17718942
                           }
    translocation_type_counts = vcf_df.Translocation_Type.value_counts()
    translocation_type_list = []
    for tc in translocation_type_counts.index:
        translocation_type_list += [int(tc)] * translocation_type_counts[tc]

    for t in range(len(vcf_df)):
        start_chr = chroms[np.random.randint(6)]
        start_pos = np.random.randint(chrom_len_dict_ce11[start_chr])
        end_chr = chroms[np.random.randint(6)]
        end_pos = np.random.randint(chrom_len_dict_ce11[end_chr])
        translocation_type = translocation_type_list[t]

        if translocation_type == 1:
            left_flank = ''.join(reversed(fasta[start_chr][start_pos - (int(region_size / 2) - 1):start_pos + int(
                region_size / 2) + 1]))
            right_flank = fasta[end_chr][end_pos - int(region_size / 2):end_pos + int(
                region_size / 2)]
        elif translocation_type == 2:
            left_flank = ''.join(reversed(fasta[start_chr][start_pos - (int(region_size / 2) - 1):start_pos + int(
                region_size / 2) + 1]))
            right_flank = complement_seq(
                ''.join(reversed(fasta[end_chr][end_pos - (int(region_size / 2) - 1):end_pos + int(
                    region_size / 2) + 1])))
        elif translocation_type == 3:
            left_flank = ''.join(reversed(fasta[start_chr][start_pos - (int(region_size / 2) - 1):start_pos + int(
                region_size / 2) + 1]))
            right_flank = fasta[end_chr][end_pos - int(region_size / 2):end_pos + int(
                region_size / 2)]
        elif translocation_type == 4:
            left_flank = complement_seq(fasta[start_chr][start_pos - int(region_size / 2):start_pos + int(
                region_size / 2)])
            right_flank = fasta[end_chr][end_pos - int(region_size / 2):end_pos + int(
                region_size / 2)]

        for i, bl in enumerate(left_flank):
            for j, br in enumerate(right_flank):
                if bl == br:
                    array_microhomology[i][j] += 1
    return array_microhomology / len(vcf_df)

# plot heatmap
def plotheat(microdf,outfile, cmap='Greys', vmin=0, maxi=0.4, fontsize=10):
    region_size = int(len(microdf)/2)
    ax = sns.heatmap(microdf, linewidth=0.5, cmap=cmap, vmin=vmin, vmax=maxi)
    labels = list(range(-region_size, 0))+list(range(region_size))
    ax.set_xticklabels(map(str, labels), fontsize=fontsize)
    ax.set_yticklabels(map(str, labels), fontsize=fontsize)
    plt.tight_layout()
    plt.savefig(outfile)
    plt.close()


def calc_microhomologies_and_plot(tmp_df, outname, region_size,fasta,  outdir,fileend, maxi):
    '''
    Calculate the microhomology of the translocations in tmp_df, do the same for a random set of the same length and plot both
    :param tmp_df: a filtered_vcf df
    :param outname: 
    :return: 
    '''

    array_microhomology_all = calc_microhomology(tmp_df,
                                                                 fasta,
                                                                 region_size=region_size)
    array_microhomology_random = calc_random_microhomology(tmp_df,
                                                                           fasta,
                                                                           region_size=region_size)
    plotheat(microdf=array_microhomology_all,
             outfile=f'{outdir}microhomology_{outname}.{fileend}',
             cmap='Greys',
             vmin=0,
             maxi=maxi)
    plotheat(microdf=array_microhomology_random,
             outfile=f'{outdir}microhomology_{outname}_random.{fileend}',
             cmap='Greys',
             vmin=0,
             maxi=maxi)



def loop_microhom(tmp_df, region_size,fasta,  outdir,fileend, maxi, outn=''):
    '''
    
    :param tmp_df: 
    :param region_size: 
    :param fasta: 
    :param outdir: 
    :param fileend: 
    :param maxi: 
    :param outn: 
    :return: 
    '''
    # Only use translocations without an insertion
    tmp_df = tmp_df.drop_duplicates(subset=['chr', 'pos', 'end_chr', 'end_pos'], keep='first')
    tmp_df = tmp_df[(tmp_df.Type == 'Translocation') & (tmp_df.Insert.isna())]
    
    # calculate and plot microhomologies for Type 1,2, and 4
    calc_microhomologies_and_plot(tmp_df=tmp_df[(tmp_df.Translocation_Type==1.0)],
                                  outname=outn+'all_type1', region_size=region_size,fasta=fasta,  outdir=outdir,fileend=fileend, maxi=maxi)
    calc_microhomologies_and_plot(tmp_df=tmp_df[(tmp_df.Translocation_Type==2.0)],
                                  outname=outn+'all_type2', region_size=region_size,fasta=fasta,  outdir=outdir,fileend=fileend, maxi=maxi)
    calc_microhomologies_and_plot(tmp_df=tmp_df[(tmp_df.Translocation_Type==4.0)],
                                  outname=outn+'all_type4', region_size=region_size,fasta=fasta,  outdir=outdir,fileend=fileend, maxi=maxi)


def find_template(tmp, regio=50, insert_size=3):
    '''
    Search the inserted sequence within +-regio of the breakpoints, additionally search for the reverse, complement and reverse complement sequence.
    :param tmp: A filtered_vcf DataFrame with Translocations containing an insertion between their breakpoints
    :param regio: region around the breakpoints to search for the template
    :param insert_size: Minimum insert size to consider, a small insertion will be found by pure chance
    :return: The same DataFrame with additional columns
    '''
    tmp = tmp[tmp.Insert.str.len() >= insert_size]
    print(f'>={insert_size} len(tmp)')
    tmp['Start_Normal'] = np.nan
    tmp['Start_Reversed'] = np.nan
    tmp['Start_Complement'] = np.nan
    tmp['Start_Reverse_Complement'] = np.nan
    tmp['End_Normal'] = np.nan
    tmp['End_Reversed'] = np.nan
    tmp['End_Complement'] = np.nan
    tmp['End_Reverse_Complement'] = np.nan

    tmp = tmp.reset_index()
    for i, tmp_row in tmp.iterrows():

        # search the left breakpoint
        fasta_region = fasta[tmp_row.chr][tmp_row.pos - regio:tmp_row.pos + regio]
        start_search = (fasta_region.find(tmp_row.Insert),
                        fasta_region.find(''.join(reversed(tmp_row.Insert))),
                        fasta_region.find(complement_seq(tmp_row.Insert)),
                        fasta_region.find(''.join(reversed(complement_seq(tmp_row.Insert)))))
        if start_search[0] > -1:
            tmp.loc[i, 'Start_Normal'] = start_search[0] - regio
        if start_search[1] > -1:
            tmp.loc[i, 'Start_Reversed'] = start_search[1] - regio
        if start_search[2] > -1:
            tmp.loc[i, 'Start_Complement'] = start_search[2] - regio
        if start_search[3] > -1:
            tmp.loc[i, 'Start_Reverse_Complement'] = start_search[3] - regio

        # search the right breakpoint
        fasta_region = fasta[tmp_row.end_chr][tmp_row.end_pos - regio:tmp_row.end_pos + regio]
        end_search = (fasta_region.find(tmp_row.Insert),
                      fasta_region.find(''.join(reversed(tmp_row.Insert))),
                      fasta_region.find(complement_seq(tmp_row.Insert)),
                      fasta_region.find(''.join(reversed(complement_seq(tmp_row.Insert)))))
        if end_search[0] > -1:
            tmp.loc[i, 'End_Normal'] = end_search[0] - regio
        if end_search[1] > -1:
            tmp.loc[i, 'End_Reversed'] = end_search[1] - regio
        if end_search[2] > -1:
            tmp.loc[i, 'End_Complement'] = end_search[2] - regio
        if end_search[3] > -1:
            tmp.loc[i, 'End_Reverse_Complement'] = end_search[3] - regio

    return tmp


def count_inserted_template_translocations(tmp_df, insert_region_limits=25, insert_size=3):
    '''
    Count translocation with and without inserts and those that have a templated one within a certain region

    '''
    translocations_with_inserts_df = tmp_df[(tmp_df.Type == 'Translocation') &
                                            (~tmp_df.Insert.isna()) ]

    total_number_of_translocations_with_inserts = len(translocations_with_inserts_df)
    total_number_of_translocations_without_inserts = len(tmp_df[(tmp_df.Type == 'Translocation') &
                                                                (tmp_df.Insert.isna())])
    total_number_of_translocations = len(tmp_df[(tmp_df.Type == 'Translocation')])

    translocations_with_inserts_bigger_insertsize_df = find_template(translocations_with_inserts_df,
                                                                     regio=insert_region_limits,
                                                                     insert_size=insert_size)
    translocations_with_inserts_of_size = len(translocations_with_inserts_bigger_insertsize_df)

    translocations_with_inserts_bigger_insertsize_df[
        'InsertPosition'] = translocations_with_inserts_bigger_insertsize_df.Start_Normal.astype(
        str) + '_' + translocations_with_inserts_bigger_insertsize_df.Start_Reversed.astype(
        str) + '_' + translocations_with_inserts_bigger_insertsize_df.Start_Complement.astype(
        str) + '_' + translocations_with_inserts_bigger_insertsize_df.Start_Reverse_Complement.astype(
        str) + '_' + translocations_with_inserts_bigger_insertsize_df.End_Normal.astype(
        str) + '_' + translocations_with_inserts_bigger_insertsize_df.End_Reversed.astype(
        str) + '_' + translocations_with_inserts_bigger_insertsize_df.End_Complement.astype(
        str) + '_' + translocations_with_inserts_bigger_insertsize_df.End_Reverse_Complement.astype(str)
    translocations_with_inserts_bigger_insertsize_df.InsertPosition = translocations_with_inserts_bigger_insertsize_df.InsertPosition.str.replace(
        'nan_', '')
    translocations_with_inserts_bigger_insertsize_df.InsertPosition = translocations_with_inserts_bigger_insertsize_df.InsertPosition.str.replace(
        '_nan', '')
    translocations_with_inserts_bigger_insertsize_df.InsertPosition = translocations_with_inserts_bigger_insertsize_df.InsertPosition.str.replace(
        'nan', '')
    translocations_with_templated_inserts = len(translocations_with_inserts_bigger_insertsize_df[
                                                    translocations_with_inserts_bigger_insertsize_df.InsertPosition != ''])
    print(translocations_with_templated_inserts)
    insertpos = translocations_with_inserts_bigger_insertsize_df.InsertPosition.str.split('_', expand=True)[0].values
    insertpos = np.append(insertpos,
                          translocations_with_inserts_bigger_insertsize_df.InsertPosition.str.split('_', expand=True)[
                              1].values)

    insertposition = []
    for i in insertpos:
        if i != '' and i != None:
            insertposition.append(int(float(i)))

    miscellaneous_inserts_of_translocations = total_number_of_translocations_with_inserts - translocations_with_templated_inserts
    data = [total_number_of_translocations_without_inserts, miscellaneous_inserts_of_translocations,
            translocations_with_templated_inserts, total_number_of_translocations]
    return insertposition, data


def plot_inserted_translocations(insertposition, data, insert_region_limits, outdir, outname=None, fileend='pdf'):
    ax = sns.kdeplot(insertposition)
    ax = sns.rugplot(insertposition)
    ax = plt.axvline(-insert_region_limits, 0, 1, color='black')
    ax = plt.axvline(insert_region_limits, 0, 1, color='black')
    plt.title(
        f'Inserts of size >= {insert_size}bp\nwithin +-{insert_region_limits}bp of the translocation break points\n')
    plt.xlabel('Distance from break points in bp')
    plt.tight_layout()
    plt.savefig(f'{outdir}Translocations_Templated_Insertions_Region_Distribution_{outname}.{fileend}')
    plt.close()

    labels = [f'No inserts\n(n={data[0]})',
              f'Miscellaneous inserts\n(n={data[1]})',
              f'Templated inserts\n(n={data[2]})']
    colors = sns.color_palette('colorblind')[0:3]
    plt.pie(data[:3], labels=labels, autopct='%.0f%%', colors=colors)
    plt.title(f'Translocations n={data[3]}')
    plt.tight_layout()
    plt.savefig(f'{outdir}Translocations_Insertions_Distribution_Pie_{outname}.{fileend}')
    plt.close()



def calc_inserted_and_plot(tmp_df, insert_region_limits, insert_size, outdir, outname=None, fileend='pdf'):

	insertposition, data = count_inserted_template_translocations(tmp_df=tmp_df, 
		                                                      insert_region_limits=insert_region_limits, 
		                                                      insert_size=insert_size)

	plot_inserted_translocations(insertposition=insertposition, 
		                     data=data, 
		                     insert_region_limits=insert_region_limits, 
		                     outdir=outdir, 
		                     outname=outname, fileend=fileend)





if __name__ == "__main__":
	# load the result of the preprocessing script
	vcf_path = ''
	all_vcf = pd.read_csv('all_vcf_with_Repeat_and_Background.csv', index_col=0)

	filtered_vcf = filter_all_vcf(all_vcf, vcf_path)



	# generate count plot
	order=['N2_hermaphrodite', 'fog-2_female', 'fog-2_male']
	new_labels = ['WT', '$\it{fog}$-$\it{2}$ female', '$\it{fog}$-$\it{2}$ male']
	count_plot(tmp_df=filtered_vcf, outdir=vcf_path, order=order, new_labels=new_labels)
	# only fog-2
	order=['fog-2_female', 'fog-2_male']
	new_labels = ['$\it{fog}$-$\it{2}$ female', '$\it{fog}$-$\it{2}$ male']
	count_plot(tmp_df=filtered_vcf[filtered_vcf.Strain=='fog-2'], outdir=vcf_path+'fog2_', order=order, new_labels=new_labels)
	# only WT
	order=['N2_hermaphrodite']
	new_labels = ['WT']
	count_plot(tmp_df=filtered_vcf[filtered_vcf.Strain=='N2'], outdir=vcf_path+'WT_', order=order, new_labels=new_labels)



	#Microhomology heatmap
	fasta = load_fasta_into_dict()
	fileend='pdf'
	maxi=1
	region_size=8
	loop_microhom(tmp_df=filtered_vcf[filtered_vcf.Strain=='N2'], 
		      region_size=region_size,
		      fasta=fasta,  
		      outdir=vcf_path,
		      fileend=fileend,
		      maxi=maxi, 
		      outn='F1_N2')
	loop_microhom(tmp_df=filtered_vcf[filtered_vcf.Strain=='fog-2'], 
		      region_size=region_size,
		      fasta=fasta, 
		      outdir=vcf_path,
		      fileend=fileend, 
		      maxi=maxi, 
		      outn='F1_fog2')


	# Templated Insertions
	insert_region_limits = 25
	insert_size = 3
	calc_inserted_and_plot(tmp_df=filtered_vcf[filtered_vcf.Strain=='N2'], 
		               insert_region_limits=insert_region_limits, 
		               insert_size=insert_size, 
		               outdir=vcf_path, 
		               outname='F1_N2', fileend=fileend)
	calc_inserted_and_plot(tmp_df=filtered_vcf[filtered_vcf.Strain=='fog-2'], 
		               insert_region_limits=insert_region_limits, 
		               insert_size=insert_size, 
		               outdir=vcf_path, 
		               outname='F1_fog2', fileend=fileend)



