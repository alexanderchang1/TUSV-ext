import gzip
import pandas as pd
from io import StringIO

def read_gridss_sv(input_file_path):
    """
    Reads the GRIDSS VCF file and processes the structural variants for each sample.
    """
    filtered_lines = []
    with gzip.open(input_file_path, 'rt') as file:
        for line in file:
            if not line.startswith('##'):  # Skipping the VCF header
                filtered_lines.append(line)
    
    filtered_content = ''.join(filtered_lines)
    df_gridss_sv = pd.read_csv(StringIO(filtered_content), sep='\t')
    
    # Adjust chromosome names if necessary
    df_gridss_sv['#CHROM'] = df_gridss_sv['#CHROM'].apply(convert_chromosome)

    return df_gridss_sv

def rewrite_gridss_results(df_gridss_sv, ascat_df, sample_name):
    """
    Processes the GRIDSS SV data, filters for relevant structural variants, and writes to a TUSV-ext compatible VCF.
    """
    # Filter for variants that pass all filters
    df_gridss_sv = df_gridss_sv[df_gridss_sv['FILTER'] == "PASS"]
    
    # Extract BND (Breakend) type variants
    bnd_df = df_gridss_sv[df_gridss_sv['INFO'].str.contains("SVTYPE=BND")]

    # Split the INFO column into a dictionary for easy access
    bnd_df['INFO_DICT'] = bnd_df['INFO'].apply(lambda x: dict(item.split('=') for item in x.split(';') if '=' in item))

    # Filter out unpaired BNDs based on MATEID
    mate_ids = bnd_df['INFO_DICT'].apply(lambda x: x.get('MATEID', None))
    paired_bnd_df = bnd_df[bnd_df['ID'].isin(mate_ids.values)]

    # Prepare new VCF columns
    vcf_columns = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', sample_name]
    vcf_df = pd.DataFrame(columns=vcf_columns)
    
    # Create a dictionary to store the paired rows
    paired_rows = {}
    
    # Renaming column heads for sample-specific access
    paired_bnd_df.rename(columns=lambda x: rename_column(x), inplace=True)
    
    for index, row in paired_bnd_df.iterrows():
        mate_id = row['INFO_DICT']['MATEID']
        paired_rows[row['ID']] = row
        if mate_id in paired_rows:
            mate_row = paired_rows[mate_id]

            # Create the ALT field
            alt = row["ALT"].replace('chr', '').replace(row['REF'], '').replace('X', '23')
            row_sv_ID = f"sv{sv_count}"
            increment_sv()

            mate_row_sv_ID = f"sv{sv_count}"
            increment_sv()

            # Create the INFO field
            info = f"MATEID={mate_row_sv_ID};SVTYPE=BND"

            # Create the FORMAT fields (e.g., GT, CNADJ, BDP, DP)
            format_field = "GT:CNADJ:BDP:DP"

            # Parse metrics from ascat_df or other quality metrics
            tumor_data, normal_data = calculate_gridss_metrics(ascat_df, row)

            # Append the row to the new DataFrame
            vcf_df = vcf_df.append({
                '#CHROM': row['#CHROM'], 'POS': row['POS'], 'ID': row_sv_ID, 'REF': '.',
                'ALT': alt, 'QUAL': row['QUAL'], 'FILTER': row['FILTER'], 'INFO': info,
                'FORMAT': format_field, sample_name: tumor_data
            }, ignore_index=True)

            # Process the mate pair similarly
            alt_mate = mate_row["ALT"].replace('chr', '').replace(mate_row['REF'], '').replace('X', '23')
            info_mate = f"MATEID={row_sv_ID};SVTYPE=BND"
            tumor_data_mate, normal_data_mate = calculate_metrics(ascat_df, mate_row)

            vcf_df = vcf_df.append({
                '#CHROM': mate_row['#CHROM'], 'POS': mate_row['POS'], 'ID': mate_row_sv_ID, 'REF': '.',
                'ALT': alt_mate, 'QUAL': mate_row['QUAL'], 'FILTER': mate_row['FILTER'], 'INFO': info_mate,
                'FORMAT': format_field, sample_name: tumor_data_mate
            }, ignore_index=True)

    return vcf_df

def calculate_gridss_metrics(ascat_df, gridss_row):
    """
    Calculates metrics for the TUMOR and NORMAL samples based on GRIDSS output.
    """
    data_columns = ['TUMOR', 'NORMAL']

    for data in data_columns:
        # Parse the FORMAT field to extract relevant information for the sample
        format_fields = gridss_row['FORMAT'].split(':')
        sample_data = gridss_row[data].split(':')

        format_dict = dict(zip(format_fields, sample_data))

        # Use read pair (RP) and split read (SR) counts if available
        PR_ref = int(format_dict.get('RP', '0,0').split(",")[0])
        PR_alt = int(format_dict.get('RP', '0,0').split(",")[1])
        PR_tot = PR_ref + PR_alt

        SR_ref = int(format_dict.get('SR', '0,0').split(",")[0])
        SR_alt = int(format_dict.get('SR', '0,0').split(",")[1])
        SR_tot = SR_ref + SR_alt

        # If there is allele fraction information (AF), it can be used instead
        AF = float(format_dict.get('AF', '0.0'))

        if data == 'TUMOR':
            GT = "1|1"
            CN_tot = get_total_copy_number(ascat_df, gridss_row['#CHROM'], gridss_row['POS'])
            
            if PR_tot + SR_tot > 0:
                CNADJ = (PR_alt + SR_alt) / (PR_tot + SR_tot) * CN_tot
            else:
                CNADJ = AF * CN_tot  # Use allele fraction if counts are not available

            CNADJ = float(CNADJ)
            CNADJ_formatted = format(CNADJ, '.2f')

            BDP = SR_alt
            DP = SR_alt + PR_alt
            variables = [GT, CNADJ_formatted, BDP, DP]
            tumor_data = ':'.join(map(str, variables))
        
        elif data == 'NORMAL':
            GT = "0|0"
            CN_tot = 2  # Normal copy number is 2

            if PR_tot + SR_tot > 0:
                CNADJ = (PR_ref + SR_ref) / (PR_tot + SR_tot) * CN_tot
            else:
                CNADJ = (1 - AF) * CN_tot  # Use allele fraction if counts are not available

            CNADJ = float(CNADJ)
            CNADJ_formatted = format(CNADJ, '.2f')

            BDP = SR_ref
            DP = SR_ref + PR_alt  # Alt depth from normal sample may still contribute
            variables = [GT, CNADJ_formatted, BDP, DP]
            normal_data = ':'.join(map(str, variables))
        else:
            raise Exception("Unknown Datatype")

    return tumor_data, normal_data


def init_vcf(output_folder, sample):
    """
    Initializes the VCF output file with the appropriate TUSV-ext format headers.
    """
    output_file_path = os.path.join(output_folder, sample + '.vcf')
    with open(output_file_path, 'w') as file_vcf:
        file_vcf.write('##fileformat=VCFv4.2\n')
        file_vcf.write('##filedate=20241011\n')
        file_vcf.write('##INFO=<ID=MATEID,Number=.,Type=String,Description="ID of mate breakends">\n')
        file_vcf.write('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n')
        file_vcf.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        file_vcf.write('##FORMAT=<ID=CNADJ,Number=.,Type=Float,Description="Copy number of adjacency">\n')
        file_vcf.write('##FORMAT=<ID=BDP,Number=1,Type=Integer,Description="Depth of split reads">\n')
        file_vcf.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">\n')
        file_vcf.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' + sample + '\n')
    
    return output_file_path

# Example usage:
# vcf_df = rewrite_gridss_results(df_gridss_sv, ascat_df, 'sample_1')
# output_file = init_vcf('/output/folder', 'sample_1')
