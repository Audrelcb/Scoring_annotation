
## open sirius file compound_identifications.tsv and add the suffix _sirius to each column of the file
df = pd.read_csv(input_sirius_path, sep='\t') 
suffix = "_sirius"

#adding _sirius suffix
for col in df.columns:
    if col != "id":
        new_col = col + suffix
        df.rename(columns={col: new_col}, inplace=True)

#for each element in column id, keep only the characters after the last '_'
df['id'] = df['id'].str.split('_').str[-1]

#Clean sirius file and save
# List of columns to keep
columns_to_keep = ['id', 'ConfidenceScore_sirius', 'CSI:FingerIDScore_sirius', 'ZodiacScore_sirius', 'SiriusScore_sirius', 'molecularFormula_sirius', 'adduct_sirius', 'InChIkey2D_sirius', 'InChI_sirius', 'name_sirius', 'smiles_sirius', 'xlogp_sirius', 'ionMass_sirius']

# Drop unwanted columns
columns_to_drop = [col for col in df.columns if col not in columns_to_keep]
df = df.drop(columns_to_drop, axis=1)

#save the new dataframe as a tsv file under the name compound_identifications_sirius_cleaned.tsv
df.to_csv(output_sirius_path, sep='\t', index=False)

# %%
#open the met_annot_enhancer file 
df2 = pd.read_csv(r'C:\Users\Lecab\gitrepo\plantes_ornementales\ornamental\cdhrc_dbgi\pos\results\met_annot_enhancer\cdhrc_dbgi_alc_pos\cdhrc_dbgi_alc_pos_spectral_match_results_repond.tsv', sep='\t')
suffix = "_metAnnot"

for col in df2.columns:
    if col != "feature_id":
        new_col = col + suffix
        df2.rename(columns={col: new_col}, inplace=True)

#save the new dataframe as a tsv file under the name compound_identifications_cleaned.tsv
df.to_csv(r'C:\Users\Lecab\gitrepo\plantes_ornementales\ornamental\cdhrc_dbgi\pos\results\met_annot_enhancer\cdhrc_dbgi_alc_pos\cdhrc_dbgi_alc_pos_spectral_match_results_repond_cleaned.tsv', sep='\t', index=False)
df2.head()

# %%
#open canopus_formula_summary.tsv 
df3 = pd.read_csv(r'C:\Users\Lecab\gitrepo\plantes_ornementales\ornamental\cdhrc_dbgi\pos\results\sirius\canopus_formula_summary.tsv', sep='\t')
suffix = "_canopus"
for col in df3.columns:
    if col != "id":
        new_col = col + suffix
        df3.rename(columns={col: new_col}, inplace=True)

# List of columns to keep
columns_to_keep = ['id', 'molecularFormula_canopus', 'adduct_canopus', 'NPC#pathway_canopus', 'NPC#superclass_canopus', 'NPC#class_canopus']

# Drop unwanted columns
columns_to_drop = [col for col in df3.columns if col not in columns_to_keep]
df3 = df3.drop(columns_to_drop, axis=1)
df3.head()

#for each element in column id, keep only the characters after the last '_'
df3['id'] = df3['id'].str.split('_').str[-1]

#save the new dataframe as a tsv file under the name canopus_formula_summary_cleaned.tsv
df3.to_csv(r'C:\Users\Lecab\gitrepo\plantes_ornementales\ornamental\cdhrc_dbgi\pos\results\sirius\canopus_formula_summary_cleaned.tsv', sep='\t', index=False)

# %%
## MERGING
# Read the SIRIUS input files
df1 = pd.read_csv(output_sirius_path, sep='\t')

# Merge the SIRIUS df with MetAnnot df based on the "id" or "feature_id" column
merged_df = pd.merge(df1, df2, left_on='id', right_on='feature_id', how='outer')

# Conserved only one column with the ID 
merged_df['merge_id'] = merged_df['id'].fillna(merged_df['feature_id'])

merged_df = merged_df.drop(['id', 'feature_id'], axis=1)

merged_df['merge_id'] = merged_df['merge_id'].astype(int)

# Saved the merged dataframe in a new file 
merged_df.to_csv(r'C:\Users\Lecab\gitrepo\plantes_ornementales\ornamental\cdhrc_dbgi\pos\results\tmp\merged_met_annot_sirius.tsv', sep='\t', index=False)  # Replace with the name of you output file 

df3 = pd.read_csv(r'C:\Users\Lecab\gitrepo\plantes_ornementales\ornamental\cdhrc_dbgi\pos\results\sirius\canopus_formula_summary_cleaned.tsv', sep='\t')
df4 = pd.read_csv(r'C:\Users\Lecab\gitrepo\plantes_ornementales\ornamental\cdhrc_dbgi\pos\results\tmp\merged_met_annot_sirius.tsv', sep='\t')  # Replace with the name of you output file 
df4.head()

# Merge the dataframes based on the "id" or "feature_id" column
merged_df = pd.merge(df3, df4, left_on='id', right_on='merge_id', how='outer')

# Conserved only one column with the ID 
merged_df['shared_id'] = merged_df['id'].fillna(merged_df['merge_id'])

merged_df = merged_df.drop(['id', 'merge_id'], axis=1)

merged_df['shared_id'] = merged_df['shared_id'].astype(int)

# Saved the merged dataframe in a new file 
merged_df.to_csv(r'C:\Users\Lecab\gitrepo\plantes_ornementales\ornamental\cdhrc_dbgi\pos\results\tmp\merged_met_annot_sirius_canopus.tsv', sep='\t', index=False)  # Replace with the name of you output file 

df5 = pd.read_csv(r'C:\Users\Lecab\gitrepo\plantes_ornementales\ornamental\cdhrc_dbgi\pos\results\tmp\merged_met_annot_sirius_canopus.tsv', sep='\t')
# %%
#open DB_result corresponding to your gnps job (csv)
merged_GNPS_df = pd.read_csv(r'C:\Users\Lecab\gitrepo\plantes_ornementales\ornamental\cdhrc_dbgi\pos\results\met_annot_enhancer\a346106e484e4c7da5f194621bd4f2d3\DB_result\e548fcbbaa394b6b95b8069f4314a6fc.tsv', sep='\t') #change the path to your respective file
suffix = "_GNPS"

for col in merged_GNPS_df.columns:
    if col != "#Scan#":
        new_col = col + suffix
        merged_GNPS_df.rename(columns={col: new_col}, inplace=True)
merged_GNPS_df.head()

# %%
#convert the dataframe csv to tsv
tsv_data = merged_GNPS_df.to_csv(sep='\t', index=False)
merged_GNPS_df = pd.read_csv(io.StringIO(tsv_data), sep='\t')

merged_GNPS_df.head()

# %%
# Read the merged GNPS with MetAnnot_sirius_canopus files dataframe 
df1 = pd.read_csv(r'C:\Users\Lecab\gitrepo\plantes_ornementales\ornamental\cdhrc_dbgi\pos\results\tmp\merged_met_annot_sirius_canopus.tsv', sep='\t')
df1.head()

# %%
# Merge the dataframes based on the "id" or "feature_id" column
merged_df = pd.merge(merged_GNPS_df, df1, left_on='#Scan#', right_on='shared_id', how='outer')

# Conserved only one column with the ID 
merged_df['ID'] = merged_df['#Scan#'].fillna(merged_df['shared_id'])

merged_df = merged_df.drop(['#Scan#', 'shared_id'], axis=1)

merged_df['ID'] = merged_df['ID'].astype(int)

# Saved the merged dataframe in a new file 
merged_df.to_csv(r'C:\Users\Lecab\gitrepo\plantes_ornementales\ornamental\cdhrc_dbgi\pos\results\tmp\merged_met_annot_sirius_canopus_GNPS.tsv', sep='\t', index=False)  # Replace with the name of you output file 


# %%
input_file = r'C:\Users\Lecab\gitrepo\plantes_ornementales\ornamental\cdhrc_dbgi\pos\results\tmp\merged_met_annot_sirius_canopus_GNPS.tsv'

# Charger uniquement les colonnes spécifiées dans un DataFrame
df6 = pd.read_csv(input_file, sep='\t')

# Afficher le DataFrame résultant
df6.head()


# %%
#Define the input file path
input_file = r'C:\Users\Lecab\gitrepo\plantes_ornementales\ornamental\cdhrc_dbgi\pos\results\tmp\merged_met_annot_sirius_canopus_GNPS.tsv'

# Define the column names
inchikey2d_column = 'InChIkey2D_sirius'
inchikey_planar_column = 'InChIKey-Planar_GNPS'
short_inchikey_column = 'short_inchikey_metAnnot'

# Define the scoring function
def calculate_score(inchikey2d_sirius, inchikey_planar_gnps, short_inchikey_metannot):
    # If all three InChIKeys are nan, empty or null, return 0
    if pd.isnull(inchikey2d_sirius) and pd.isnull(inchikey_planar_gnps) and pd.isnull(short_inchikey_metannot):
        return 0
    
    # If all three InChIKeys are similar, but not null, empty or nan, return 3
    elif inchikey2d_sirius == inchikey_planar_gnps == short_inchikey_metannot and pd.notnull(inchikey2d_sirius):
        return 3
    
    # If two of the InChIKeys are similar, but not null, empty or nan, return 2 
    elif (inchikey2d_sirius == inchikey_planar_gnps and pd.notnull(inchikey2d_sirius)) or (inchikey2d_sirius == short_inchikey_metannot and pd.notnull(inchikey2d_sirius)) or (inchikey_planar_gnps == short_inchikey_metannot and pd.notnull(inchikey_planar_gnps)):
        return 2
    
    # If at least one of the InChIKeys is not nan, empty or null, return 1
    else:
        return 1
 
    

# Apply the scoring function to create a new 'score' column
df6['score'] = df6.apply(lambda row: calculate_score(row[inchikey2d_column], row[inchikey_planar_column], row[short_inchikey_column]), axis=1)


df6.to_csv(r'C:\Users\Lecab\gitrepo\plantes_ornementales\ornamental\cdhrc_dbgi\pos\results\sirius\merged_met_annot_sirius_canopus_scre.tsv', sep='\t', index=False)  # Replace with the name of you output file 

# %%
df6.head()