import numpy as np
import pandas as pd
from IPython.display import display

resfinder = pd.read_table('C:/Users/diogo/Desktop/ARG_search/WasteWater/resfinder/SRR1237782_rf/ResFinder_results_tab.txt')
card = pd.read_table('C:/Users/diogo/Desktop/ARG_search/WasteWater/card_rgi/SRR1237782_card/SRR1237782_card.txt')

# resfinder
resfinder_clean = resfinder.drop(columns=['Alignment Length/Gene Length', 'Position in reference'])
resfinder_new = resfinder_clean.rename(columns={'Accession no.': 'Accession_rf', 'Position in contig': 'Position_in_contig',
                                                'Resistance gene': 'Resistance_gene', 'Coverage': 'Nucleotide_Coverage'})
resfinder_standard = resfinder_new.assign(AMR_Gene_Family=np.nan, Accession_card=np.nan, Resistance_Mechanism=np.nan, Protein_Coverage=np.nan)
resfinder_final = resfinder_standard.reindex(sorted(resfinder_standard.columns), axis=1)

# card
card_clean = card.drop(columns=['ORF_ID', 'Orientation', 'Cut_Off', 'Pass_Bitscore', 'Best_Hit_Bitscore', 'Model_type', 'SNPs_in_Best_Hit_ARO',
                                'Other_SNPs', 'Predicted_DNA', 'Predicted_Protein', 'CARD_Protein_Sequence',
                                'ID', 'Model_ID', 'Nudged', 'Note'])
card_new = card_clean.rename(columns={'AMR Gene Family': 'AMR_Gene_Family', 'Best_Hit_ARO': 'Resistance_gene', 'Drug Class': 'Phenotype',
                                      'Best_Identities': 'Identity', 'ARO': 'Accession_card', 'Resistance Mechanism': 'Resistance_Mechanism',
                                      'Percentage Length of Reference Sequence': 'Protein_Coverage'})
card_standard = card_new.assign(Accession_rf=np.nan, Nucleotide_Coverage=np.nan)
card_last = card_standard.reindex(sorted(card_standard.columns), axis=1)
card_last["Position_in_contig"] = card_last.Start.astype(str).str.cat(card_last.Stop.astype(str), sep="..")
card_final = card_last.drop(columns=["Start", "Stop"])

final_df = pd.concat([resfinder_final, card_final])
#print(final_df)
header=list(final_df.columns)
#print(header)
final_df.to_csv(r'/Users/diogo/Desktop/final_df.csv', columns=header, index=False)
display(final_df)




