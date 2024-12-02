import conorm
import pandas as pd

def RNA_norm(rawRNAdata, normalizedRNAdata):
    df=pd.read_csv(rawRNAdata, sep='\t', header=1)

    for label,row in df.iterrows():
        df.loc[label,'Gene'] = row['Geneid'].split(':')[1]

    df=df.set_index('Gene')

    #tmm normalizing
    df_tmm = conorm.tmm(df)
    #taking the average of the normalized samples
    df_tmm['Average'] = df_tmm.mean(axis=1)
    #output txt file to call in loadshallowdata.py
    df_tmm.to_csv(normalizedRNAdata, sep='\t')
