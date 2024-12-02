	import pandas as pd
	import pickle
	import numpy as np
	from sklearn.model_selection import train_test_split
	from sklearn.preprocessing import StandardScaler
	from tensorflow import keras
	from keras.utils import to_categorical #one-hot encode target column

def loaddeepdata(species1, species2, n_classes, bin_num,
				 species1cov, species1rna, species1mods,
				 species2cov, species2rna, species2mods):
	#want to add gene name capture

	if species1cov:
		data = pickle.load(open(species1cov,"rb"))
		data = pd.DataFrame(data)
		for label,row in data.iterrows():
			data.loc[label,'Gene'] = row['Gene'].split(':')[1]
        
		#join to match up genes
		df=pd.read_csv(species1rna,sep='\t', header=0)

		species1data = pd.merge(data, df, how="left", left_on='Gene', right_on = 'Gene') #**

		X = pd.DataFrame(species1data.loc[:,species1mods])

		y = species1data.loc[:,["Average"]]        
		X1 = X
		y1 = y

		## Species 2
		data = pickle.load(open(species2cov,"rb"))
		data = pd.DataFrame(data)
		for label,row in data.iterrows():
			data.loc[label,'Gene'] = row['Gene'].split(':')[1]
        
		#join to match up genes
		df=pd.read_csv(species2rna,sep='\t', header=0)

		species2data = pd.merge(data, df, how="left", left_on='Gene', right_on = 'Gene') #**

		X = pd.DataFrame(species2data.loc[:,species2mods])

		y = species2data.loc[:,["Average"]]        
		X2 = X
		y2 = y
	else: 
		if 'Ncrassa' in (species1, species2):
			#path to all species pickle coverage files. Here using most recent 5kb files
			data = pickle.load(open("./data/Ncrassa_promoter5.0kb_genes_coverage20230718.p","rb"))
			# data = pickle.load(open("/projects/wg-feeds/ldw_downloads/Ncrassa_promoter40.0kb_genes_coverage20230718.p","rb"))
			data = pd.DataFrame(data)
			for label,row in data.iterrows():
				data.loc[label,'Gene'] = row['Gene'].split(':')[1]
			
			#join to match up genes
			df=pd.read_csv('./data/NcrassaRNACountsNorm.txt',sep='\t', header=0)

			Ncrassadata = pd.merge(data, df, how="left", left_on='Gene', right_on = 'Gene') #**

			#features: H3K27me2/3=14, H3K36me3=09, H3K4me1=06, H3K4me2=07, H3K4me3=08, H3K9me3=10

			if 'Fram' in (species1, species2):
				X = pd.DataFrame(Ncrassadata.loc[:,["SRR12229314_signal","SRR12229309_signal",
									"SRR12229307_signal", "SRR12229308_signal"]])
			elif 'LmacL' in (species1, species2):
				X = pd.DataFrame(Ncrassadata.loc[:,["SRR12229314_signal",
								"SRR12229307_signal", "SRR12229310_signal"]])  
			elif 'Anid' in (species1, species2):
				X = pd.DataFrame(Ncrassadata.loc[:,["SRR12229308_signal", "SRR12229309_signal", "SRR12229310_signal"]])
			else: 
				X = pd.DataFrame(Ncrassadata.loc[:,["SRR12229314_signal","SRR12229309_signal",
									"SRR12229307_signal", "SRR12229308_signal", "SRR12229306_signal", "SRR12229310_signal"]])

			y = Ncrassadata.loc[:,["Average"]]        
			if species1 == 'Ncrassa':
				X1 = X
				y1 = y
			else:
				X2 = X
				y2 = y

		if 'Fgram' in (species1, species2):
			Fgramdata = pickle.load(open("./data/Fgram_promoter5.0kb_alldepths_genesOnly_coverage20230722.p", "rb"))
			Fgramdata = pd.DataFrame(Fgramdata)
			for label,row in Fgramdata.iterrows():
				Fgramdata.loc[label,'Gene'] = row['Gene'].split(':')[1]

			df=pd.read_csv("./data/FgramRNACountsNorm.txt", sep='\t', header=0)
			FgramSigExprData = pd.merge(Fgramdata, df, how="left", left_on='Gene', right_on = 'Gene') #**
			#features: H3K27me3=08, H3K36me3=11, H3K4me2=13, H3K4me3=17

			if 'Ncrassa' in (species1, species2):
				X = pd.DataFrame(FgramSigExprData.loc[:,["SRR999608_signal", "SRR999611_signal",
											"SRR999613_signal","SRR999617_signal"]])
			elif 'LmacL' in (species1, species2):
				X = pd.DataFrame(FgramSigExprData.loc[:,["SRR999608_signal", "SRR999613_signal"]])
			elif 'Anid' in (species1, species2):
				X = pd.DataFrame(FgramSigExprData.loc[:,["SRR999611_signal","SRR999617_signal"]])
			else:
				X = pd.DataFrame(FgramSigExprData.loc[:,["SRR999608_signal", "SRR999611_signal",
											"SRR999613_signal","SRR999617_signal"]])
			
			y = FgramSigExprData.loc[:,["Average"]]

			if species1 == 'Fgram':
				X1 = X
				y1 = y
			else:
				X2 = X
				y2 = y

		if  'LmacL' in (species1, species2): 
			Lmacepidata = pickle.load(open("./data/Lmac_promoter5.0kb_alldepths_genesOnly_coverage20240103.p","rb"))
			Lmacepidata = pd.DataFrame(Lmacepidata)
			for label,row in Lmacepidata.iterrows():
				Lmacepidata.loc[label,'Gene'] = row['Gene'].split(':')[1]
				
			df = pd.read_csv("./data/LmacRNACountsNorm.txt",sep = '\t', header=0)
			LmacSigExprData = pd.merge(Lmacepidata, df, how="left", left_on='Gene', right_on = 'Gene') #**
			
			if 'Ncrassa' in (species1, species2): 
				X = pd.DataFrame(LmacSigExprData.loc[:,['H3K27me3_LmlA_signal', "H3K4me2_LmlA_signal",
														"H3K9me3_LmlA_signal"]])
			elif 'Fgram' in (species1, species2):
				X = pd.DataFrame(LmacSigExprData.loc[:,['H3K27me3_LmlA_signal', "H3K4me2_LmlA_signal"]])		
			elif 'Anid' in (species1, species2):
				X = pd.DataFrame(LmacSigExprData.loc[:,['H3K9me3_LmlA_signal']])
			else:
				X = pd.DataFrame(LmacSigExprData.loc[:,['H3K27me3_LmlA_signal', "H3K4me2_LmlA_signal",
														"H3K9me3_LmlA_signal"]])

			y = LmacSigExprData.loc[:,["Average"]]
			if species1 == 'LmacL':
				X1 = X
				y1 = y
			else:
				X2 = X
				y2 = y

		if 'Anid' in (species1, species2):
			data = pickle.load(open("./data/Anidulans_promoter5.0kbcoverage20240413.p","rb"))
			# data = pickle.load(open("/projects/wg-feeds/Anidulans/Anid_promoter40.0kb_alldepths_genesOnly_coverage20240628.p","rb"))
			
			data = pd.DataFrame(data)
			for label,row in data.iterrows():
				data.loc[label,'Gene'] = row['Gene'].split(':')[1]
			
			#join to match up genes
			df=pd.read_csv('./data/AnidulansRNACountsNorm.txt',sep='\t', header=0)

			Aniddata = pd.merge(data, df, how="left", left_on='Gene', right_on = 'Gene') #**
			#features: SRR2170287=H3K36me3, SRR2170288=H3K4me3, SRR2170289=H3K9me3

			if 'Ncrassa' in (species1, species2):
				X = pd.DataFrame(Aniddata.loc[:,["H3K4me3_signal", "H3K36me3_signal", "H3K9me3_signal"]])
			elif 'Fgram' in (species1, species2):
				X = pd.DataFrame(Aniddata.loc[:,["H3K36me3_signal","H3K4me3_signal"]])
			elif 'LmacL' in (species1, species2):
				X = pd.DataFrame(Aniddata.loc[:,["H3K9me3_signal"]])   
			else:
				X = pd.DataFrame(Aniddata.loc[:,["H3K4me3_signal", "H3K36me3_signal", "H3K9me3_signal"]])
			
			y = Aniddata.loc[:,["Average"]]        
			if species1 == 'Anid':
				X1 = X
				y1 = y
			else:
				X2 = X
				y2 = y

	#Data1 processing
	channel_count = X1.shape[1]
	print(f'this is the channel count: {channel_count}')
	
	X1.reset_index(inplace=True, drop=True)
	y1.reset_index(inplace=True, drop=True)

	a = X1[X1.columns[0]]
	idx = [i for i in range(len(a)) if len(a.iloc[i]) != bin_num]
	X1.drop(idx, inplace=True)
	y1.drop(idx, inplace=True)
    
	X_train1, X_test1, y_train1, y_test1 = train_test_split(X1, y1, test_size=0.2, random_state=1)
    
	X_train1 = X_train1.explode(X_train1.columns.to_list(), ignore_index=True).to_numpy()
	X_train1 = np.reshape(X_train1, (-1, bin_num, channel_count))
	X_test1 = X_test1.explode(X_test1.columns.to_list(), ignore_index=True).to_numpy()
	X_test1 = np.reshape(X_test1, (-1, bin_num, channel_count))

	scaler = StandardScaler()

	X_train1 = scaler.fit_transform(X_train1.reshape(-1, X_train1.shape[-1])).reshape(X_train1.shape)
	X_test1 = scaler.transform(X_test1.reshape(-1, X_test1.shape[-1])).reshape(X_test1.shape)

	y_train1 = pd.qcut(y_train1.iloc[:,0], q=n_classes, labels = np.arange(n_classes))
	y_test1 = pd.qcut(y_test1.iloc[:,0], q=n_classes, labels = np.arange(n_classes))
	y_train1 = to_categorical(y_train1, num_classes = n_classes,  dtype = "float32")
	y_test1 = to_categorical(y_test1, num_classes = n_classes,  dtype = "float32")

	## Data2 processing ##
	if species2 == "":
		X_train2 = []
		X_test2 = []
		y_train2 = []
		y_test2 = []
	elif species2 != "":
		X2.reset_index(inplace=True, drop=True)
		y2.reset_index(inplace=True, drop=True)
        
		a = X2[X2.columns[0]]
		idx = [i for i in range(len(a)) if len(a.iloc[i]) != bin_num]        
        
		X2.drop(idx, inplace=True)
		y2.drop(idx, inplace=True)
		
		X_train2, X_test2, y_train2, y_test2 = train_test_split(X2, y2, test_size=0.75, random_state=1)
		print(f'species2: {species2}')
		X_train2 = X_train2.explode(X_train2.columns.to_list(), ignore_index=True).to_numpy()
		
		X_train2 = np.reshape(X_train2, (-1, bin_num, channel_count))
		X_test2 = X_test2.explode(X_test2.columns.to_list(), ignore_index=True).to_numpy()
		X_test2 = np.reshape(X_test2, (-1, bin_num, channel_count))

		scaler = StandardScaler()    
		X_train2 = scaler.fit_transform(X_train2.reshape(-1, X_train2.shape[-1])).reshape(X_train2.shape)
		X_test2 = scaler.transform(X_test2.reshape(-1, X_test2.shape[-1])).reshape(X_test2.shape)
		
		y_train2 = pd.qcut(y_train2.iloc[:,0], q=n_classes, labels = np.arange(n_classes))
		y_test2 = pd.qcut(y_test2.iloc[:,0], q=n_classes, labels = np.arange(n_classes))
		y_train2 = to_categorical(y_train2, num_classes = n_classes,  dtype = "float32")
		y_test2 = to_categorical(y_test2, num_classes = n_classes,  dtype = "float32")

	return X_train1, X_test1, y_train1, y_test1, X_train2, X_test2, y_train2, y_test2

