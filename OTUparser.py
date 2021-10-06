#/root/miniconda3/bin/python
import subprocess
import pandas as pd
import scipy
from skbio.diversity import alpha_diversity
import seaborn as sns
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser(description='Calculate alpha and beta diversities from normalized or unormalized OTU tables')
parser.add_argument("-i", "--input", type=str,help='OTU table with taxa by row and sample by column extension .tsv')
parser.add_argument("-m", "--method", type=str, help='Type "cum" for cumulative sum scalling normalization')
parser.add_argument("-d", "--distance", type=str,help='Use normalized "norm" or unormalized "un" OTU table')
args = parser.parse_args()


def calculate_normalized_alpha():
	df=pd.read_csv("normalized.tsv", sep="\t")
	obs_otus=pd.DataFrame(alpha_diversity('observed_otus', df.transpose(), df.columns), columns=['Observed OTUs'])
	chao1=pd.DataFrame(alpha_diversity('chao1', df.transpose(),df.columns), columns=['Chao1'])
	simpson=pd.DataFrame(alpha_diversity('simpson',df.transpose(),df.columns), columns=['Simpson'])
	shannon=pd.DataFrame(alpha_diversity('shannon', df.transpose(),df.columns), columns=['Shannon'])
	alpha_div=pd.concat(([obs_otus, chao1, simpson, shannon]), axis=1, join='inner', sort=True).reset_index()
	return alpha_div.to_csv("alpha_diversity_normalized.tsv", index=False, sep="\t", header=True)

def plot_normalized():
	normalized=pd.read_csv("alpha_diversity_normalized.tsv", sep="\t", index_col=0)
	fig,ax=plt.subplots(2,2)
	for i,f in enumerate(list(normalized.columns.values)):
		graph=normalized.boxplot(f, ax=ax.flatten()[i])
	return  plt.savefig("alpha_metrics_normalized.pdf")

def calculate_unormalized_alpha():
	print("Skip normalization","\n","Calculating alpha-diversity indexes")
	df=pd.read_csv(args.input, sep="\t",  index_col=0)
	obs_otus=pd.DataFrame(alpha_diversity('observed_otus', df.transpose(), df.columns), columns=['Observed OTUs'])
	chao1=pd.DataFrame(alpha_diversity('chao1', df.transpose(),df.columns), columns=['Chao1'])
	simpson=pd.DataFrame(alpha_diversity('simpson',df.transpose(),df.columns), columns=['Simpson'])
	shannon=pd.DataFrame(alpha_diversity('shannon', df.transpose(),df.columns), columns=['Shannon'])
	alpha_div=pd.concat(([obs_otus, chao1, simpson, shannon]), axis=1, join='inner', sort=True).reset_index()
	return alpha_div.to_csv("alpha_diversity_unormalized.tsv", index=False, sep="\t", header=True)

def plot_unormalized():
	unormalized=pd.read_csv("alpha_diversity_unormalized.tsv", sep="\t",  index_col=0)
	fig,ax=plt.subplots(2,2)
	for i,f in enumerate(list(unormalized.columns.values)):
		graph=unormalized.boxplot(f, ax=ax.flatten()[i])
	return plt.savefig("alpha_metrics_unormalized.pdf")


def plot_beta_div():
	if ((args.distance=='norm') and args.method=='cum'):
		beta=pd.read_csv("normalized.tsv", sep="\t")
		chartBeta=sns.clustermap(beta)
		return plt.savefig("Beta_diversity_normalized.pdf")
	elif ((args.distance=='un') and args.method==None):
		beta=pd.read_csv(args.input, sep="\t",  index_col=0)
		chartBeta=sns.clustermap(beta)
		return plt.savefig("Beta_diversity_unormalized.pdf")

if args.method =='cum':
	subprocess.call(["Rscript", "/root/workshop_UNICAMP/normalize.R", "-i", args.input])
	calculate_normalized_alpha()
	plot_normalized()
elif args.method == None:
	calculate_unormalized_alpha()
	plot_unormalized()

if args.distance is not None:
	plot_beta_div()
