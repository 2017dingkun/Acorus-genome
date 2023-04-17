import argparse
import pandas as pd
from rpy2.robjects import pandas2ri
import rpy2.robjects as robj
pandas2ri.activate()

parser = argparse.ArgumentParser()
parser.add_argument('genePair')
parser.add_argument('RsemReadsCounts')
args = parser.parse_args()

geneExpression = {}
with open(args.RsemReadsCounts, 'r') as RRC:
	for idx, line in enumerate(RRC):
		if idx == 0:
			continue
		
		linelist = line.rstrip().split('\t')
		gene, count = linelist[0], linelist[4]
		geneExpression[gene] = float(count)


expressionA = {}
expressionB = {}
with open(args.genePair, 'r') as RB, open('homolog_pair_expression.txt', 'w') as HPE:
	for line in RB:
		geneA, geneB = line.rstrip().split('\t')[:2]
		HPE.write('{}\t{}\t{}\n'.format(geneA+'#'+geneB, geneExpression[geneA], geneExpression[geneB]))
		expressionA[geneA+'#'+geneB] = geneExpression[geneA]
		expressionB[geneA+'#'+geneB] = geneExpression[geneB]

expressionData = pd.DataFrame({'A':expressionA, 'B':expressionB}) #这里必须把基因名放在行名上，不然edgeR会出问题

def EdgeR(indata):
	Rcode = '''
			library('limma')
			library('edgeR')

			NoRepDifferentailExpression <- function(data){
				group=factor(c("A","B"))
				d<-DGEList(count=data,group=group)
				d <- calcNormFactors(d,method="TMM")
				bcv <- 0.2
				et <- exactTest(d, dispersion=bcv^2)
				padj <- p.adjust(et$table$PValue,method = "BH")
				allRes <- cbind(et$table[,-2], padj)
				#b<-subset(allRes,(logFC>2|-logFC>2)&padj<0.05&PValue<0.01)
				write.table(allRes, file="diff.norep.edgeR.xls", sep="\t", quote=F,row.name=T)
				}
			'''
	data = indata[(indata['A']>0) | (indata['B']>0)]
	print(data)
	data = pandas2ri.py2rpy_pandasdataframe(indata)
	robj.r(Rcode)
	NoRepDifferentailExpression = robj.globalenv['NoRepDifferentailExpression']
	NoRepDifferentailExpression(data)

EdgeR(expressionData)
