
library(GEOquery)
gse <- getGEO("GSE42842")
class(gse)
show(gse)
show(pData(phenoData(gse[[1]]))[1:4,])
class(gse[[1]])

# The expression set is a first element of the list 
eset <- gse[[1]]
head(pData(eset))
class(eset)
# Save expression set
save( eset, file="GSE42842_ESet.RData" )


# Save the expression matrix
write.csv(exprs(eset), file="GSE42842_expressions.csv")

# Get sample data
write.csv(pData(eset), file="GSE42842_phenodata.csv")

# Get annotation
gpl <- getGEO(eset@annotation)
class(gpl)
head(gpl@dataTable@table)
annotation <- gpl@dataTable@table
head(annotation)[1:3,]
# Save annotation
write.csv(annotation, file="GSE42842_annotation.csv")
