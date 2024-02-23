# BT1013 Análisis de biología computacional
# RETO INTEGRADOR
# Libraries
library(dplyr)
library(seqinr)
library(ggplot2)
# library(ggmsa) #####INSTALAR Y PROGRAMAR URGENTE
cat("\014") 

# Constants
translate = c(UUU="FPhe", UUC="FPhe", UUA="LLeu", UUG="LLeu", UCU="SSer", 
              UCC="SSer", UCA="SSer", UCG="SSer", UAU="YTyr", UAC="YTyr", 
              UAA='-', UAG='-', UGU="CCys", UGC="CCys", UGA='-', UGG="WTrp", 
              CUU="LLeu", CUC="LLeu", CUA="LLeu", CUG="LLeu", CCU="PPro", 
              CCC="PPro", CCA="PPro", CCG="PPro", CAU="HHis", CAC="HHis", 
              CAA="QGln", CAG="QGln", CGU="RArg", CGC="RArg", CGA="RArg", 
              CGG="RArg", AUU="IIle", AUC="IIle", AUA="IIle", AUG="MMet", 
              ACU="TThr", ACC="TThr", ACA="TThr", ACG="TThr", AAU="NAsn", 
              AAC="NAsn", AAA="KLys", AAG="KLys", AGU="SSer", AGC="SSer", 
              AGA="RArg", AGG="RArg", GUU="VVal", GUC="VVal", GUA="VVal", 
              GUG="VVal", GCU="AAla", GCC="AAla", GCA="AAla", GCG="AAla", 
              GAU="DAsp", GAC="DAsp", GAA="EGlu", GAG="EGlu", GGU="GGly", 
              GGC="GGly", GGA="GGly", GGG="GGly");

cat("\014") 
# Main
# Leemos nuestras secuencias
original = read.fasta("wuHan.txt")
omicronVariant = read.fasta("omicron.txt")
deltaVariant = read.fasta("delta.txt")
# Guardamos la secuencia del gen S en un vector
dnaOriginal = toupper(original[[2]])
dnaOmicron = toupper(omicronVariant [[2]])
dnaDelta = toupper(deltaVariant [[2]])
# Convertimos de RNAm a ADN
dnaOriginal[dnaOriginal == "T"] = "U"
dnaOmicron[dnaOmicron == "T"] = "U"
dnaDelta[dnaDelta == "T"] = "U"
# Alineamos nuestras secuencias 
# Alineamos Original v Omicron
A = dnaOriginal
B = dnaOmicron
m = matrix(data=0,nrow=length(B)+1, ncol=length(A)+1)
m[1, ] = seq(0,length(A))*-2
m[ ,1] = seq(0,length(B))*-2
m

CalcularPeso = function(fila, col){
  if (A[col-1]==B[fila-1]) diagonal = m[fila-1,col-1] + 1
  else diagonal = m[fila-1,col-1] - 1
  up = m[fila-1,col] - 2
  left = m[fila,col-1] - 2
  peso = max (diagonal, up, left)
  return (peso)
}

for (fila in seq(2,length(B)+1)){
  for (col in seq(2,length(A)+1)){
    m[fila, col] = CalcularPeso(fila, col) 
  }
}
m

fila = length(B)+1
col = length(A)+1
solA = solB = c()
while(fila>1 || col>1){
  if (col>1 && fila>1 && A[col-1]==B[fila-1]){
    solA = c(A[col-1], solA)
    solB = c(B[fila-1], solB)
    col = col - 1
    fila = fila -1
  }else{
    if (col>1 && fila>1 && m[fila-1, col-1]>m[fila, col-1] && m[fila-1, col-1]>m[fila-1, col]){
      solA = c(A[col-1], solA)
      solB = c(B[fila-1], solB)
      col = col - 1
      fila = fila -1
    }else if (fila==1 || (col>1 && m[fila, col-1] > m[fila-1, col])){
      solA = c(A[col-1], solA)
      solB = c("_", solB)
      col = col - 1
    }else if (col==1 || (fila>1 && m[fila, col-1] <= m[fila-1, col])){
      solA = c("_", solA)
      solB = c(B[fila-1], solB)
      fila = fila - 1 
    }
  }
}
sol = matrix(data=c(solA,solB),nrow = 2, ncol = length(solA), byrow = TRUE)
sol
#Calculamos similitud en %
counter = length(solA)
counter2 = 0
acc = 0
for(i in 1:counter) {
  if(solA[i] == solB[i]) {
    acc = acc + 1
  }
}
acc = acc / counter
acc = acc * 100
cat("Obtuvimos una similitud del", acc, "%")

#

# INSTALAR GGMSA Y PROGRAMAR PARA GRAFICAR SECUENCIAS ALINEADAS

#


# Alineamos Original vs Delta
A = dnaOriginal
B = dnaDelta
m = matrix(data=0,nrow=length(B)+1, ncol=length(A)+1)
m[1, ] = seq(0,length(A))*-2
m[ ,1] = seq(0,length(B))*-2
m

CalcularPeso = function(fila, col){
  if (A[col-1]==B[fila-1]) diagonal = m[fila-1,col-1] + 1
  else diagonal = m[fila-1,col-1] - 1
  up = m[fila-1,col] - 2
  left = m[fila,col-1] - 2
  peso = max (diagonal, up, left)
  return (peso)
}

for (fila in seq(2,length(B)+1)){
  for (col in seq(2,length(A)+1)){
    m[fila, col] = CalcularPeso(fila, col) 
  }
}
m

fila = length(B)+1
col = length(A)+1
solA = solB = c()
while(fila>1 || col>1){
  if (col>1 && fila>1 && A[col-1]==B[fila-1]){
    solA = c(A[col-1], solA)
    solB = c(B[fila-1], solB)
    col = col - 1
    fila = fila -1
  }else{
    if (col>1 && fila>1 && m[fila-1, col-1]>m[fila, col-1] && m[fila-1, col-1]>m[fila-1, col]){
      solA = c(A[col-1], solA)
      solB = c(B[fila-1], solB)
      col = col - 1
      fila = fila -1
    }else if (fila==1 || (col>1 && m[fila, col-1] > m[fila-1, col])){
      solA = c(A[col-1], solA)
      solB = c("_", solB)
      col = col - 1
    }else if (col==1 || (fila>1 && m[fila, col-1] <= m[fila-1, col])){
      solA = c("_", solA)
      solB = c(B[fila-1], solB)
      fila = fila - 1 
    }
  }
}
sol = matrix(data=c(solA,solB),nrow = 2, ncol = length(solA), byrow = TRUE)
sol
#Calculamos similitud en %
counter = length(solA)
counter2 = 0
acc2 = 0
for(i in 1:counter) {
  if(solA[i] == solB[i]) {
    acc2 = acc2 + 1
  }
}
acc2 = acc2 / counter
acc2 = acc2 * 100
cat("Obtuvimos una similitud del", acc2, "%")

#

# INSTALAR GGMSA Y PROGRAMAR PARA GRAFICAR SECUENCIAS ALINEADAS

#


# Alineamos Delta vs Omicron
A = dnaDelta
B = dnaOmicron
m = matrix(data=0,nrow=length(B)+1, ncol=length(A)+1)
m[1, ] = seq(0,length(A))*-2
m[ ,1] = seq(0,length(B))*-2
m

CalcularPeso = function(fila, col){
  if (A[col-1]==B[fila-1]) diagonal = m[fila-1,col-1] + 1
  else diagonal = m[fila-1,col-1] - 1
  up = m[fila-1,col] - 2
  left = m[fila,col-1] - 2
  peso = max (diagonal, up, left)
  return (peso)
}

for (fila in seq(2,length(B)+1)){
  for (col in seq(2,length(A)+1)){
    m[fila, col] = CalcularPeso(fila, col) 
  }
}
m

fila = length(B)+1
col = length(A)+1
solA = solB = c()
while(fila>1 || col>1){
  if (col>1 && fila>1 && A[col-1]==B[fila-1]){
    solA = c(A[col-1], solA)
    solB = c(B[fila-1], solB)
    col = col - 1
    fila = fila -1
  }else{
    if (col>1 && fila>1 && m[fila-1, col-1]>m[fila, col-1] && m[fila-1, col-1]>m[fila-1, col]){
      solA = c(A[col-1], solA)
      solB = c(B[fila-1], solB)
      col = col - 1
      fila = fila -1
    }else if (fila==1 || (col>1 && m[fila, col-1] > m[fila-1, col])){
      solA = c(A[col-1], solA)
      solB = c("_", solB)
      col = col - 1
    }else if (col==1 || (fila>1 && m[fila, col-1] <= m[fila-1, col])){
      solA = c("_", solA)
      solB = c(B[fila-1], solB)
      fila = fila - 1 
    }
  }
}
sol = matrix(data=c(solA,solB),nrow = 2, ncol = length(solA), byrow = TRUE)
sol
#Calculamos similitud en %
counter = length(solA)
counter2 = 0
acc3 = 0
for(i in 1:counter) {
  if(solA[i] == solB[i]) {
    acc3 = acc3 + 1
  }
}
acc3 = acc3 / counter
acc3 = acc3 * 100
cat("Obtuvimos una similitud del", acc3, "%")
#

# INSTALAR GGMSA Y PROGRAMAR PARA GRAFICAR SECUENCIAS ALINEADAS

#

# Inicializamos nuestros data frames
# donde estaremos comparando y guardando nuestros aminoacidos
# 1. Original vs Omicron
dfOriginalOmicron = data.frame(
  mutation = character(), 
  nucleo = numeric(), 
  codon = character(), 
  protein = character(), 
  index = numeric()
)
dfOriginalOmicron
str(dfOriginalOmicron)
# 2. Original vs Delta
dfOriginalDelta = data.frame(
  mutation = character(), 
  nucleo = numeric(), 
  codon = character(), 
  protein = character(), 
  index = numeric()
)
dfOriginalDelta
str(dfOriginalDelta)
# 3. Delta vs Omicron
dfDeltaOmicron = data.frame(
  mutation = character(), 
  nucleo = numeric(), 
  codon = character(), 
  protein = character(), 
  index = numeric()
)
dfDeltaOmicron
str(dfDeltaOmicron)
cat("\014") 
# Comparamos nuestras variantes de interes con la original
# E incluso entre ellas
OvO = which(dnaOriginal != dnaOmicron);
OvD = which(dnaOriginal != dnaDelta);
DvO = which(dnaDelta != dnaOmicron);
# Comenzamos nuestra comparacion y almacenaje de Original vs Omicron
if(length(OvO) > 0) {
  print(paste(length(OvO), "Nucleoditos"))
  for (i in OvO) {
    mutation = paste(dnaOriginal[i], "to", dnaOmicron[i], sep = "")
    start = i - (i - 1) %% 3
    startG = 21563 + start
    index = as.integer(i/3 + 1)
    codonOriginal = paste(dnaOriginal[start], dnaOriginal[start + 1], dnaOriginal[start + 2], sep = "")
    codonOmicron = paste(dnaOmicron[start], dnaOmicron[start + 1], dnaOmicron[start + 2], sep = "")
    newCodon = paste(codonOriginal, startG, codonOmicron, sep = "")
    newAmino = paste(translate[codonOriginal], index, translate[codonOmicron], sep = "")
    print(paste(mutation, startG, newCodon, newAmino, index))
    newRow = list(mutation, startG, newCodon, newAmino, index)
    dfOriginalOmicron[nrow(dfOriginalOmicron) + 1, ] = newRow
  }
}
dfOriginalOmicron
# Comenzamos nuestra comparacion y almacenaje de Original vs Delta
if(length(OvD) > 0) {
  print(paste(length(OvD), "Nucleoditos"))
  for (i in OvD) {
    mutation = paste(dnaOriginal[i], "to", dnaDelta[i], sep = "")
    start = i - (i - 1) %% 3
    startG = 21563 + start
    index = as.integer(i/3 + 1)
    codonOriginal = paste(dnaOriginal[start], dnaOriginal[start + 1], dnaOriginal[start + 2], sep = "")
    codonDelta = paste(dnaDelta[start], dnaDelta[start + 1], dnaDelta[start + 2], sep = "")
    newCodon = paste(codonOriginal, startG, codonDelta, sep = "")
    newAmino = paste(translate[codonOriginal], index, translate[codonDelta], sep = "")
    print(paste(mutation, startG, newCodon, newAmino, index))
    newRow = list(mutation, startG, newCodon, newAmino, index)
    dfOriginalDelta[nrow(dfOriginalDelta) + 1, ] = newRow
  }
}
dfOriginalDelta
# Comenzamos nuestra comparacion y almacenaje de Delta vs Omicron
if(length(DvO) > 0) {
  print(paste(length(DvO), "Nucleoditos"))
  for (i in DvO) {
    mutation = paste(dnaDelta[i], "to", dnaOmicron[i], sep = "")
    start = i - (i - 1) %% 3
    startG = 21563 + start
    index = as.integer(i/3 + 1)
    codonDelta = paste(dnaDelta[start], dnaDelta[start + 1], dnaDelta[start + 2], sep = "")
    codonOmicron = paste(dnaOmicron[start], dnaOmicron[start + 1], dnaOmicron[start + 2], sep = "")
    newCodon = paste(codonDelta, startG, codonOmicron, sep = "")
    newAmino = paste(translate[codonDelta], index, translate[codonOmicron], sep = "")
    print(paste(mutation, startG, newCodon, newAmino, index))
    newRow = list(mutation, startG, newCodon, newAmino, index)
    dfDeltaOmicron[nrow(dfDeltaOmicron) + 1, ] = newRow
  }
}
dfDeltaOmicron

# Guardamos todo en un solo data frame para su graficacion
df = data.frame(
  mutation = character(), 
  nucleo = numeric(), 
  codon = character(), 
  protein = character(), 
  index = numeric()
)
df
str(df)

df = rbind(dfOriginalOmicron, dfOriginalDelta) 
df = rbind(df, dfDeltaOmicron)
df

# Graficamos
m = ggplot(dfOriginalOmicron)
m = m + aes(x=mutation, fill=mutation, label=..count..)
m = m + ggtitle("Mutaciones detectadas, gen. Original vs variante Omicron")
m = m + labs(x= "Mutacion", y= "Frecuencia", fill= "Mutaciones")
m = m + geom_bar(stat = "count") 
m = m + geom_text(stat = "count", vjust= 0)
m

n = ggplot(dfOriginalDelta)
n = n + aes(x=mutation, fill=mutation, label=..count..)
n = n + ggtitle("Mutaciones detectadas, gen. Original vs variante Delta")
n = n + labs(x= "Mutacion", y= "Frecuencia", fill= "Mutaciones")
n = n + geom_bar(stat = "count") 
n = n + geom_text(stat = "count", vjust= 0)
n

o = ggplot(dfDeltaOmicron)
o = o + aes(x=mutation, fill=mutation, label=..count..)
o = o + ggtitle("Mutaciones detectadas, variante Omicron vs variante Delta")
o = o + labs(x= "Mutacion", y= "Frecuencia", fill= "Mutaciones")
o = o + geom_bar(stat = "count") 
o = o + geom_text(stat = "count", vjust= 0)
o

p = ggplot(df)
p = p + aes(x=mutation, fill=mutation, label=..count..)
p = p + ggtitle("Mutaciones detectadas, gen. Original vs variante Omicron vs variante Delta")
p = p + labs(x= "Mutacion", y= "Frecuencia", fill= "Mutaciones")
p = p + geom_bar(stat = "count") 
p = p + geom_text(stat = "count", vjust= 0)
p

cat("Obtuvimos una similitud del", acc, "%")
cat("Obtuvimos una similitud del", acc2, "%")
cat("Obtuvimos una similitud del", acc3, "%")
