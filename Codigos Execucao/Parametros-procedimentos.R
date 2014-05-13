setwd("~/Dropbox/R/NPFDA")
load(".RData")


## Carrega as bibliotecas necessárias
library(termstrc)
source("Experimentos reproduziveis/biblioteca-diebold-li.R")
source("Experimentos reproduziveis/biblioteca-npfda-ettj.R")
source("Experimentos reproduziveis/biblioteca simulacoes.R")
source("biblioteca-npfda.r")

# O vetor maturidade contém o tamanho, em anos, das maturidades com os meses abaixo:
# tempo.maturidades <- c(0.25 , 0.5,  0.75,  1,  1.25,  1.5,  1.75,	2,	2.5,	3,	4,
#                        5,	6,	7,	8,	9,	10,	12, 15, 20) ## Usar para Base nefasta
tempo.maturidades <- c(0.25 , 0.5,  0.75,  1,  1.25,  1.5,  1.75,	2,	2.5,	3,	4,
                       5,	6,	7,	8,	9,	10) ## Usar para Base mensal


############## Escolher UMA das BASES DE DADO abaixo: ######################
## 1.a) BASE NEFASTA: A base de dados pode ter todas as observações ...
# taxas.juro <- as.matrix(read.csv("Dados/base_nefasta_mat.csv", sep=",")
#                         [, 2:(length(tempo.maturidades)+1)])
# datas <- as.Date(as.character(read.csv("Dados/base_nefasta_mat.csv", sep=",")
#                               [ , 1]), format = "%Y-%m-%d")
# ## 1.b)  ... ou conter as observações de 5 em 5 dias
# taxas.juro <- as.matrix(read.csv("Dados/base_nefasta_mat.csv", sep=",")
#                         [1:1307*5, 2:(length(tempo.maturidades)+1)])
# datas <- as.Date(as.character(read.csv("Dados/base_nefasta_mat.csv", sep=",")
#                               [1:1307*5 , 1]), format = "%Y-%m-%d")
## 2) Base menor
taxas.juro <- as.matrix(read.table("~/Dropbox/R/NPFDA/Dados/dados_mensais.csv", 
                                   header=T, quote="\""))
datas <- as.Date(as.character(rownames(taxas.juro)), format = "%Y%m%d")
##################################

# intervalos: uma lista com os intervalos nos quais se fará a análise
# vetor.intervalos <- list(c(1,120),c(11,130),c(201,1000),c(301,1100),
#                          c(401,1200),c(501,1300)) # para base reduzida
# vetor.intervalos <- list(c(1,120),c(11,130),c(21,140)) # para base reduzida
tamanho.janela <- c(120) # tamanho da janela rolante para fazer a estimação
limites.estimacao <- c(169,349) # o início e fim da primeira observação da janela rolante. 

# METODOS ABAIXO SAO OS ESTIMADOS ATÉ AGORA
# metodos.estimar <- c("valores.reais","rw","crt.p5","ar","crt.p5.rmat6","crt.d1.rmat6",
#                      "crt.p5.rb1","crt.d1.rb1","crt.p5.rmat17", "crt.p1.rpmat",
#                      "crt.d1.rmat17","crt.p5.rpmat","crt.d1.rpmat","diebold.02","diebold.05")
metodos.estimar <- c(
"ftsa.p5_rw.rmean", "ftsa.p5_ets.rmat7", "ftsa.p3_ets.rmat17", 
"ftsa.p5_ets.rmat1", "ftsa.p5_rw.rmat7", "ftsa.p5_ets.rb1",
"crt.d1.rmat1", "crt.p5.rmat1", "crt.p5.rmean", "crt.d1.rmean", "crt.d1"
)

vetor.knn  <- 1:10*2
vetor.retirar <- list("rpmat", "rb1", "rmat9", "")
vetor.pca <- c(1,2,5)
metodos.estimar <- NULL
for (i in vetor.knn) {
  for (j in vetor.retirar) {
    for (k in vetor.pca) {
      metodos.estimar <- c(metodos.estimar,str_c("crt.k",i,"_",k,".",j, sep=""))
    }
  }
}
(metodos.estimar <- sub(pattern="[.]$",replacement="",x=metodos.estimar))


# vetor.maturidade: contém os índices das maturidades que serão estimadas 
# e previstas
vetor.maturidade <- 1:17

# vetor.horizontes
vetor.horizontes <- c(1,3,6,12)

# janela movel (moving) ou expandindo (expanding)
window <- "moving" #window pode ser "expanding" ou "moving"



############# EXECUTAR AS ESTIMAÇÕES --------------------------------------------------
#  antes verificar o arquivo "Procedimento-estimacao.r"!!!!

#### Carregar os dados
# simulacoes <- vector("list", 0) # RESETA os dados de simulacoes
load("expanding_window.dad")
load("moving_window.dad") 

# load("simulacoes_reduzida_base_mensal2.dad")

####  simulacoes <- nome.variavel ####
if (window == "moving") {
  simulacoes <- moving; print("moving")
} else if (window == "expanding") {
  simulacoes <- expanding; print("expanding")
}


source("Procedimento-estimacao.r")

# salva o objeto simulações no lugar desejado. Aconselha-se salvar com
# frequência para não correr o risco de sobrescrever os experimentos e perder
# todo o trabalho

simulacoes[["metodos"]] <- union(simulacoes[["metodos"]],metodos.estimar)
simulacoes[["horizontes"]]<- union(simulacoes[["horizontes"]],vetor.horizontes)
simulacoes[["maturidades"]] <- union(simulacoes[["maturidades"]],vetor.maturidade)
simulacoes[["janela"]] <- window

# nome.variavel <- simulacoes
if (window == "moving") {
  moving  <- simulacoes
} else if (window == "expanding") {
  expanding <- simulacoes
} else {
  print("O valor de window é inválido")
}

knn.movel <- simulacoes

# save(knn.movel,file="simulacoes_knn")
#save(moving,file="moving_window.dad")
#save(expanding,file="expanding_window.dad")
