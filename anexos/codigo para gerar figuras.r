# retira as margens dos arquivos
ajeita.janela <- function() par(mar=c(2.1, 2.1, 1.1, 2.1))
setwd("~/Dropbox/dissertacao/anexos")


# ETTJ completa -----------------------------------------------------------

taxas.juro <- as.matrix(read.table("~/Dropbox/R/NPFDA/Dados/dados_mensais.csv", 
                                    quote="\""))
pdf(width=6,height=4,file="taxas_juro.pdf")
ajeita.janela()
ts.plot(taxas.juro,lty = (1:dim(taxas.juro)[2] - 1) %/% 8 + 1, col=1:dim(taxas.juro)[2], lwd=0.3)
dev.off()


# Bases: b-splines e fourier ----------------------------------------------

require(fda)
pdf(width=6,height=4,file="base_bsplines1.pdf")
ajeita.janela()
plot(tempbasis <- create.bspline.basis(c(0,10),nbasis=8))
dev.off()

pdf(width=6,height=4,file="base_bsplines2.pdf")
ajeita.janela()
plot(tempbasis <- create.bspline.basis(c(0,10),nbasis=6,norder=3))
dev.off()


pdf(width=6,height=4,file="base_fourier.pdf")
ajeita.janela()
plot(monthbasis <- create.fourier.basis(c(0,12),dropind=1,nbasis=5))
dev.off()


# Tabelas de resultados ---------------------------------------------------

load("~/Dropbox/R/NPFDA/moving_window.dad")
load("~/Dropbox/R/NPFDA/expanding_window.dad")
source("~/Dropbox/R/NPFDA/Experimentos reproduziveis/biblioteca simulacoes.R")
require(xtable)
tempo.maturidades <- c(0.25 , 0.5,  0.75,  1,  1.25,  1.5,  1.75,  2,	2.5,	3,	4,
                       5,	6,	7,	8,	9,	10)
require(xtable)
for (horizonte in c(1,3,6,12)) {
  nome.arquivo <- paste0("movel-hor",horizonte,".tex")
  legenda = paste("Janela móvel - horizonte", horizonte)
  tabela.moving <- TabelaEQM(moving,1:17, horizonte = horizonte,normaliza=TRUE)
  tmp <- apply(tabela.moving,c(1,2),as.character)
  dados <- apply(tmp, c(1,2), substr, start=1, stop=4)
  #tabela.gw <- GiacominiWhite(moving,benchmark="rw",tamanho=120,horizonte=horizonte,vetor.maturidade)
  tabela.dm <- DieboldMariano(moving,benchmark="rw",horizonte=horizonte,vetor.maturidade)
  signif <- t(apply(tabela.dm,c(1,2), Significancia))
  nomes.metodos <- row.names(dados)
  require(stringr)
  compilado <- str_c(dados,signif)
  tabela.final <- matrix(compilado, ncol=17)
  row.names(tabela.final) <- nomes.metodos
  colnames(tabela.final)<- tempo.maturidades
#   View(tabela.final)
#   View(dados)
#   View(signif)
  tab.latex <- xtable(tabela.final, caption=legenda)
  print(tab.latex,type = "latex",file=nome.arquivo)
  linhas.tabela <- readLines(con=nome.arquivo)
  linhas.tabela <- append(linhas.tabela, values = "\\scalebox{0.7}{", after=4)
  linhas.tabela <- append(linhas.tabela, values = "}", after=46)
  write.table(linhas.tabela, file=nome.arquivo, quote=FALSE, row.names=FALSE)
}
for (horizonte in c(1,3,6,12)) {
  nome.arquivo <- paste0("expand-hor",horizonte,".tex")
  legenda = paste("Janela expansível - horizonte", horizonte)
  tabela.expanding <- TabelaEQM(expanding,1:17, horizonte = horizonte,normaliza=TRUE)
  tmp <- apply(tabela.expanding,c(1,2),as.character)
  dados <- apply(tmp, c(1,2), substr, start=1, stop=4)
  #tabela.gw <- GiacominiWhite(moving,benchmark="rw",tamanho=120,horizonte=horizonte,vetor.maturidade)
  tabela.dm <- DieboldMariano(expanding,benchmark="rw",horizonte=horizonte,vetor.maturidade)
  signif <- t(apply(tabela.dm,c(1,2), Significancia))
  nomes.metodos <- row.names(dados)
  require(stringr)
  compilado <- str_c(dados,signif)
  tabela.final <- matrix(compilado, ncol=17)
  row.names(tabela.final) <- nomes.metodos
  colnames(tabela.final)<- tempo.maturidades
  #   View(tabela.final)
  #   View(dados)
  #   View(signif)
  tab.latex <- xtable(tabela.final, caption=legenda)
  print(tab.latex,type = "latex",file=nome.arquivo)
  linhas.tabela <- readLines(con=nome.arquivo)
  linhas.tabela <- append(linhas.tabela, values = "\\scalebox{0.7}{", after=4)
  linhas.tabela <- append(linhas.tabela, values = "}", after=46)
  write.table(linhas.tabela, file=nome.arquivo, quote=FALSE, row.names=FALSE)
}

Significancia <- function(p)  {
 if (is.nan(p)) return("")
 if ((p < 0) || (p>1)) stop("Número do p-valor inválido!")
 if (p < 0.01) return("***") 
 if (p < 0.05) return("**")
 if (p < 0.1) return("*")
 return("")  
}

# Significancia <- function(p)  {
#   if (is.nan(p)) return("")
#   if ((p < 0) || (p>1)) {
#     stop("Número do p-valor inválido!")
#   } 
#   if (p < 0.01) { 
#     return("***") 
#   } else if (p < 0.05) {
#     return("**")
#   } else if (p < 0.1) {
#     return("*")
#   } else {
#     return("")  
#   }
# }

tabela.gw <- GiacominiWhite(moving,benchmark="rw",tamanho=120,horizonte=1,vetor.maturidade)
require(plyr)
(tabela.sign <- aaply(tabela.gw,c(1,2), Significancia))


View(tabela.expanding)
View(tabela.gw)

tabela.dm <- DieboldMariano(moving,benchmark="rw",horizonte=1,vetor.maturidade)
View(tabela.dm)

x <- as.character(1.293829831)
y <- substr(x,start=1, stop=4)
str_c(y,"**")
xtable()

# Gráficos de retirar curvas ----------------------------------------------

pdf(width=6,height=4,file="perfil_curvas_retira0.pdf")
ajeita.janela()
curvas  <- PreparaCurvasCorte(taxas.juro,maturidade,intervalo.passado=c(1,300),intervalo.futuro=c(312,312))
tmp2 <- fts(x=tempo.maturidades, y = t(curvas$passado.learn), xname = "maturidade", yname = "taxas")
plot(tmp2)
dev.off()

pdf(width=6,height=4,file="perfil_curvas_retira7.pdf")
ajeita.janela()
curvas  <- PreparaCurvasCorte(taxas.juro,maturidade,intervalo.passado=c(1,300),intervalo.futuro=c(312,312), retirar = taxas.juro[, 7])
tmp2 <- fts(x=tempo.maturidades, y = t(curvas$passado.learn), xname = "maturidade", yname = "taxas")
plot(tmp2)
dev.off()

pdf(width=6,height=4,file="perfil_curvas_retira17.pdf")
ajeita.janela()
curvas  <- PreparaCurvasCorte(taxas.juro,maturidade,intervalo.passado=c(1,300),intervalo.futuro=c(312,312), retirar = taxas.juro[, 17])
tmp2 <- fts(x=tempo.maturidades, y = t(curvas$passado.learn), xname = "maturidade", yname = "taxas")
plot(tmp2)
dev.off()

pdf(width=6,height=4,file="perfil_curvas_retira1.pdf")
ajeita.janela()
curvas  <- PreparaCurvasCorte(taxas.juro,maturidade,intervalo.passado=c(1,300),intervalo.futuro=c(312,312), retirar = taxas.juro[, 1])
tmp2 <- fts(x=tempo.maturidades, y = t(curvas$passado.learn), xname = "maturidade", yname = "taxas")
plot(tmp2)
dev.off()


# Teste estacionariedade --------------------------------------------------

require("tseries")
for (i in 1:17) {
  serie <- taxas.juro[, i]
  for (j in c(1,7,17)) {
    if (i == j) next
    serie.ret <- serie - taxas.juro[, j]
    cat(i,",",j,",",sep="")
    cat(adf.test(x=ts(serie.ret),5,alternative="s")$p.value)
    cat("\n")
  }
}



# Heatmap -----------------------------------------------------------------

#tabela.moving <- TabelaEQM(moving,vetor.maturidade, horizonte = horizonte,normaliza=TRUE,subconjunto=intervalo)
exibir <- c("rw",
"ar",
"crt.p5.rmat6",
"crt.d1.rmat6",
"crt.p5.rb1",
"crt.d1.rb1",
"crt.p5.rmat17",
"crt.p1.rpmat",
"crt.d1.rmat17",
"crt.p5.rpmat",
"crt.d1.rpmat",
"diebold.05",
"ftsa.p3_ets.rmat17",
"ftsa.p5_rw.rmat7",
"ftsa.p5_ets.rmat1")


for (horizonte in c(1,3,6,12)) {
  nome.arquivo <- paste0("heatmap",horizonte,".pdf")
  pdf(width=6,height=4,file=nome.arquivo)
  par(mar=c(0.1, 0.1, 0.1, 0.1))
  tabela.expanding <- TabelaEQM(expanding,vetor.maturidade, horizonte = horizonte,normaliza=FALSE,subconjunto=intervalo, exibir=exibir)
  heatmap(tabela.expanding[-2, ],scale="column", Colv=NA, Rowv=NA,main=paste0("Horizonte ",horizonte))
  dev.off()
}


# Gráficos com erro cumulativo --------------------------------------------

#tabela.moving <- TabelaEQM(moving,vetor.maturidade, horizonte = horizonte,normaliza=TRUE,subconjunto=intervalo)
exibir <- c("rw",
            "ar",
            "crt.p5.rb1",
            "crt.d1.rb1",
            "crt.d1.rmat17",
            "crt.p5.rpmat",
            "crt.d1.rpmat",
            "diebold.05",
            "ftsa.p3_ets.rmat17",
            "ftsa.p5_rw.rmat7",
            "ftsa.p5_ets.rmat1")

# GRÁFICO DA
horizonte <- 1; maturidade <- 1 ####
nome.arquivo <- paste0("EQAP",horizonte,"-",maturidade,".pdf")
pdf(width=6,height=4,file=nome.arquivo)
par(mar=c(2.1, 2.1, 1.1, 1.1))
PlotarSeriesRw(expanding,maturidade=maturidade,horizonte=horizonte, exibir=exibir)
dev.off()

horizonte <- 3; maturidade <- 5 ####
nome.arquivo <- paste0("EQAP",horizonte,"-",maturidade,".pdf")
pdf(width=6,height=4,file=nome.arquivo)
par(mar=c(2.1, 2.1, 1.1, 1.1))
PlotarSeriesRw(expanding,maturidade=maturidade,horizonte=horizonte, exibir=exibir)
dev.off()

horizonte <- 6; maturidade <- 10 ####
nome.arquivo <- paste0("EQAP",horizonte,"-",maturidade,".pdf")
pdf(width=6,height=4,file=nome.arquivo)
par(mar=c(2.1, 2.1, 1.1, 1.1))
PlotarSeriesRw(expanding,maturidade=maturidade,horizonte=horizonte, exibir=exibir)
dev.off()

exibir <- c("rw",
            "ar",
            "crt.p5.rb1",
            "crt.d1.rb1",
            "crt.d1.rmat17",
            "crt.d1.rpmat",
            "diebold.05",
            "ftsa.p3_ets.rmat17",
            "ftsa.p5_rw.rmat7",
            "ftsa.p5_ets.rmat1")

horizonte <- 12; maturidade <- 17 ####
nome.arquivo <- paste0("EQAP",horizonte,"-",maturidade,".pdf")
pdf(width=6,height=4,file=nome.arquivo)
par(mar=c(2.1, 2.1, 1.1, 1.1))
PlotarSeriesRw(expanding,maturidade=maturidade,horizonte=horizonte, exibir=exibir)
dev.off()



tabela.expanding <- TabelaEQM(expanding,vetor.maturidade, horizonte = horizonte,normaliza=FALSE,subconjunto=intervalo, exibir=exibir)
