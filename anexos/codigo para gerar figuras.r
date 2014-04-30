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
for (horizonte in c(1,3,6,12)) {
  nome.arquivo <- paste0("movel-hor",horizonte,".tex")
  legenda = paste("Janela móvel - horizonte", horizonte)
  tabela.moving <- TabelaEQM(moving,1:17, horizonte = horizonte,normaliza=TRUE)
  tab.latex <- xtable(tabela.moving, caption=legenda)
  print(tab.latex,type = "latex",file=nome.arquivo)
}
for (horizonte in c(1,3,6,12)) {
  nome.arquivo <- paste0("expand-hor",horizonte,".tex")
  legenda = paste("Janela expandindo - horizonte", horizonte)
  tabela <- TabelaEQM(expanding,1:17, horizonte = horizonte,normaliza=TRUE)
  tab.latex <- xtable(tabela, caption=legenda)
  print(tab.latex,type = "latex",file=nome.arquivo)
}
for (horizonte in c(1,3,6,12)) {
  nome.arquivo <- paste0("gw.movel-hor",horizonte,".tex")
  legenda = paste("Teste de Giacomini-White para janela móvel - horizonte", horizonte)
  tabela.moving <- TabelaEQM(moving,1:17, horizonte = horizonte,normaliza=TRUE)
  tab.latex <- xtable(tabela.moving, caption=legenda)
  print(tab.latex,type = "latex",file=nome.arquivo)
}

tabela.gw <- GiacominiWhite(moving,benchmark="rw",tamanho=120,horizonte=1,vetor.maturidade)
View(tabela.gw)

tabela.dm <- DieboldMariano(moving,benchmark="rw",horizonte=1,vetor.maturidade)
View(tabela.dm)



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

for (horizonte in c(1,3,6,12)) {
  nome.arquivo <- paste0("heatmap",horizonte,".pdf")
  pdf(width=6,height=4,file=nome.arquivo)
  par(mar=c(1.1, 1.1, 0.5, 1.1))
  tabela.expanding <- TabelaEQM(expanding,vetor.maturidade, horizonte = horizonte,normaliza=TRUE,subconjunto=intervalo)
  heatmap(tabela.expanding[-2, ],scale="column", Colv=NA, Rowv=NA)
  dev.off()
}

heatmap(tabela.expanding,scale="column", Colv=NA, Rowv=NA)

