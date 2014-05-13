###########################################################
# Conjunto de funções 
#
#
###########################################################


AnimacaoJuros <- function(taxas.juro, tempo.maturidades, file) {
  
  image <- array(0, c(2,17, dim(taxas.juro)[1]))
  for (i in 1:(dim(taxas.juro)[1])) {
    image[,,i] <- as.matrix(cbind(tempo.maturidades,taxas.juro[i,]))
  }
  View(image[,,i])
}

Ar.Previsao <- function(ts, intervalo.passado, intervalo.futuro) {
  #########################################################
  # Realiza previsões para o intervalo futuro. O horizonte é definido pela
  # diferença entre os intervalos e o modelo é um AR(1) por default
  # ARGS
  #   ts: a série temporal para prever
  #   intervalo.passado: vetor de tamanho 2 com o início e fim do intervalo de treino
  #   intervalo.futuro: vetor de tamanho 2 com o início e fim do intervalo para previsão
  #   grau.arima: vetor de ordem 3 que especifica os graus do modelo arima c(p,i,q)
  # RETORNA
  # Um vetor com as previsões
  #########################################################
  previstos <- rep(0,times = diff(intervalo.futuro) + 1)
  horizonte  <- intervalo.futuro[1] - intervalo.passado[2]
  for (i in intervalo.futuro[1]:intervalo.futuro[2]) {
    estimacao.arima <- ar(ts[intervalo.passado[1]:(i - horizonte)],order.max=1,method="ols")
    previstos[i - intervalo.futuro[1] + 1] <- predict(estimacao.arima,
                                                      n.ahead=horizonte, aic = FALSE)$pred[horizonte]
  }
  return(previstos)
}



CombinarPrevisoes <- function(df, ...) {
  ####################
  df <- as.data.frame(df)
  metodos <- c(...)
  tmp <- apply(sapply(metodos, function (x) df[, which(names(df)==x)]),1,sum)/
    length(metodos)
  df <- cbind(tmp, df)
  names(df)[1] <- paste(metodos,collapse=' + ')
  return(as.data.frame(as.matrix(df)))
}


cummean <- function(x) {
  ############################################################ 
  #Calculates the cumulative mean of a vector. Works like the function cumsum, but
  #taking its mean
  ###########################################################
  if (is.vector(x)) output <- cumsum(x)/1:length(x)
  return(output)
}

EQM <- function(d1,d2) {
  ##########################################################
  # Calcula o erro quadrático médio entre dois vetores
  #########################################################
  mean((d1-d2)^2)
}


ExtraiSeries <- function(arquivo,horizonte,maturidade) {
  #######################################################
  # Exibe uma tabela com o erro quadrático médio de cada método para cada conjunto
  # de maturidade x horizonte x intervalo.
  ######################################################
  series <- as.data.frame(arquivo[[as.character(horizonte)]][[as.character(maturidade)]])
  return(series)
}


GiacominiWhite <- function(arquivo, benchmark, tamanho, horizonte, vetor.maturidade) {
  ############################################################
  # Aplica o teste de Giacomini-White para uma sequência de previsões, comparando
  # todas com a capacidade de previsão de um mesmo método como benchmarking
  # ARGS
  #   tabela: tabela com os valores previstos para cada um dos métodos, incluindo
  #                 o benchmarking e os valores.reais. Séries nas colunas.
  #   benchmarking: o nome do método benchmark
  #   tamanho: o tamanho da amostra que gerou as previsões
  #   horizonte: o horizonte das previsões
  #   vetor.maturidade: vetor com o índice das maturidades para testar
  # RETORNO
  #   um vetor com os valores da estatística de Giacomini-White, para cada  método
  ############################################################# 
  #install.packages(c("polynom", "fracdiff", "hypergeo", "longmemo")) 
  #install.packages("~/Dropbox/R/NPFDA/afmtools_0.1.8.tar.gz", repos = NULL, type
  #= "source")
  tabela.gw <- NULL
  for (maturidade in vetor.maturidade) {
    tabela <- as.data.frame(arquivo[[as.character(horizonte)]][[as.character(maturidade)]])
    i.valores.reais <- which(names(tabela) == "valores.reais")
    i.benchmark <- which(names(tabela) == benchmark)
    serie.benchmark <- tabela[, i.benchmark]; serie.real <- tabela[, i.valores.reais]
    tabela.tmp  <- tabela[, c(- i.benchmark,- i.valores.reais)] ; tabela.tmp <- tabela.tmp[, order(names(tabela.tmp))] # coloca as colunas em ordem alfabética
    if (!require(afmtools)) stop("Por favor, instale a biblioteca afmtools. Instrucoes na descricao da funcao.")
    tmp <- apply(tabela.tmp, 2, gw.test, y = serie.benchmark, p = serie.real,
                 T = tamanho, tau = horizonte, method="HAC", alternative = "two.sided")
    resultados <- unlist(lapply(tmp,function(x) x$p.value))
    #gw.test(x = tmp[,7], y = tmp[,3], p = tmp[,12], T = 800, tau = 5,method="HAC")
    tabela.gw <- cbind(tabela.gw,resultados)
  }
  row.names(tabela.gw) <- tempo.maturidades
  return(tabela.gw)
}

DieboldMariano <- function(arquivo, benchmark, horizonte, vetor.maturidade) {
  ############################################################
  # Aplica o teste de Diebold-Mariano para uma sequência de previsões, comparando
  # todas com a capacidade de previsão de um mesmo método como benchmarking
  # ARGS
  #   arquivo: lista contendo as simulações
  #   benchmarking: o nome do método benchmark
  #   horizonte: o horizonte das previsões
  #   vetor.maturidade: vetor com o índice das maturidades para testar
  # RETORNO
  #   uma matriz com os valores da estatística de Giacomini-White, para cada método e maturidade
  ############################################################# 
  tabela.est <- NULL
  if (!require(forecast)) stop("Por favor, instale a biblioteca forecast.")
  for (maturidade in vetor.maturidade) {
    series <- as.data.frame(arquivo[[as.character(horizonte)]][[as.character(maturidade)]])
    series <- series[, order(names(series))] # coloca as colunas em ordem alfabética
    i.valores.reais <- which(names(series) == "valores.reais")
    i.benchmark <- which(names(series) == benchmark)
    series.erros <- as.data.frame((series[, c(-i.valores.reais)] - series[, i.valores.reais])^2)    
    tmp <- apply(series.erros[, -i.benchmark], 2, dm.test, 
                 e2 = series.erros[, i.benchmark], alternative = "less",
                 h = horizonte, power=2)
    resultados <- unlist(lapply(tmp,function(x) x$p.value))
    #gw.test(x = tmp[,7], y = tmp[,3], p = tmp[,12], T = 800, tau = 5,method="HAC")
    tabela.est <- rbind(tabela.est,resultados)
  }
  row.names(tabela.est) <- tempo.maturidades
  return(tabela.est)
}



IntervaloFuturo <- function(intervalo,percentual.testar, horizonte) {
  #############################################################
  # ARGS
  #   intervalo: vetor tipo c(x,y)
  #   percentual.testar: um número entre 0 e 1
  # RETORNO
  #   um vetor c(z,y), em que z é o primeiro valor que será previsto  
  #   Ex: IntervaloFuturo(c(1001,2000),0.8,h) retorna c(1801 + h,2000)
  #############################################################
  
  tamanho.serie.temporal  = intervalo[2]-intervalo[1]+1
  output <- c(trunc((1-percentual.testar)*tamanho.serie.temporal
  ) + intervalo[1]-1 + horizonte, tamanho.serie.temporal + 
    intervalo[1] - 1)
  return(output)
}

IntervaloPassado <- function(intervalo,percentual.testar) {
  #############################################################
  # ARGS
  #   intervalo: vetor tipo c(x,y)
  #   percentual.testar: um número entre 0 e 1
  # RETORNO
  #   um vetor c(x,z), em que z corresponde ao (1-percentual.testar)
  #   do intervalo inicial. 
  #   Ex: Intervalo.base.previsao(c(1001,2000),0.8) retorna c(1001,1800)
  #############################################################
  output <- intervalo
  output[2] <- trunc((intervalo[2]-intervalo[1]+1)*
                       (1-percentual.testar) + intervalo[1] - 1)
  return(output)
}




Multiplot.ts <- function(data, ... , range, dates){
  ############################################################
  # This function plots many time series in the same graph
  #
  # ARGS  
  #    Data is a matrix with each serie in its lines
  #    Range is an optional vector with the indexes of the time series
  #       to be exibited. If it is absent, then every line is ploted
  ############################################################
  if (missing(range)) 
    range <- 1:dim(data)[2]
  size.series <- dim(data)[2]
  if (missing(dates)) {
    plot(data[, range[1]],col=1,ylim=c(min(data[, range]), 
                                       max(data[, range])), type="l",  ...)    
    for (i in range[-1]) {
      lines(data[, i],col=i, lty=((i - 1) %/% 8 + 1), ...)   
    }
  } else {
    plot(dates, data[, range[1]],col=1,ylim=c(min(data[, range]), 
                                              max(data[, range])), type="l",  ...)
    for (i in range[-1]) {
      lines(dates,data[, i],col=i, lty=((i - 1) %/% 8 + 1), ...)   
    }
  } 
}  

Nome.vetor <- function(object) {
  return(paste0(object[1],":",object[2]))
}

Plot.list <- function(object) {
  tmp <- as.data.frame(object)
  Multiplot.ts(tmp)
}

PlotarSeries <- function(arquivo,maturidade,horizonte,subconjunto = NULL, exibir = NULL, graficos.erros = FALSE) {
  ####################################################### 
  # Exibe 3 gráficos: I) séries previstas por cada método e a realizada; II) erro
  # quadrático de cada série prevista, por unidade de tempo; III) erro quadrático
  # médio de cada série, até a unidade de tempo t. 
  ######################################################
  series <- as.data.frame(arquivo[[as.character(horizonte)]][[as.character(maturidade)]])
  if (!is.null(subconjunto)) series <- series[py.range(subconjunto), ]
  series <- series[, order(names(series))] # coloca as colunas em ordem alfabética
  if (!is.null(exibir)) series <- series[, (which(names(series) %in% c(exibir,"rw","valores.reais")))]
  #series <- Combinar.Previsoes(series,"corte.pca.5","diebold.0.5") # para combinar previsões
  series.erros <- as.data.frame((series[, -which(names(series)=="valores.reais")] - series[, which(names(series)=="valores.reais")])^2)
  series.erros.medios <- as.data.frame(apply(series.erros,2,FUN=cummean))
  #eqm <- apply(series[, -which(names(series)=="valores.reais")],2,FUN=EQM,series[, which(names(series)=="valores.reais")])
  parte.titulo <- paste("horizonte ", as.character(horizonte), " | maturidade", as.character(tempo.maturidades[maturidade]))
  # par(xpd=T, mar=par()$mar+c(0,0,0,4)) # Diminui o espaço do gráfico à direita para ter mais espaço para a legenda
  if (graficos.erros == TRUE) { # exibe as séries de erros
    Multiplot.ts(series.erros.medios); title(paste("Series dos erros médios: ",parte.titulo),cex.main=0.8)
    legend("bottomleft", legend=names(series.erros.medios), lty = (1:length(names(series.erros.medios)) - 1) %/% 8 + 1, col=1:length(names(series.erros.medios)), title="Métodos",cex=0.6, y.intersp=0.5)
    Multiplot.ts(series.erros); title(paste("Series dos erros: ",parte.titulo),cex.main=0.8)
    legend("bottomleft", legend=names(series.erros), lty = (1:length(names(series)) - 1) %/% 8 + 1, col=1:length(names(series.erros)), title="Métodos",cex=0.6, y.intersp=0.5)
  }
  Multiplot.ts(series,ylab="taxas de juro"); title(paste("Series: ",parte.titulo),cex.main=0.8)
  legend("bottomleft", legend=names(series), lty = (1:length(names(series)) - 1) %/% 8 + 1, col=1:length(names(series)), title="Métodos",cex=0.6, y.intersp=0.5)
  # par(mar=c(5, 4, 4, 2) + 0.1) # retorna 'par' para seu valor original
  # print(eqm)
}


PlotarSeriesRw <- function(arquivo,maturidade,horizonte,subconjunto = NULL, exibir = NULL) {
  ####################################################### 
  # Exibe 3 gráficos: I) séries previstas por cada método e a realizada; II) erro
  # quadrático de cada série prevista, por unidade de tempo; III) erro quadrático
  # médio de cada série, até a unidade de tempo t. 
  ######################################################
  series <- as.data.frame(arquivo[[as.character(horizonte)]][[as.character(maturidade)]])
  if (!is.null(subconjunto)) series <- series[py.range(subconjunto), ]
  series <- series[, order(names(series))] # coloca as colunas em ordem alfabética
  series.erros <- as.data.frame((series[, -which(names(series)=="valores.reais")] - series[, which(names(series)=="valores.reais")])^2)
  series.erros.cum <- as.data.frame(apply(series.erros,2,FUN=cumsum))
  series.erros.rw <- series.erros.cum[, which(names(series.erros.cum)=="rw")] - series.erros.cum[, -which(names(series.erros.cum)=="rw")]
  if (!is.null(exibir)) series.erros.rw <- series.erros.rw[, (which(names(series.erros.rw) %in% exibir))]
  #series.erros.rw <- series.erros.medios[, which(names(series.erros.medios)=="rw")] - series.erros.medios[, -which(names(series.erros.medios)=="rw")]
  parte.titulo <- paste("horizonte ", as.character(horizonte), " | maturidade", as.character(tempo.maturidades[maturidade]))
#  par(xpd=T, mar=par()$mar+c(0,0,0,4)) # Diminui o espaço do gráfico à direita para ter mais espaço para a legenda
  Multiplot.ts(series.erros.rw); title(paste("Erros quadráticos acumulados de previsão: ",parte.titulo),cex.main=0.8)
  abline(a=0,b=0)
  #legend("topright", inset=c(-0.3,0), legend=names(series.erros.rw), lty = (1:length(names(series.erros.rw)) - 1) %/% 8 + 1, col=1:length(names(series.erros.rw)), title="Métodos",cex=0.6, y.intersp=0.5)
  legend("bottomleft", legend=names(series.erros.rw), lty = (1:length(names(series.erros.rw)) - 1) %/% 8 + 1, col=1:length(names(series.erros.rw)), title="Métodos",cex=0.8, y.intersp=0.5)
#  par(mar=c(5, 4, 4, 2) + 0.1) # retorna 'par' para seu valor original
}

Tabelao <- function(arquivo,vetor.maturidades,vetor.horizontes,vetor.intervalos, exibir = NULL) {
  #######################################################
  # Exibe uma tabela com o erro quadrático médio de cada método para cada conjunto
  # de maturidade x horizonte x intervalo.
  ######################################################
  tabela <- NULL
  for (i in vetor.maturidades)
    for (j in vetor.horizontes)
      for(k in vetor.intervalos) {
        series <- as.data.frame(arquivo[[Nome.vetor(k)]][[as.character(j)]][[as.character(i)]])
        series <- 
        eqm <- apply(series[, -which(names(series)=="valores.reais")],2,FUN=EQM,series[, which(names(series)=="valores.reais")])
        if (!is.null(exibir)) eqm <- eqm[, (which(names(eqm)) %in% exibir)]
        nome.linha <- paste(j,"|",tempo.maturidades[i],"|",k[1],"-",k[2],sep='')
        tabela <- rbind(tabela,eqm)
        rownames(tabela)[dim(tabela)[1]] <- nome.linha
      }
  return(tabela)
}

TabelaoAgrupado <- function(arquivo,vetor.maturidade,horizonte,vetor.intervalos) {
  #######################################################
  # Exibe uma tabela com o erro quadrático médio de cada método para cada conjunto
  # de maturidade x horizonte x intervalo.
  ######################################################
  tabela <- NULL
  for (i in vetor.maturidade) {
    tabela.tmp <- Tabelao(arquivo, i, horizonte, vetor.intervalos)
    tabela <- cbind(tabela,colMeans(tabela.tmp))
  }
  colnames(tabela) <- tempo.maturidades[vetor.maturidade]
  return(t(tabela))
}

TabelaEQM <- function(arquivo,vetor.maturidade, horizonte, normaliza = TRUE, subconjunto = NULL, exibir = NULL, gw.test = FALSE, ...) {
  #######################################################
  # Exibe uma tabela com o erro quadrático médio de cada método para cada
  # de maturidade x método, para um determinado horizonte.
  ######################################################
  tabela <- NULL
  for (maturidade in vetor.maturidade) { 
    # arquivo <- simulacoes
    series <- as.data.frame(arquivo[[as.character(horizonte)]][[as.character(maturidade)]])
    serie.real <- series[, which(names(series)=="valores.reais")]
    series <- series[, -which(names(series)=="valores.reais")]
    series <- series[, order(names(series))] # coloca as colunas em ordem alfabética
    if (!is.null(subconjunto)) series <- series[py.range(subconjunto), ]
    eqm <- apply(series,2,FUN=EQM,serie.real)
    tabela <- cbind(tabela,eqm)
  }
  colnames(tabela) <- tempo.maturidades[vetor.maturidade]
  if (normaliza) tabela <- t(t(tabela[-which(rownames(tabela)=="rw"), ]) / tabela[which(rownames(tabela)=="rw"), ])
  tabela <- sqrt(tabela)
  if (!is.null(exibir)) tabela <- tabela[(which(rownames(tabela) %in% exibir)), ]
  tabela <- tabela[order(rownames(tabela)), ] # coloca as colunas em ordem alfabética  
  return(tabela)
}
