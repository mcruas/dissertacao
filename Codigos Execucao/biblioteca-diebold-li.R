# Um pequeno tutorial de como utilizar as funções pode ser lido em:
# http://rpubs.com/mcruas/diebold-li

library("termstrc")

DieboldLi.EstimaBetas <- function (maturidades, taxas.juro, datas, lambda) {
  #################################################################
  # Esta função estima os valores de beta para uma dada maturidade, taxa de
  # juros,datas e lambda.
  # Precisa da biblioteca termstrc!
  # ARGS
  #   betas: matriz com a série temporal dos betas nas colunas
  #   maturidades: vetor com a duração das maturidades
  #   lambda: o valor de lambda utilizado na estimação
  # RETORNA
  #   uma matriz com as os coeficientes beta dos fatores em função do tempo
  ################################################################  
  
  # cria um objeto zeroyields e depois estima-se 
  datazeroyields1 <- zeroyields(maturidades, taxas.juro, datas)
  dl_res <- estim_nss(dataset=datazeroyields1, method="dl", lambda=lambda)
  output <- as.data.frame(param(dl_res)[[1]])
  names(output) <- c("beta_1","beta_2","beta_3")
  return(as.matrix(output))
  
}


DieboldLi.PreveBetas <- function(betas, intervalo.passado, intervalo.futuro, metodo = "ar") {
  ############################################################
  # Utiliza 
  # ARGS
  #   betas: matriz com a série temporal dos betas nas colunas
  #   intervalo.passado: vetor c(x,y), sendo x o início e y o fim
  #                     do intervalo do treino
  #   intervalo.futuro: vetor c(x,y), sendo x o início e y o fim
  #                     do intervalo a ser previsto
  #   metodo: o método para modelar e prever
  #               "ar" (default) para usar um AR(1) e "arima" para usar a função
  #               auto.arima, do pacote forecast
  # RETORNA
  #   uma matriz com as taxas de juros em função do tempo e 
  #   das maturidades
  ###########################################################
  
  # estima previsões para cada beta
  betas.previstos <- matrix(data=0,nrow=(intervalo.futuro[2]-
                                           intervalo.futuro[1]+1),ncol=3)
  horizonte  <- intervalo.futuro[1] - intervalo.passado[2]
  for (i in intervalo.futuro[1]:intervalo.futuro[2]) {
    # betas.arima é uma lista com o arima aplicado para cada um dos betas.
    # O beta de h passos à frente é obtido através de um AR(1) ou ARIMA
    if (metodo == "arima") {
      betas.arima <- apply(betas[intervalo.passado[1]:(i - horizonte), ],
                           2,auto.arima)
    } else {
      betas.arima.1 <- ar(betas[intervalo.passado[1]:(i - horizonte), 1],
                          order.max = 1, method = "ols")            
      betas.arima.2 <- ar(betas[intervalo.passado[1]:(i - horizonte), 2],
                          order.max = 1, method = "ols")            
      betas.arima.3 <- ar(betas[intervalo.passado[1]:(i - horizonte), 3],
                          order.max = 1, method = "ols")            
    }                   
    betas.previstos[i - intervalo.futuro[1] + 1, ] <- 
      c(predict(betas.arima.1, n.ahead=horizonte)$pred[horizonte],
        predict(betas.arima.2, n.ahead=horizonte)$pred[horizonte],
        predict(betas.arima.3, n.ahead=horizonte)$pred[horizonte])
  }
  return(betas.previstos)
}


DieboldLi.BetasRW <- function(betas, intervalo.passado, intervalo.futuro, horizonte) {
  ############################################################
  # Utiliza 
  # ARGS
  #   betas: matriz com a série temporal dos betas nas colunas
  #   intervalo.passado: vetor c(x,y), sendo x o início e y o fim
  #                     do intervalo do treino
  #   intervalo.futuro: vetor c(x,y), sendo x o início e y o fim
  #                     do intervalo a ser previsto
  #   horizonte: o horizonte de previsão
  # RETORNA
  #   uma matriz com as taxas de juros em função do tempo e 
  #   das maturidades
  ###########################################################
  betas.previstos <- betas[py.range(intervalo.futuro) - horizonte, ]
  return(t(betas.previstos))
}


DieboldLi.BetasParaTaxas <- function(betas,tempo.maturidades,lambda) {
  ############################################################
  # Transforma a série temporal dos coeficientes dos fatores em
  # série temporal das taxas de juros para as maturidades desejadas
  # ARGS
  #   betas: matriz com a série temporal dos betas nas colunas
  #   maturidades: vetor com a duração das maturidades
  #   lambda: o valor de lambda utilizado na estimação
  # RETORNA
  #   uma matriz com as taxas de juros em função do tempo e 
  #   das maturidades
  ###########################################################
  output <-  t(apply(betas,1,spr_dl,tempo.maturidades,lambda))
  colnames(output) <- tempo.maturidades * 12
  return(output)
}

py.range <- function(range){
  # implements python function range, but starting on 1 instead of 0
  return(range[1]:range[2])
}
