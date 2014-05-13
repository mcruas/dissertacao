# Código principal --------------------------------------------------------



# betas.02.a <- DieboldLi.EstimaBetas(tempo.maturidades, taxas.juro + 3,
#                                   datas, 0.2)
# betas.05.a <- DieboldLi.EstimaBetas(tempo.maturidades, taxas.juro + 3,
#                                   datas, 0.5)
# betas.02 <- cbind(beta_1=betas.02.a[, 1] - 3,betas.02.a[, 2:3]); rm(betas.02.a)
# betas.05 <- cbind(beta_1=betas.05.a[, 1] - 3,betas.05.a[, 2:3]); rm(betas.05.a)

tamanho.previsao <- diff(limites.estimacao) + 1
for (horizonte in vetor.horizontes) {
  for (i in 1:length(vetor.maturidade)) {
     maturidade <- vetor.maturidade[i]
     cat("[",horizonte,",",maturidade,"]\n")
     for (metodo in metodos.estimar) {
       simulacoes[[as.character(horizonte)]][[as.character(maturidade)]][[metodo]] <- rep(0,tamanho.previsao)
     }
      ## Valores reais ##
      if ("valores.reais" %in% metodos.estimar) { 
        valores.reais <- taxas.juro[py.range(limites.estimacao)+tamanho.janela-1+horizonte, maturidade]
        simulacoes[[as.character(horizonte)]][[as.character(maturidade)]][["valores.reais"]] <- valores.reais
      } 
      ## Prevê RW sem drift ##
      if ("rw" %in% metodos.estimar) { 
      taxas.previstas.rw <- taxas.juro[py.range(limites.estimacao)+tamanho.janela-1, maturidade]
      simulacoes[[as.character(horizonte)]][[as.character(maturidade)]][["rw"]] <- taxas.previstas.rw
      }
      for (j in 1:tamanho.previsao) {
        cat("I")
        if (window == "moving") {
          intervalo.passado <- limites.estimacao[1] + j - 1 +  c(0,tamanho.janela - 1) 
        } else { #window == "expanding"
          intervalo.passado <- limites.estimacao[1] +  c(0,tamanho.janela - 2 + j) 
        }
        intervalo.futuro <- intervalo.passado[2] + c(horizonte,horizonte)
        janela <- taxas.juro[py.range(intervalo.passado), maturidade]
        ## Prevê DieboldLi ##
        if ("diebold.05" %in% metodos.estimar) {
        betas.previstos.05 <- DieboldLi.PreveBetas(betas.05, intervalo.passado,
                                                     intervalo.futuro)
        taxas.previstas.diebold.05 <- DieboldLi.BetasParaTaxas(
            betas.previstos.05, tempo.maturidades, 0.5)
        simulacoes[[as.character(horizonte)]][[as.character(maturidade)]][["diebold.05"]][j] <- taxas.previstas.diebold.05[, maturidade]
        }
        if ("diebold.02" %in% metodos.estimar) {
        betas.previstos.02 <- DieboldLi.PreveBetas(betas.02, intervalo.passado,
                                                     intervalo.futuro)
        taxas.previstas.diebold.02 <- DieboldLi.BetasParaTaxas(
            betas.previstos.02, tempo.maturidades, 0.2)
        simulacoes[[as.character(horizonte)]][[as.character(maturidade)]][["diebold.02"]][j] <- taxas.previstas.diebold.02[, maturidade]
        }
        ## Prevê Ar(1) ##
        if ("ar" %in% metodos.estimar) {
        taxas.previstas.ar <- predict(ar(janela,order.max=1,method="ols"),n.ahead=horizonte)$pred[horizonte]
        simulacoes[[as.character(horizonte)]][[as.character(maturidade)]][["ar"]][j] <- taxas.previstas.ar
        }
        ## Prevê Corte ##
        metodos.corte <- str_split_fixed(metodos.estimar, "[.]",3)[
          which(str_split_fixed(metodos.estimar, "[.]",3)[, 1] == "crt"), ,drop=FALSE]
        if (nrow(metodos.corte) == 0 )   next
        for (k in 1:nrow(metodos.corte)) { # laço para explorar as diferentes estimações com o método corte
          if (metodos.corte[k,3] == "") { # caso não haja nada para retirar
            retirar <- rep(0,nrow(taxas.juro))
          } else if (grepl(pattern="rmat",metodos.corte[k,3])) { # caso retire alguma maturidade
            mat <- as.integer(sub(pattern="rmat","",metodos.corte[k,3]))
            retirar <- taxas.juro[, mat]
          } else if (grepl(pattern="rb1",metodos.corte[k,3])) { # caso retire o beta 1
            retirar <- betas.05[, 1] 
          } else if (grepl(pattern="rpmat",metodos.corte[k,3])) { # caso retire o valor da própria maturidade a estimar
            retirar <- taxas.juro[, maturidade]
          } else if (grepl(pattern="rmean",metodos.corte[k,3])) { # caso retira a média das taxa de juros de todas as maturidades
            retirar <- apply(taxas.juro,1,mean)
          }
          curvas <- PreparaCurvasCorte(taxas.juro,maturidade,intervalo.passado, intervalo.futuro,retirar = retirar)
          if (grepl(pattern="^d",metodos.corte[k,2])) { # caso a semimétrica seja da derivada
            grau <- as.integer(sub(pattern="d","",metodos.corte[k,2]))
            semimetricas <- SemimetricasClasse(curvas,q=grau,tipo="deriv")
            taxas.previstas.corte <- predict(curvas,semimetricas)
          } else if (grepl(pattern="^p",metodos.corte[k,2])) { # caso a semimétrica seja pca
            grau <- as.integer(sub(pattern="p","",metodos.corte[k,2]))
            semimetricas <- SemimetricasClasse(curvas,q=grau,tipo="pca")
            taxas.previstas.corte <- predict(curvas,semimetricas)
          } else if (grepl(pattern="^k",metodos.corte[k,2])) { # espera algo da forma kAA_BB, em que AA é o numero do knn e BB p grau da pca
            grau <- as.integer(sub("k.*_","",metodos.corte[k,2]))
            knn <- sub("_.*$","",metodos.corte[k,2]); knn <- as.integer(sub("^k","",knn))
            semimetricas <- SemimetricasClasse(curvas,q=grau,tipo="pca")
            taxas.previstas.corte <- predict(curvas,semimetricas,n.vizinhos = knn, cv = "knn")
          }
          nome.metodo <- paste(metodos.corte[k,], sep = "", collapse = "."); nome.metodo <- sub("[.]$","",nome.metodo)
          simulacoes[[as.character(horizonte)]][[as.character(maturidade)]][[nome.metodo]][j] <- taxas.previstas.corte
        }
      }
    }
}

## Métodos de séries temporais funcionais
metodos.ftsa.inteiro <- metodos.estimar[grep(pattern="^ftsa",metodos.estimar)]
if (nrow(metodos.ftsa) != 0 ) { # testa se há algum método de ftsa para ser estimado. Se falso, pula todo o código.
  metodos.ftsa <- str_split_fixed(metodos.ftsa.inteiro, "[.]",3)
  for (k in 1:nrow(metodos.ftsa)) { # laço para explorar as diferentes estimações com o método ftsa
    print(metodos.ftsa.inteiro[k])
    for (j in 1:tamanho.previsao) { # laço para cada intervalo em que a estimação será feita. obs: o ftsa estima para todos horizontes e maturidades ao mesmo tempo
      if (window == "moving") { 
        intervalo.passado <- limites.estimacao[1] + j - 1 +  c(0,tamanho.janela - 1) 
      } else { #window == "expanding"
        intervalo.passado <- limites.estimacao[1] +  c(0,tamanho.janela - 2 + j) 
      }
      if (metodos.ftsa[k,3] == "") { # caso não haja nada para retirar
        retirar <- rep(0,nrow(taxas.juro))
      } else if (grepl(pattern="rmat",metodos.ftsa[k,3])) { # caso retire alguma maturidade
        mat <- as.integer(sub(pattern="rmat","",metodos.ftsa[k,3]))
        retirar <- taxas.juro[, mat]
      } else if (grepl(pattern="rb1",metodos.ftsa[k,3])) { # caso retire o beta 1
        retirar <- betas.05[, 1] 
      } else if (grepl(pattern="rmean",metodos.ftsa[k,3])) { # caso retire o valor da própria maturidade a estimar
          retirar <- apply(taxas.juro,1,mean)
     } 
      curvas.fts <- PreparaCurvasFts(taxas.juro,intervalo.passado,tempo.maturidades, retirar = retirar)
      if (grepl(pattern="^p",metodos.ftsa[k,2])) {
        metodo.prever.fatores <- sub("p.*_","",metodos.ftsa[k,2])
        tmp <- sub("_.*$","",metodos.ftsa[k,2]); grau.pca <- as.integer(sub("^p","",tmp))
        previsoes <- predict(curvas.fts,vetor.horizontes, ordem.pca = grau.pca, metodo.previsao = metodo.prever.fatores)
      }
      nome.metodo <- metodos.ftsa.inteiro[k]
#           simulacoes[[as.character(horizonte)]][[as.character(maturidade)]][[nome.metodo]][j] <- taxas.previstas.corte
      #previsoes <- predict(curvas.fts,vetor.horizontes)
      for (mat in 1:length(vetor.maturidade)) {
        for (hor in 1:length(vetor.horizontes)) {
          simulacoes[[as.character(vetor.horizontes[hor])]][[as.character(vetor.maturidade[mat])]][[nome.metodo]][j] <- previsoes[mat,hor, drop=TRUE]
        }
      }
    }
  }
}
