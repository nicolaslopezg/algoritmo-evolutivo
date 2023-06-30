#======================
#Bibliotecas
#======================
library("ggplot2")
library("irace")
source("/Users/kimn-it/Desktop/Uni/Codigo/algoritmo_sa.R")

#======================
#Configuración de irace
#======================

#=====================
# Definición del runner
#=====================

target.runner = function(experiment, scenario){
  
  set.seed(1)
  
  #Ajuste de entrada
  #entrada=experiment$instance
  entrada=strsplit(entrada,"/")
  entrada=entrada[[1]][length(entrada[[1]])]
  
  #Otros parámetros
  Ciclos_t=as.numeric(experiment$configuration[["Ciclos_t"]])
  Ciclos_s=as.numeric(experiment$configuration[["Ciclos_s"]])
  T=as.numeric(experiment$configuration[["T"]])
  alpha=as.numeric(experiment$configuration[["alpha"]])
  modificacion=as.numeric(experiment$configuration[["modificacion"]])
  
  
  resultado=simmulated_annealing(Ciclos_t,Ciclos_s,T,alpha,modificacion)
  
  return(list(cost =resultado))
}

# Lectura de scenario
escenario = readScenario(filename = "/Users/kimn-it/Desktop/Uni/Codigo/Tuning/escenario.txt", scenario = defaultScenario())

# Lectura de parámetros
parametros = readParameters(file =  "/Users/kimn-it/Desktop/Uni/Codigo/Tuning/parametros.txt")

escenario$targetRunner=target.runner

irace(scenario = escenario, parameters = parametros)
