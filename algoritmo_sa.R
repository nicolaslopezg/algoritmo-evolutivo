#======================
#Bibliotecas
#======================

library(ggplot2)
library(phylotools)
library(phytools)
library(phangorn)
library(ape)
library(ecr)
library(data.table)
library(apex)

# #======================
# #Parámetros de entrada
# #======================


# #======================
# #Datos globales
# #======================

cantidad_secuencias = 1011
cantidad_genes = 6015
max_genes = 3
cantidad_maxima = 10
archivo_ref = "/Users/kimn-it/Desktop/Uni/Arbol_ref/REF.tree"
arbol_ref = read.tree(archivo_ref)
arbol_ref_midpoint = midpoint(arbol_ref)
dir = "/Users/kimn-it/Desktop/Uni/Secuencias_alineadas/"
archivos = list.files(dir)

# #===================================
# #Definición de funciones
# #===================================

simmulated_annealing=function(Ciclos_t,Ciclos_s,T,alpha,modificacion){
  #Función para obtener la matriz de solucion de una lista
  obtener_matriz = function(cadena_actual,interseccion,indices,matriz_actual){
    nueva_matriz = NA
    matrices=NA
    largo = length(cadena_actual)-(length(interseccion)+1)
    if(length(interseccion) == 0){
      matrices = vector(mode = "list", length = length(cadena_actual))
    }
    else{
      largo = length(cadena_actual)-length(interseccion)+1
      matrices = vector(mode = "list", length = largo)
    }
    
    #Si no hay genes en comun se crea desde 0
    if(length(interseccion) == 0){
      for(i in 1:length(cadena_actual)){
        if(i == 1){
          archivo_actual = paste(dir,cadena_actual[i],sep = "")
          sec_actual = read.phyDat(archivo_actual, format = "phylip", type = "DNA")
          #nueva_matriz = as.data.frame(sec_actual)
          matrices[[i]] = sec_actual
        }
        else{
          archivo_actual = paste(dir,cadena_actual[i],sep = "")
          sec_actual = read.phyDat(archivo_actual, format = "phylip", type = "DNA")
          #sec_frame_data = as.data.frame(sec_actual)
          #nueva_matriz = concatenar_secs(nueva_matriz,sec_actual)
          matrices[[i]] = sec_actual
        }
      }
      matriz_actual = concatenate(new("multiphyDat",matrices))
      #nueva_matriz = concatenar_secs(matrices)
      return(matriz_actual)
    }
    #En caso contrario se deben guardar los genes en la matriz
    else{
      #Agregar luego para otros indices (recordar
      cont_temp = 0
      for (i in 1:length(interseccion)){
        archivo_actual = paste(dir,interseccion[i],sep = "")
        n_filas=as.integer((strsplit(readLines(archivo_actual,n=1)," ")[[1]][2]))
        cont_temp = cont_temp + n_filas
      }
      matriz_actual = matriz_actual[-(cont_temp+1:n_filas),]
      matrices[[1]] = matriz_actual
      cadena_restante = setdiff(cadena_actual,interseccion)
      for(i in 1:length(cadena_restante)) {
        archivo_actual = paste(dir,cadena_restante[i],sep = "")
        sec_actual = read.phyDat(archivo_actual, format = "phylip", type = "DNA")
        matrices[[i+1]] = sec_actual
      }
    }
    matriz_actual = concatenate(new("multiphyDat",matrices))
    #matriz_actual = concatenar_secs(matrices)
    return(matriz_actual)
  }
  #Funcion para normalizar la cantidad
  norm = function(cantidad){
    return(cantidad/100)
  }
  #Funcion para devolver el valor antes de normalizar
  desnorm = function(num){
    num = 100-(num*100)
    return(num)
  } 
  #Función que entrega el puntaje de la solucion
  puntaje_sol = function(matriz_solucion){
    hamming=dist.hamming(matriz_solucion)
    arbol_msa_midpoint = midpoint(NJ(as.matrix(hamming)))
    return(RF.dist(arbol_msa_midpoint, arbol_ref_midpoint,normalize = T))
  }
  #Función que entrega solución random
  sol_random = function(max_genes, cantidad_genes, archivos) {
    archivos_disponibles = archivos
    matrices = vector(mode = "list", length = max_genes)
    matrices[[1]] = read.phyDat(paste(dir, archivos[1], sep = ""), format = "phylip", type = "DNA")
    archivos_disponibles = archivos_disponibles[-1]
    for (i in 2:max_genes) {
      posicion = sample(length(archivos_disponibles), 1)
      matrices[[i]] = read.phyDat(paste(dir, archivos_disponibles[posicion], sep = ""), format = "phylip", type = "DNA")
      archivos_disponibles = archivos_disponibles[-posicion]
    }
    nueva_matriz = concatenate(new("multiphyDat",matrices))
    solucion = list(matriz = nueva_matriz, archivos = archivos[-which(archivos_disponibles %in% archivos)])
    return(solucion)
  }
  #Función que modifica la solucion segun un porcentaje
  mod_solucion = function(solucion_actual, modificacion, largo_genes){
    numero_genes = largo_genes
    cambio_genes = round(numero_genes * modificacion)
    cambio_cadena = round(cambio_genes * modificacion)
    archivos_usados = c()
    
    cadena_actual = solucion_actual$archivos
    archivos_diff = setdiff(archivos, cadena_actual)
    
    # Obtener la diferencia entre archivos y cadena_actual solo una vez
    if (length(archivos_usados) == 0) {
      archivos_diff = setdiff(archivos, cadena_actual)
    } else {
      archivos_diff = setdiff(archivos, archivos_usados)
    }
    
    # Modificación cadena
    posicion_temporal = numero_genes
    for (i in 1:cambio_genes) {
      posicion = sample(1:length(archivos_diff), 1)
      cadena_actual[posicion_temporal] = archivos_diff[posicion]
      posicion_temporal = posicion_temporal - 1
      archivos_usados = c(archivos_usados, archivos_diff[posicion])
    }
    
    #Agregar o quitar nuevo elemento a la cadena
    for (j in 1:cambio_cadena){
      archivos_diff = setdiff(archivos,archivos_usados)
      posicion = sample(1:(length(archivos_diff)),1)
      agregar_quitar = round(runif(1,0,4))
      if(agregar_quitar <= 1){
        if(length(cadena_actual) >= 3){
          cadena_actual = cadena_actual[-1]
        }
      }
      else{
        cadena_actual = c(cadena_actual,archivos_diff[posicion])
      }
    }
    
    #Acá se regula el tamaño
    if(length(cadena_actual) > cantidad_maxima){
      cadena_actual = cadena_actual[1:cantidad_maxima]
    }
    
    # Usar match en lugar de which para encontrar los índices de los elementos que coinciden en dos vectores
    interseccion = intersect(cadena_actual, solucion_actual$archivos)
    indices = match(interseccion, cadena_actual)
    
    nueva_matriz = obtener_matriz(cadena_actual, interseccion, indices, solucion_actual$matriz)
    solucion = list(matriz = nueva_matriz, archivos = cadena_actual)
    return(solucion)
  }
  #Funcion para reemplazar una solucion en el conjunto de soluciones
  reemplazar_solucion = function(vector_solucion,soluciones,posicion){
    soluciones[[posicion]] = vector_solucion[[1]]
    return(soluciones)
  }
  
  #######   Simulated annealing   ########
  
  soluciones_encontradas = c()
  solucion_actual = NULL
  nueva_solucion = NULL
  soluciones_encontradas = c()
  cadenas_encontradas = c()
  cadenas_genes=c()
  
  #Solucion inicial
  solucion_actual = sol_random(max_genes,cantidad_genes,archivos)
  puntaje_actual = puntaje_sol(solucion_actual$matriz)
  largo_actual = length(solucion_actual$archivos)
  vector_actual = c(puntaje_actual,norm(largo_actual))

  for(i in 1:ciclos_temperatura){
    gc()
    #Variar maximo numero de genes
    for(j in 1:ciclos_solucion){
      gc()
      #cadena_anterior = solucion_actual$archivos
      #largo_genes = length(solucion_actual$archivos)
      solucion_actual = mod_solucion(solucion_actual,modificacion,vector_actual[2]*100)
      #puntaje_actual = puntaje_sol(solucion_actual$matriz)
      puntaje_nueva = puntaje_sol(solucion_actual$matriz)
      #vector_actual = c(puntaje_actual,norm(length(solucion_actual$archivos)))
      vector_nueva = c(puntaje_nueva,norm(length(solucion_actual$archivos)))
      #Se consulta si la solucion nueva domina a la actual
      if(dominates(vector_nueva,vector_actual)){
        #solucion_actual = nueva_solucion
        vector_actual = vector_nueva
        #cadenas_genes = c(cadenas_genes,solucion_actual$archivos)
        #puntaje_actual = puntaje_nueva
        #Guarda info
        soluciones_encontradas = append(soluciones_encontradas,vector_nueva[1])
        #largo = length(solucion_actual$archivos)
        cadenas_encontradas = append(cadenas_encontradas,vector_nueva[2])
      }
      #Se consulta si la solucion actual domina a la nueva
      else if(dominates(vector_actual,vector_nueva)){
        delta=desnorm(vector_actual[1])-desnorm(vector_nueva[1])
        if((exp(-delta/T)) >  runif(1)){
          vector_actual = vector_nueva
          #cadenas_genes = c(cadenas_genes,cadenas_genes[length(cadenas_genes)])
          #solucion_actual = nueva_solucion
          #puntaje_actual = puntaje_nueva
          #Guarda info
          soluciones_encontradas = append(soluciones_encontradas,vector_nueva[1])
          #largo = length(solucion_actual$archivos)
          cadenas_encontradas = append(cadenas_encontradas,vector_nueva[2])
        }
        else{
          #Guarda info
          #cadenas_genes = c(cadenas_genes,cadena_anterior)
          soluciones_encontradas = append(soluciones_encontradas,vector_actual[1])
          #largo = length(solucion_actual$archivos)
          cadenas_encontradas = append(cadenas_encontradas,vector_actual[2])
        }
      }
      #Entonces ninguna domina si se entra en este punto
      else{
        vector_hv = rbind(vector_actual,vector_nueva)
        hipervolumen = computeHVContr(t(vector_hv), ref.point = c(2,2))
        #Si nueva solucion contribuye mas que la actual, se toma
        if(hipervolumen[2] > hipervolumen[1]){
          #solucion_actual = nueva_solucion
          vector_actual = vector_nueva
          #cadenas_genes = c(cadenas_genes,solucion_actual$archivos)
          #puntaje_actual = puntaje_nueva
          #Guarda info
          soluciones_encontradas = append(soluciones_encontradas,vector_nueva[1])
          #largo = length(solucion_actual$archivos)
          cadenas_encontradas = append(cadenas_encontradas,vector_nueva[2])
        }
        else{
          delta= abs(desnorm(hipervolumen[1]) - desnorm(hipervolumen[2]))
          if((exp(-delta/T)) >  runif(1)){
            vector_actual = vector_nueva
            #cadenas_genes = c(cadenas_genes,solucion_actual$archivos)
            #solucion_actual = nueva_solucion
            #puntaje_actual = puntaje_nueva
            #Guarda info
            soluciones_encontradas = append(soluciones_encontradas,vector_nueva[1])
            #largo = length(solucion_actual$archivos)
            cadenas_encontradas = append(cadenas_encontradas,vector_nueva[2])
          }
          else{
            #Guarda info
            #cadenas_genes = c(cadenas_genes,cadenas_genes[length(cadenas_genes)])
            soluciones_encontradas = append(soluciones_encontradas,vector_actual[1])
            largo = length(solucion_actual$archivos)
            cadenas_encontradas = append(cadenas_encontradas,vector_actual[2])
          }
        }
      }
    }
    T = T * alpha
  }
  #endTime = Sys.time()
  #tiempo = endTime - startTime
  #todas_soluciones_finales = append(todas_soluciones_finales,soluciones_encontradas[length(soluciones_encontradas)])
  #todas_cadenas_finales = append(todas_cadenas_finales,cadenas_encontradas[length(cadenas_encontradas)])
  
  #}
  
  #soluciones_sa = rbind(todas_soluciones_finales,todas_cadenas_finales)
  resultados =c()
  ciclos=c(1:length(cadenas_encontradas))
  
  for (i in 1:length(cadenas_encontradas)) {
    vector_prueba=rbind(c(soluciones_encontradas[i],cadenas_encontradas[i]))
    resultado=computeHV(t(vector_prueba), ref.point = c(2,2))
    resultados=append(resultados,resultado)
  }
  
  return(resultados)
  
  #Gráfica soluciones
  #plot(ciclos,resultados,type="l",col="black",pch=16,xlab="Iteraciones",ylab="Hipervolumen",main="Gráfica")
  
  ####### FIN   ########
}


