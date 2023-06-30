##################################################### 
#  Algoritmo evolutivo para minimización de genes  #
#                  Nicolás López                   #
####################################################

library(ShortRead)
library(ggplot2)
library(phylotools)
library(phytools)
library(phangorn)
library(ape)
library(ecr)
library(data.table)

#Definición de variables
cantidad_secuencias = 1011
cantidad_genes = 6015
max_genes = 5
cantidad_maxima = 50
archivo_ref = "/Users/kimn-it/Desktop/Uni/Arbol_ref/REF.tree"
#archivo_ref = "/Users/Nicolás López/Desktop/TESIS/Arbol_ref/REF.tree"
arbol_ref = read.tree(archivo_ref)
arbol_ref_midpoint = midpoint(arbol_ref)

#########################################################
##################### FUNCIONES #########################
#########################################################

#Funcion concatenar
concatenar_secs = function(matrices){
  for (i in 1:length(matrices)) {
    #nueva_sec[[i]] = as.data.frame(matrices[[i]])
    matrices[[i]] = as.data.frame(matrices[[i]])
  }
  nueva_sec = do.call(rbind,matrices)
  #sec_1 = as.data.frame(sec_1)
  #sec_2 = as.data.frame(sec_2)
  #nueva_sec = rbind(sec_1,sec_2)
  return(nueva_sec)
}

#Función para obtener la matriz de solucion de una lista
#(EDITADA)
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
    nueva_matriz = concatenar_secs(matrices)
    return(nueva_matriz)
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
  matriz_actual = concatenar_secs(matrices)
  return(matriz_actual)
}

#Funcion para normalizar la cantidad
norm = function(cantidad){
  return(cantidad/100)
}

#Función que entrega el puntaje de la solucion
#(EDITADA)
puntaje_sol = function(matriz_solucion){
  matriz_phyDat = as.phyDat(matriz_solucion)
  hamming=dist.hamming(matriz_phyDat)
  arbol_msa = NJ(as.matrix(hamming))
  arbol_msa_midpoint = midpoint(arbol_msa)
  dif_arboles = RF.dist(arbol_msa_midpoint, arbol_ref_midpoint,normalize = T)
  #puntaje = 100 -(dif_arboles*100)
  return(dif_arboles)
}

#Función que entrega solución random
#(EDITADA)
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
  nueva_matriz = concatenar_secs(matrices)
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
      cadena_actual = cadena_actual[-1]
    }
    else{
      cadena_actual = c(cadena_actual,archivos_diff[posicion])
    }
  }
  
  #Acá se regula el tamaño
  if(length(cadena_actual) > cantidad_maxima){
    cadena_actual = cadena_actual[1:50]
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

#########################################################################
#########################################################################

#Dirección archivos
dir = "/Users/kimn-it/Desktop/Uni/Secuencias_alineadas/"
#dir =  "/Users/Nicolás López/Desktop/TESIS/Secuencias_alineadas/"

#Lectura todos los archivos
archivos = list.files(dir)

#Grafica
#plot(arbol_msa,cex=0.6,type = "phylogram",use.edge.length =T)

#Impresion arbol de referencia
#plot(arbol_ref_midpoint,cex=0.6,type = "phylogram",use.edge.length =T)


######################### SA ############################

ciclos_temperatura = 10 #Ciclos externos
ciclos_solucion = 4    #Ciclos para solucion en especifico
modificacion = 0.6    #Porcentaje modificacion de solucion
aceptacion = 0.95
T = 100

#todas_soluciones_finales = c()
#todas_cadenas_finales = c()

#for (ejecucion in 1:10) {
  
  #max_genes = as.integer(runif(1,20,50))
  
  #Reiniciando variables
  soluciones_encontradas = c()
  solucion_actual = NULL
  nueva_solucion = NULL
  soluciones_encontradas = c()
  cadenas_encontradas = c()
  
  #Luego descomentar esta linea
  #startTime = Sys.time()
  solucion_actual = sol_random(max_genes,cantidad_genes,archivos)
  puntaje_actual = puntaje_sol(solucion_actual$matriz)
  for(i in 1:ciclos_temperatura){
    #Variar maximo numero de genes
    for(j in 1:ciclos_solucion){
      largo_genes = length(solucion_actual$archivos)
      nueva_solucion = mod_solucion(solucion_actual,modificacion,largo_genes)
      #puntaje_actual = puntaje_sol(solucion_actual$matriz)
      puntaje_nueva = puntaje_sol(nueva_solucion$matriz)
      vector_actual = c(puntaje_actual,norm(length(solucion_actual$archivos)))
      vector_nueva = c(puntaje_nueva,norm(length(nueva_solucion$archivos)))
      #Se consulta si la solucion nueva domina a la actual
      if(dominates(vector_nueva,vector_actual)){
        solucion_actual = nueva_solucion
        puntaje_actual = puntaje_nueva
        #Guarda info
        soluciones_encontradas = append(soluciones_encontradas,puntaje_nueva)
        largo = length(solucion_actual$archivos)
        cadenas_encontradas = append(cadenas_encontradas,largo)
      }
      #Se consulta si la solucion nueva domina a la actual
      else if(dominates(vector_actual,vector_nueva)){
        delta=puntaje_actual - puntaje_nueva
        if((exp(-delta/T)) >  aceptacion){
          solucion_actual = nueva_solucion
          puntaje_actual = puntaje_nueva
          #Guarda info
          soluciones_encontradas = append(soluciones_encontradas,puntaje_nueva)
          largo = length(solucion_actual$archivos)
          cadenas_encontradas = append(cadenas_encontradas,largo)
        }
        else{
          #Guarda info
          soluciones_encontradas = append(soluciones_encontradas,puntaje_actual)
          largo = length(solucion_actual$archivos)
          cadenas_encontradas = append(cadenas_encontradas,largo)
        }
      }
      #Entonces ninguna domina si se entra en este punto
      else{
        vector_hv = rbind(vector_actual,vector_nueva)
        hipervolumen = computeHVContr(t(vector_hv), ref.point = c(2,2))
        #Si nueva solucion contribuye mas que la actual, se toma
        if(hipervolumen[2] > hipervolumen[1]){
          solucion_actual = nueva_solucion
          puntaje_actual = puntaje_nueva
          #Guarda info
          soluciones_encontradas = append(soluciones_encontradas,puntaje_nueva)
          largo = length(solucion_actual$archivos)
          cadenas_encontradas = append(cadenas_encontradas,largo)
        }
        else{
          delta= hipervolumen[1] - hipervolumen[2]
          if((exp(-delta/T)) >  aceptacion){
            solucion_actual = nueva_solucion
            puntaje_actual = puntaje_nueva
            #Guarda info
            soluciones_encontradas = append(soluciones_encontradas,puntaje_nueva)
            largo = length(solucion_actual$archivos)
            cadenas_encontradas = append(cadenas_encontradas,largo)
          }
          else{
            #Guarda info
            soluciones_encontradas = append(soluciones_encontradas,puntaje_actual)
            largo = length(solucion_actual$archivos)
            cadenas_encontradas = append(cadenas_encontradas,largo)
          }
        }
      }
      gc()
    }
    gc()
  }
  endTime = Sys.time()
  tiempo = endTime - startTime
  todas_soluciones_finales = append(todas_soluciones_finales,soluciones_encontradas[length(soluciones_encontradas)])
  todas_cadenas_finales = append(todas_cadenas_finales,cadenas_encontradas[length(cadenas_encontradas)])
  
#}

#soluciones_sa = rbind(todas_soluciones_finales,todas_cadenas_finales)
resultados =c()
ciclos=c(1:20)

for (i in 1:length(cadenas_encontradas)) {
  vector_prueba=rbind(c(soluciones_encontradas[i],norm(cadenas_encontradas[i])))
  resultado=computeHV(t(vector_prueba), ref.point = c(2,2))
  resultados=append(resultados,resultado)
}

#Gráfica soluciones
plot(ciclos,resultados,type="p",col="red",pch=16,xlab="Iteraciones",ylab="Hipervolumen",main="Gráfica")
lines(ciclos,resultados,type="l",col="red")

#Hipervolumenes SA
hvs_sa = c()
for (i in 1:length(soluciones_encontradas)) {
  v_hv = rbind(c(soluciones_encontradas[i],norm(cadenas_encontradas[i])))
  hv_actual = computeHV(t(v_hv), ref.point = c(2,2))
  hvs_sa = c(hvs_sa,hv_actual)
}

####################### FIN SA ##########################


###################### GENETICO #########################

#AGREGAR CICLO FOR PARA GENERACIONES!

#Funcion que hace mutacion de genes
mutar = function(solucion){
  genes = length(solucion[-1:-2])
  cadena = solucion[-1:-2]
  genes_a_mutar = round(genes/2)
  posiciones = c(1:genes)
  posiciones_usadas = c()
  
  archivos_diff = setdiff(archivos,cadena)
  posicion_archivo = sample(1:(length(archivos_diff)),1)
  
  for (i in 1:genes_a_mutar){
    posiciones_diff = setdiff(posiciones,posiciones_usadas)
    pos_random = sample(posiciones_diff,1)
    archivos_diff = setdiff(archivos,cadena)
    cadena[pos_random] = archivos_diff[posicion_archivo]
    posiciones_usadas = c(posiciones_usadas,pos_random)
  }
  
  sol = obtener_matriz(cadena)
  puntaje = puntaje_sol(sol)
  n_grupo = length(cadena)
  
  vector_solucion = list(c(puntaje,n_grupo,cadena))
  return(vector_solucion)
}

#Funcion que hace el crossover de genes
crossover = function(s1,s2){
  largo_s1 = length(s1[-1:-2])
  largo_s2 = length(s2[-1:-2])
  if(largo_s1 > largo_s2){
    genes = length(s2[-1:-2])
    genes_a_cruzar = round(genes/2)
  }
  else{
    genes = length(s1[-1:-2])
    genes_a_cruzar = round(genes/2)
  }
  
  cadena_s1 = s1[-1:-2]
  cadena_s2 = s2[-1:-2]
  
  for (i in 1:genes_a_cruzar){
    cadena_s1[i] = cadena_s2[i]
  }
  
  sol = obtener_matriz(cadena_s1)
  puntaje = puntaje_sol(sol)
  n_grupo = length(cadena_s1)
  
  vector_solucion = list(c(puntaje,n_grupo,cadena_s1))
  return(vector_solucion)
}

#Funcion que cambia random los genes
random_gen = function(cadena,posicion){
  archivos_diff = setdiff(archivos,cadena)
  posicion = sample(1:(length(archivos_diff)),1)
  cadena_actual[posicion_temporal] = archivos_diff[posicion]
}

prueba_soluciones = list()
#Crear conjunto inicial de soluciones P random
generaciones = 4
tam_P = 10
soluciones = list()
soluciones_r = list()
soluciones_p = list()
soluciones_q = list()
modificacion_genetico = 0.8
sol = sol_random(max_genes,cantidad_genes,archivos)
for (i in 1:tam_P) {
  sol = mod_solucion(sol,modificacion_genetico,max_genes)
  puntaje = puntaje_sol(sol$matriz)
  n_grupo = length(sol$archivos)
  vector_solucion = list(c(puntaje,n_grupo,sol$archivos))
  soluciones = append(soluciones,vector_solucion)
}

for(j in 1:generaciones){
  
  soluciones_p = soluciones
  soluciones_q = soluciones
  
  #Crear conjunto de soluciones Q a partir de P
  #soluciones_q = list()
  #40% mutacion
  cantidad_mutar = round(length(soluciones) * 0.4)
  for (i in 1:cantidad_mutar){
    sol_temporal = soluciones[[i]]
    soluciones_q[[i]] = mutar(sol_temporal)[[1]]
  }
  
  #40% crossover
  cantidad_cruzar = round(length(soluciones) * 0.4)
  contador = cantidad_cruzar + 1
  for (i in 1:round(cantidad_cruzar/2)){
    s1 = soluciones[[contador]]
    s2 = soluciones[[contador+1]]
    crossover_1 = crossover(s1,s2)
    crossover_2 = crossover(s2,s1)
    soluciones_q = reemplazar_solucion(crossover_1,soluciones_q,contador)
    soluciones_q = reemplazar_solucion(crossover_2,soluciones_q,contador+1)
    contador = contador + 2
  }
  
  #20% random
  cantidad_random = round(length(soluciones) * 0.2)
  contador = cantidad_mutar + cantidad_cruzar + 1
  for (i in contador:length(soluciones_p)) {
    s_random = sol_random(max_genes,cantidad_genes,archivos)
    puntaje = puntaje_sol(s_random$matriz)
    n_grupo = length(s_random$archivos)
    vector_solucion = list(c(puntaje,n_grupo,s_random$archivos))
    soluciones_q = reemplazar_solucion(vector_solucion,soluciones_q,contador)
  }
  
  soluciones_r = append(soluciones_p,soluciones_q)
  
  
  #Se escojen las mejores soluciones y se recorta el tamaño
  #En caso de cortar un conjunto se realiza el survival selection
  #Se deberían tener las soluciones
  
  soluciones_ordenamiento = c()
  for (i in 1:length(soluciones_r)) {
    v = c(as.double(soluciones_r[[i]][1]),as.integer(soluciones_r[[i]][2]))
    soluciones_ordenamiento = rbind(soluciones_ordenamiento,v)
  }
  
  ordenamiento_nodominado=doNondominatedSorting(t(soluciones_ordenamiento))
  
  #Tomar los mejores rankings pertenecientes a F y escoger los mejores con
  #Crowding distance (ARREGLAR ESTA PARTE YA QUE NO SE ESTÁ VIENDO BIEN)
  cont = 0
  rank = 1
  soluciones = list()
  while(cont < tam_P){
    indices = which(ordenamiento_nodominado$ranks %in% rank)
    if(length(indices) != 0){
      for (i in indices) {
        if(cont < tam_P){
          soluciones = c(soluciones,soluciones_r[i])
          cont = cont + 1
        }
      }
    }
    rank = rank + 1
  }
  
  prueba_soluciones = append(prueba_soluciones,soluciones)
}

#Teniendo el tamaño P se guarda este conjunto de soluciones para una nueva
#generacion

#####################################





ss = rbind(c(0,1),c(1,0),c(0.5,0.5),c(0.7,0.6),c(0.3,0.5))

ordenamiento_nodominado=doNondominatedSorting(t(ss))

#Ver donde corta
corta_en = 1
soluciones_rank1=soluciones_ordenamiento[which((ordenamiento_nodominado$ranks==corta_en)),]

distancia_aglomeracion=computeCrowdingDistance(t(soluciones_rank1))

#distancia_aglomeracion=computeCrowdingDistance(t(soluciones))


#Contribución por punto
vector_prueba=rbind(c(0.5,0.5),c(0.0,0.0),c(0.8,0.3))
resultado=computeHVContr(t(vector_prueba), ref.point = c(2,2))
print(resultado)

#Hypervolumen para un punto
vector_prueba=rbind(c(0.688,0.4))
resultado=computeHV(t(vector_prueba), ref.point = c(2,2))
print(resultado)

#Hypervolumen para una frontera
vector_prueba=rbind(c(0.714,0.23),c(0.765,0.16),c(0.743,0.20),c(0.714,0.23),c(0.714,0.23))
resultado=computeHV(t(vector_prueba), ref.point = c(2,2))
print(resultado)



vector_prueba = rbind(c(0.7490,0.27),c(0.7470,0.27),c(0.7272,0.28),c(0.7221,0.30),c(0.7112,0.33),c(0.6854,0.40),c(0.6805,0.42))

eje_x = c(0.7490,0.7470,0.7272,0.7221,0.7112,0.6854,0.6805)
eje_y = c(0.27,0.27,0.28,0.30,0.33,0.40,0.42)
plot(eje_x, eje_y, type="o", col="black", ylab="Largo", xlab="Similitud")
