#############################################################################
#
# PRACTICA 1 - Megan Trayling Jimenez 130218
#
# Expresión diferencial de genes de ratón
# Microarray de Affymetrix (Affymetrix Murine Genome U74A version 2 MG_U74Av2
# Origen de los datos: GEO GSE5583 (http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE5583)
# Publicación: Mol Cell Biol 2006 Nov;26(21):7913-28.  16940178 (http://www.ncbi.nlm.nih.gov/pubmed/16940178)
#
# Muestras: 3 Wild Type x 3 Histone deacetylase 1 (HDAC1)
#
# R código original (credits): Ahmed Moustafa
#
## ENTREGA EL 01 OCTUBRE 23:59
## Se requiere la entrega de este script completado con los códigos más las imágenes y las respuestas a las preguntas
## Adjuntar en la entrega el PDF final y el archivo con los genes
#
##############################################################################

# Instalar RCurl

if (!requireNamespace("BiocManager"))
    install.packages("BiocManager")
BiocManager::install("RCurl")

# Si esto falla, que seguro lo hace tratar de instalarlo usando el menú, Paquetes, Servidor Spain A Coruña, RCurl

# Cargamos el paquete y los datos
library(RCurl)
url = getURL ("http://bit.ly/GSE5583_data", followlocation = TRUE)
data = as.matrix(read.table (text = url, row.names = 1, header = T))

# Chequeamos las dimensiones de los datos, y vemos las primeras y las últimas filas
dim(data) #Para mirar dimensiones de la tabla
head(data) #Para mirar las primeras filas
tail(data) #Para mirar las últimas filas


# Hacemos un primer histograma para explorar los datos
hist(data) #Para crear histograma

# Transformamos los datos con un logaritmo  -> solo para crear imágenes
data_log=log2(data) #Creamos una transformación logarítmica y la guardamos en otra variable
hist(data_log) #Crear histograma de la transformación logarítmica 

# ¿Qué pasa si hacemos una transformación logarítima de los datos? ¿Para qué sirve? Pasa de tener una distribucion sesgada a mas normal, como una campana de Gauss (mas bonito)

# Hacemos un boxplot con los datos transformados. ¿Qué significan los parámetros que hemos empleado?
boxplot(data_log) #para crear un boxplot de nuesta variable de logaritmicos

boxplot(data_log, col=c("blue","blue","blue","orange","orange","orange"), main="GSE5583-boxplots", las=2) #Para poder acceder a cada elemento, y cambiar el color de cada uno, los tres primero wildtypes de azul y los 3 siguientes knockout de naranja. Y luego cambiarle el nombre

#col=c es para poder cambiar el color de cada elemento por separado
#main es para cambiar el título
#las=2 es para poner los ejes en vertical

# ¿Qué es un boxplot? Diagrama de cajas que representa la mediana, los cuartiles

# Hacemos un hierarchical clustering de las muestras basándonos en un coeficiente de correlación
#Tenemos que tener un cluster con los wildtype y otro con los knockout 

hc=hclust(as.dist(1-cor(data_log))) 
plot(hc, main="GSE5583-Hierarchical Clustering")

# de los valores de expresión. ¿Es correcta la separación? Si, es correcta, en un cluster nos aparecen nos knowckout y en el otro los wildtype

#######################################
# Análisis de Expresión Diferencial 
#######################################

# Primero separamos las dos condiciones (wildtype 1-3, knockout 4-6). ¿Qué tipo de datos has generado? He generado una matiz
wt <- data[,1:3]
ko <- data[,4:6]
class(wt)
head(wt) #Para ver las primeras lineas

#Para acceder a una tabla usar corchetes, la coma es para escoger solo columnas, los dos puntos para que sea del 1 al 3 

# Calcula las medias de las muestras para cada condición. Usa apply -> funcion que calcula donde yo quiera lo que yo quiera (sumas, medias, desviación estandar...)Para ello crear una variable para guardar las medias. En nuestro caso es la media de cada fila de cada condición, por eso ponemos un 1, si fuese un 2 sería para cada columna.
wt.mean=apply(wt, 1, mean)
wt.mean
head(wt.mean)

ko.mean=apply(ko, 1, mean)
ko.mean
head(ko.mean)

# ¿Cuál es la media más alta? Con max podemos saber la media más alta de cada condición. La media más alta de los WT es 35375.53 y la de los KO es 37460.5
max(wt.mean) 
max(ko.mean)

# Ahora hacemos un scatter plot (gráfico de dispersión)-> compara medias del eje X con el eje Y. Vamos a enfre
plot(ko.mean ~ wt.mean)
plot(ko.mean ~ wt.mean, xlab="WT", ylab="KO", main="GSE5583 - Scatter")

#xlab para cambiar nombre eje x y ylab para nombre eje Y

# Añadir una línea diagonal con abline, SIEMPRE CUANDO EL PLOT ESTÉ ABIERTO, la b (de donde tiene que salir la recta) va primero, y la a va segundo (que es el valor de la pendiente, si queremos que vaya por el medio hay que poner 1)
abline(0,1,col="red")
#Para poner una horizontal: 
abline(h=2, col="blue")
#Para poner una vertical: 
abline(v=5, col="green")


# ¿Eres capaz de añadirle un grid?
grid()

# Calculamos la diferencia entre las medias de las condiciones
diff.means = wt.mean - ko.mean

# Hacemos un histograma de las diferencias de medias
hist(diff.means)
hist(diff.means, col="purple")


# Calculamos la significancia estadística con un t-test. Para que sea significativo tiene que ser al menos 0.05
# Primero crea una lista vacía para guardar los p-values (va a haber un p-value para cada gen)
# Segundo crea una lista vacía para guardar las estadísticas del test.
# OJO que aquí usamos los datos SIN TRANSFORMAR. ¿Por qué? Porque al hacerlos con logaritmos los datos que nos da no son fiables, porque están manipulados
# ¿Cuántas valores tiene cada muestra? 12488
pvalue = NULL #para crear la lista 
tstat = NULL
for(i in 1 : nrow(data)){
	x=wt[i,] #creamos variable wt solo cogiendo columnas
	y=ko[i,] #creamos variable ko solo cogiendo columnas

#Hacemos el test
t=t.test(x,y) #varianle t para guardar resultado del test

#Añadimos el p-value a la lista
pvalue[i]=t$p.value
#Añadimos las estadisticas a la linea
tstat[i]=t$statistic
}
head(pvalue) 
length(pvalue)

# Ahora comprobamos que hemos hecho TODOS los cálculos -> salir de la sesión de R y poner el script sin error


# Hacemos un histograma de los p-values. 
hist(pvalue)

# ¿Qué pasa si le ponemos con una transformación de -log10? Lo transformamos porque en el anterior histograma el pvalue estaba en una zona muy pequeña, y así cambia la distribución, así la mayoria está en el 0
hist(-log10(pvalue), col="pink")

# Hacemos un volcano plot. Aquí podemos meter la diferencia de medias y la significancia estadística -> ponemos pvalue en logaritmo porque es lo que queremos representar, nos representa la diferencia de medias en contra del logaritmo
plot(diff.means, -log10(pvalue), main = "GSE5583 - Volcano")
#los valores significativos en nuestro volcano plot están del 4 hacia arriba

# Queremos establecer que el mínimo para considerar una diferencia significativa, es con una diferencia de 2 y un p-value de 0.01
# ¿Puedes representarlo en el gráfico?
diff.means_cutoff = 2 #diferencia de medias en 2 por arriba y en -2 por abajo
pvalue_cutoff = 0.01 #cutoff del pvalue
abline(v = diff.means_cutoff, col = "blue", lwd = 3) #nos marca una recta en el volcano plot, para saber a partir de donde hay significancia, lwd es el ancho de linea
abline(h = -log10(pvalue_cutoff), col = "green", lwd = 3) #ponemos el logaritmo porque sino no nos lo reconoce
#entre ambas franjas se encuentran los valores significativos

# Ahora buscamos los genes que satisfagan estos criterios
# Primero hacemos el filtro para la diferencia de medias (fold)- creamos variables para cada filtro y las combinamos, la diferencia de media es 2 y -2, abs es poner todo en valor absoluto
filter_by_diff.means = abs(diff.means) >= diff.means_cutoff #guardar todos los genes con valores mayores o iguales al cutoff
dim(data[filter_by_diff.means, ])

# Ahora el filtro de p-value
filter_by_pvalue = pvalue <= pvalue_cutoff #lo mismo pero con el filtro de menores o iguales al cutoff
dim(data[filter_by_pvalue, ]) #para que me de el numero de genes que sobrepasan el cutoff

# Ahora las combinamos. ¿Cuántos genes cumplen los dos criterios?
filter_combined = filter_by_diff.means & filter_by_pvalue #y solo me quedo con los genes comunes que hayan pasado el filtro
filtered = data[filter_combined,] #extraer los datos 
dim(filtered) #salen los mismos que los del pvalue, porque los demas son más de un filtro
head(filtered)

# Ahora generamos otro volcano plot con los genes seleccionados marcados en rojo
plot(diff.means, -log10(pvalue), main = "GSE5583 - Volcano #2")
points(diff.means[filter_combined], -log10(pvalue[filter_combined]), col="red")

# Ahora vamos a marcar los que estarían sobreexpresados (rojo) y reprimidos (azul). ¿Por qué parece que están al revés? Porque hemos hecho wt-ko, los sobreepresados son los ko que tienen 2 veces la expresion que el wt, y los reprimidos los que tienen la mitad del wt, al hacer la diferencia de medias nos sale negativo, por ello nos sale que los sobreexpresados estan en el lado negativo 
plot(diff.means, -log10(pvalue), main = "GSE5583 - Volcano #3")
points(diff.means[filter_combined & diff.means < 0],
	-log10(pvalue[filter_combined & diff.means < 0]), col = "red")
points(diff.means[filter_combined & diff.means > 0],
	-log10(pvalue[filter_combined & diff.means > 0]), col = "blue")

# Ahora vamos a generar un mapa. Para ello primero tenemos que hacer un cluster de las columnas y los genes 
# ¿Qué es cada parámetro que hemos usado dentro de la función heatmap? (no lo pregunta) pero labRow=FALSE es para quitar el nombre de los genes
# ¿Eres capaz de cambiar los colores del heatmap? Pista: usar el argumento col y hcl.colors
heatmap(filtered) #para crear el heatmap
#Para agruparlo y colocarlo, agrupa los wt por un lado y ko por otro, y entre las dos mitades de arriba y abajo nos distingue entre sobreexpresados y reprimidos
rowv = as.dendrogram(hclust(as.dist(1-cor(t(filtered)))))
colv = as.dendrogram(hclust(as.dist(1-cor(filtered))))
heatmap(filtered, Rowv=rowv, Colv=colv, cexCol=0.7,labRow=FALSE)

# Ahora vamos a crear un heatmap más chulo. Para ello necesitamos dos paquetes: gplots y RcolorBrewer
#if (!requireNamespace("BiocManager"))
#    install.packages("BiocManager")
#BiocManager::install(c("gplots","RColorBrewer"))
install.packages("gplots")		
install.packages("RColorBrewer")	

library(gplots)

# Hacemos nuestro heatmap


# Lo guardamos en un archivo PDF


# Guardamos los genes diferencialmente expresados y filtrados en un fichero
write.table (filtered, "GSE5583_DE.txt", sep = "\t",
	quote=FALSE)
