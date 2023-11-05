###### TRABAJO R - MEGAN TRAYLING JIMENEZ

library(RCurl)
#CARGAR LOS DATOS Y EXAMINARLOS
data <- read.table ("datos-trabajoR.txt", header = TRUE, sep = "\t") #Variable data creada con los datos del archivo 
head(data)#Usado para ver las primeras filas y columnas de la tabla
summary(data)#Para ver un resumen general sobre las variables del data frame
dim(data)#Para ver las dimensiones de la tabla
str(data)#Muestra todos los detalles de los objetos en memoria
#¿Cuántas variables hay? 2
#¿Cuántos tratamiento hay? 5

#Separación datos por tratamiento y variable creando una variable difetente para cada uno
Tratamiento <- data[,1] #Le pedimos que solo me coja la primera columna, lo mismo con las otras variables, pero diciendo que me coja la 2 y 3 columna respectivamente
Variable1 <- data[,2]
Variable2 <- data[,3]
#Para comprobar si me ha creado las variables con los datos correctos, antes de continuar
Tratamiento
Variable1
Variable2

#Creación de BOXPLOT para cada variable y coloreandolos de colores diferentes
boxplot(Variable1~Tratamiento, col=c("light pink"))
boxplot(Variable2~Tratamiento, col=c("light green"))

#GRAFICO DE DISPERSION de las dos variables, cada tratamiento tiene un color distinto 
plot(x = Variable1, y = Variable2, col=Tratamiento)
legend(x = "bottomright", legend = c("Tratamiento 1", "Tratamiento 2", "Tratamiento 3", "Tratamiento 4", "Tratamiento 5"), #Comando usado para que aparezca una leyenda en la parte inferior derecha con los nombres indicados de cada tratamiento
	fill = c("black", "red", "green", "cyan", "blue"), title = "Tratamientos")#Comando para indicar los colores de cada tratamiento en la leyenda y el nombre que le queremos poner

#HISTOGRAMA de cada variable manteniendo los colores anteriores  
hist(Variable1, col=c("light pink"))
hist(Variable2, col=c("light green"))

#FACTOR de la columna de tratamiento guardandolo en una variable 
FactorTratamiento<-factor(Tratamiento)
FactorTratamiento #Para comprobar si me ha creado la variable con los datos que quiero

#MEDIA y DESVIACIÓN ESTANDAR de cada tratamiento usando aggregate y el factor anterior, he creado una variable de cada resultado
resultados_aggregateVar1 <- aggregate(Variable1 ~ FactorTratamiento, data = data, 
	FUN = function(x) c(media = mean(x), desviacion_estandar = sd(x)))
resultados_aggregateVar1 #Para poder comprobar si me ha creado bien la variable 

resultados_aggregateVar2 <- aggregate(Variable2 ~ FactorTratamiento, data = data, 
	FUN = function(x) c(media = mean(x), desviacion_estandar = sd(x)))
resultados_aggregateVar2

#CUANTOS ELEMENTOS TIENE CADA TRATAMIENTO 
table(FactorTratamiento)

#EXTRAER DATOS de Tratamiento1 y Tratamiento4 guardados en una variable diferente
Tratamiento1<-data[1:10,1:3] #He pedido que me guarde en una variable de la fila 1 a la 10(que corresponden con el tratamiento 1) de las columnas 1, 2 y 3
Tratamiento1 #Para comprobar que ha cogido el tratamiento que quiero


Tratamiento4<-data[31:40,1:3] #He pedido que me guarde en una variable de la fila 31 a la 40 (que corresponde con el tratamiento 4) de las columnas 1, 2 y 3
Tratamiento4

#DISTRIBUCIÓN NORMAL usando shapiro test, ya que la muestra no tiene más de 5000 observaciones
Tratamiento1V1<-data[1:10,2] #He creado esta variable para poder calcular la distribución normal del tratamiento que quiero en la variable 1, diciendo que coja solo de la linea 1 a la 10 los datos de la columna 2
Tratamiento4V1<-data[31:40,2] #Esta variable es parecida a la anterior, solo que en este caso le he pedido que me agrupe en una variable los valores de las lineas 31 a 40 de la columna 2

shapiro.test(Tratamiento1V1)  
shapiro.test(Tratamiento4V1)

#De acuerdo a los resultados obtenidos el p-value de Tratamiento 1 (0.06434) y Tratamiento 4 (0.1564) son mayor a 0.05 por lo que la variable1 presenta una distribución normal en ambos valores
#En función del resultado de la prueba de normalidad, como se distribuye de forma normal haria un t-test para COMPRROBAR la HIPOTESIS NULA

t.test(Tratamiento1V1, Tratamiento4V1) 

#En función a los resultados obtenidos, podemos afirmar que hay una diferencia significativa de los datos por lo que la hipotesis nula se rechaza

#Para saber si sus VARIANZAS son iguales hay que hacer Fisher’s F-test

var.test(Tratamiento1V1, Tratamiento4V1)

#Los resultados muestran que hay diferencia significativa, como el resultado es diferente a 1 se rechaza y se puede afirmar que las varianzas NO son iguales 

