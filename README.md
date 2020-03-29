# SERIM 

Fuentes:

https://web.stanford.edu/~jhj1/teachingdocs/Jones-on-R0.pdf
http://indico.ictp.it/event/7960/session/3/contribution/19/material/slides/0.pdf
https://en.wikipedia.org/wiki/Compartmental_models_in_epidemiology

# Resumen:
 Hice un modelito SEIR , compartimente por simplicidad en 3 grupos etarios: 
   
 -  [0] = 0-19 
 -  [1] = 20-49 
 -  [2] = 50+ 
  
 (el modelo tiene mas modificaciones ver comments en el codigo)
 
En el sistema de ecuaciones pienso que:
  
  - la poblacion 0 y 1 estan en contacto y no con los viejos
  - la poblacion 0 no muere
  - la poblacion 1 esta en contacto con todos
  - la poblacion 2 esta en contacto solo con la 1
 
detalles: 
 - de los contagiados asintomaticos que pasan a recuperarse solos
 - los contagiados, contagian 1/4 del rate que los sanos
 - los contagiados viejos no se sanan, (si o si se infectan)

parametros:
  
  t = entra la secuencia de tiempo en dias, en este caso no es usado, pero uno puede cambiar 
       cosas en funcion del tiempo en que se llama a la funcion. 
       Idea original: cambiar la porpagacion o las muertes segun las camas o aislamiento social en el tiempo
  
  x = vector de poblaciones,ver notacion mia:
  
      S: Sanos. es decir poblacion original. a t-> inf esto tiene que ser cero
      E: Contagiados. poblacion que tiene el virus asintomaticamente
      I: Infectados. Gente que se le detecto el virus
      R: recuperados. Gente que sobrevino a la epidemia
      M: Muertos. gente que NO sobrevino a la epidemia

  params = parametros generales del modelo:
  
      bet: 	beta en wikipedia, multiplicado por los contagiados es el rate de contagio, unidades 1/dias. Paso de S -> E
      a: 	1/tiempo de incuvacion, esto es propio de la sepa de virus (y del cuerpo de la persona).
          yo lo fije 1.0/5.2 dias (Whutan), como propio del virus. Paso de E -> I
      gam:    1/tiempo de reuperacion (con esto tmb calculamos el R0), esto si varia en el tiempo y grupo poblacional.
          En mi "modelo", permito que tanto los contagiados como los infectados se recuperen. Pasa de E -> R y I-> R
      mu:	1/Tiempo de muertos. Paso de I-> M
      N: 	Numero total de la poblacion. 3.5E6 es como la pobacion de la provincia de cordoba
  
  ** Ejemplo de modo de uso **

R> ddo2()
